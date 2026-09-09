//
// Created by Nikhil Chandra Admal on 11/19/24.
//

#ifndef OILAB_LAMMPS_H
#define OILAB_LAMMPS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <sstream>
#ifdef _WIN32
    #include <io.h>
    #include <windows.h>
#else
    #include <unistd.h>
    #include <sys/stat.h>
#endif
#include <Eigen/Eigen>
#include <iomanip>
#ifdef OILAB_HAVE_LAMMPS_LIBRARY
#include "library.h"
#include <vector>
#endif

/*! The separation below which the minimisation treats two atoms as one and deletes one of them.
 *
 *  A mesostate brings atoms of the two grains exactly on top of one another, so this only has to
 *  clear floating-point noise; it is named here, beside the command that acts on it, so that code
 *  measuring how many atoms a state fuses measures the same thing LAMMPS will fuse.  The two
 *  drifting apart would make the count of coincidences in a configuration disagree with the count
 *  of atoms the relaxation removes, for no reason visible in either. */
constexpr double lammpsOverlapCutoff = 1.0e-2;

/*! Half-thickness of the region the boundary energy is summed over, in Angstrom.
 *
 *  The construction puts the boundary at the middle of the box, so the region is centred there
 *  rather than measured in from the free surfaces: an energy defined relative to the box rather
 *  than to the boundary changes when the structure sits a little to one side, and the relaxation
 *  is free to slide it.  Ten Angstrom either side is several times the range over which a
 *  boundary disturbs its neighbourhood, so the excess is wholly inside it.
 *
 *  Twelve, with a four Angstrom gap before the reference slab, was chosen by scanning both
 *  against the one thing that can be measured without knowing the right answer: take a single
 *  relaxed structure, slide the regions along x, and see how much the reported energy moves.  It
 *  should not move at all.  At heightScaling 4 that scan bottoms out at 6e-4 eV over a whole
 *  angstrom of travel, and the choice is insensitive -- anything putting the region edge in the
 *  plateau and the slab short of the second interface does about as well. */
constexpr double gbHalfThickness = 12.0;

/*! Thickness of the slab the cohesive energy per atom is measured from, in Angstrom.
 *
 *  Taken from just outside the boundary region, on both sides, so it samples material that is
 *  neither disturbed by the boundary nor by a free surface, and samples both grains equally.
 *
 *  This reference matters more than its size suggests: it is multiplied by the number of atoms
 *  in the boundary region, so an error of 4e-4 eV/atom in it moves a boundary energy of a few eV
 *  by 0.1 eV.  That is why it is a slab of real thickness placed away from both the boundary and
 *  the surfaces, rather than the first few atoms that happen to lie a fixed distance from a
 *  box face. */
constexpr double bulkSlabThickness = 5.0;

/*! Gap left between the boundary region and the reference slab, in Angstrom.
 *
 *  Zero puts the slab immediately outside the boundary region, which is only sound if the
 *  boundary's disturbance has died out by the edge of that region.  Where it has not, the slab
 *  measures strained material and calls it bulk -- and since the result multiplies the number of
 *  atoms in the boundary region, a small error there moves the answer a long way.  The gap buys
 *  distance from the boundary; what limits it is the second interface, which a periodic cell
 *  puts at the far edge. */
constexpr double bulkSlabGap = 4.0;

/*! Empty space left beyond the crystal along the non-periodic direction, in Angstrom.
 *
 *  Without it the cell is periodic in x and the two grains meet again across the far face, so
 *  every configuration carries a second grain boundary there -- one nobody asked for, which
 *  relaxes along with the real one and whose atoms sit a fixed distance from the reference slab.
 *  Opening a gap replaces it with two free surfaces, which are not boundaries and are further
 *  from anything that is measured.
 *
 *  Five Angstrom is past the potential's cutoff, so the two surfaces do not see each other. */
constexpr double vacuumThickness = 5.0;

/*! Whether the tethered and the free relaxation share one LAMMPS invocation.
 *
 *  They remain two independent runs either way: with this set, the second begins with its own
 *  clear and read_data, so it starts from the configuration the construction produced and not
 *  from where the tether left the atoms.  What is shared is the process, the parse of the
 *  potential file and the neighbour-list build -- together most of what an invocation costs, and
 *  none of it physics.  Setting it false runs them as two invocations, which is the same
 *  calculation the slow way, and is how to check that claim. */
constexpr bool bothRelaxationsInOneInvocation = true;

/*! Whether the energies go through the LAMMPS library rather than the `lmp` executable.
 *
 *  The two agree: over a 143-state sweep the same states survive, with the same atom counts and
 *  the same fused counts, and every energy within 6e-6 eV -- roughly 2e-5 J/m^2, and smaller than
 *  the 1.3e-4 eV that writing the data file at eight digits was costing.  The residue is
 *  summation order: create_atoms and read_data hand LAMMPS the same atoms in a different
 *  internal arrangement.  It reorders states in the sorted tables, but only ones tied to within
 *  5e-6 eV of each other.
 *
 *  What it buys is a quarter of the LAMMPS time -- 574 against 762 thread-seconds on that sweep,
 *  and LAMMPS is over ninety percent of a sweep.  No process is launched, and the atoms are
 *  handed over in memory rather than through a data file one side writes and the other parses.
 *
 *  Holding one instance open across states was measured too, and saves a further 0.5%: not worth
 *  the questions it raises about what the previous state left behind, so each state opens and
 *  closes its own.  False falls back to the executable, which is also what happens when the
 *  library was not found at configure time. */
constexpr bool useLammpsLibrary = true;

std::tuple<Eigen::MatrixXd,
           Eigen::Matrix3d,
           Eigen::Vector3d> read_oILAB_output(const std::string& path)
{
    Eigen::MatrixXd atoms;
    Eigen::Matrix3d box;
    box.setZero();
    Eigen::Vector<double,9> boxVectorized;
    boxVectorized.setZero();
    Eigen::Vector3d origin;
    origin.setZero();

    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error opening file for reading oilab config file: " << path << std::endl;
    }

    std::string line;
    int number_atoms = 0;

    int lineCount= 0;
    int atomCount= 0;
    int boxCount= -1;
    int originCount= -1;

    while (std::getline(file, line)) {
        lineCount++;
        std::stringstream ss(line);
        std::string field;
        std::vector<std::string> fields;

        while (ss >> field) {
            fields.push_back(field);
        }

        if (lineCount == 1) {
            number_atoms = std::stoi(fields[0]);
            atoms.resize(number_atoms,5);
            atoms.setZero();
        }
        else if (lineCount == 2)
        {
            for (size_t i = 0; i < fields.size(); ++i) {
                if (fields[i].find("Lattice=") != std::string::npos) {
                    boxCount++;
                    continue;
                }
                if (boxCount>=0 && boxCount<9) {
                    boxVectorized(boxCount) = std::stod(fields[i]);
                    boxCount++;
                    continue;
                }
                if (fields[i].find("origin=") != std::string::npos)
                {
                    originCount++;
                    continue;
                }
                if (originCount>-1 && originCount<3) {
                    origin(originCount) = std::stod(fields[i]);
                    originCount++;
                    continue;
                }
            }
        }
        else if (lineCount > 2) {
            // Reading atom data (type, x, y, z, radius)
            Eigen::VectorXd atom_data(5);
            int coordCount= 0;
            for (const auto &field : fields) {
                atom_data(coordCount)= std::stod(field);
                coordCount++;
            }
            atoms.row(atomCount)= atom_data;
            atomCount++;
        }
    }
    box= boxVectorized.reshaped(3,3);

    return {atoms, box, origin};

}

/*-----------------------*/
//std::string write_lammps_datafile(const std::string &filename, const std::vector<std::vector<double>> &box, const std::vector<std::vector<double>> &atom_data, int atom_types) {
void write_lammps_datafile(const std::string &filename,
                           const std::vector<std::vector<double>> &box,
                           const Eigen::MatrixXd& atom_data,
                           int atom_types)
{
    std::ofstream file(filename);
    if (!file.is_open())
        std::cerr << "Error opening file for writing lammps configuration file: " << filename << std::endl;

    file << "# LAMMPS data file via write_data\n";
    file << "\n";
    file << atom_data.rows() << " atoms\n";
    file << atom_types << " atom types\n";
    file << "\n";
    file << std::setprecision(8) << box[0][0] << " " << std::setprecision(8) << box[0][1] << " xlo xhi\n";
    file << std::setprecision(8) << box[1][0] << " " << std::setprecision(8) << box[1][1] << " ylo yhi\n";
    file << std::setprecision(8) << box[2][0] << " " << std::setprecision(8) << box[2][1] << " zlo zhi\n";
    file << "0 0 0 xy xz yz\n";
    file << "\n";
    file << "Atoms # atomic\n";
    file << "\n";
    for (size_t i = 0; i < atom_data.rows(); ++i) {
        //file << atom_data(i,0) << " " << atom_data(i,1) << " " << atom_data(i,2) << " " << atom_data(i,3) << " " << atom_data(i,4) << " 0 0 0\n";
        // Seventeen digits, which is what a double survives a round trip through decimal at.
        // Eight was enough to read back a structure but not to reproduce its energy: against the
        // library path, which is handed the doubles themselves, the truncation moved the
        // as-constructed energy by up to 1.3e-4 eV, and cutting it to 1.6e-6 was almost entirely
        // a matter of writing the digits out.
        file << std::setprecision(17) << atom_data.row(i) << " 0 0 0\n";
    }
}

/*-----------------------*/
/*! Writes the LAMMPS input script that evaluates a configuration.
 *
 *  @param minimize when true a conjugate-gradient relaxation is run before the energy is read
 *         out, so the reported energy is that of the relaxed configuration.  The default keeps
 *         the original behaviour -- a static evaluation of the configuration as written.
 */
void write_lammps_input_script(const std::string &filename,
                               const std::string &infile,
                               const std::string &outfile,
                               const std::string &potential_file_path,
                               const std::string &output_dump_file,
                               bool minimize = false,
                               const std::string &minimized_dump_file = "",
                               double tether_half_width = 0.0,
                               double tether_stiffness = 1.0,
                               bool chain_free_minimization = false,
                               const std::string &free_dump_file = "",
                               bool both_relaxations = false) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for writing lammps input script: " << filename << std::endl;
        return;
    }

    file << "# Find energy of a config and write in " << outfile << " file\n\n";
    file << "variable potential_path string " << potential_file_path << "\n";

    // Everything from the reset to the energy of the as-constructed configuration.  It is a
    // lambda because a run that wants both relaxations emits it twice: `clear` and a fresh
    // read_data put the second relaxation back at the configuration the construction produced,
    // rather than wherever the first one left the atoms.
    const auto emitSetup = [&]()
    {
        file << "clear\n";
        file << "units metal\n";
        // Periodic in the boundary plane, open along the normal.  `m` rather than `f`: the box
        // shrink-wraps to the atoms but never inside the bounds the data file gives, so the
        // vacuum is kept and an atom the deformation pushes outward is followed rather than lost.
        file << "boundary m p p\n";
        file << "atom_style atomic\n";
        file << "neighbor 1.0 bin\n";
        file << "neigh_modify every 1 delay 2 check yes\n";
        file << "read_data " << infile << "\n";
        file << "pair_style      eam/alloy\n";
        // Three species: grain A, grain B, and the atoms the mesostate brings into coincidence.
        // All are the same element, so every type maps to the same entry of the potential; the
        // third exists only so the boundary the construction built can be picked out of the
        // output.  The overlap removal below fuses each coincident pair, and the survivor keeps
        // the species.
        file << "pair_coeff      * * ${potential_path} Cu Cu Cu\n";
        file << "delete_atoms overlap " << lammpsOverlapCutoff << " all all\n";
        file << "variable area equal ly*lz\n";
        // Everything is measured from the middle of the box, which is where the construction puts
        // the boundary, so that a structure sitting a little to one side is measured the same way.
        file << "variable        xmid equal (xlo+xhi)/2\n";
        file << "variable        xlogb equal ${xmid}-" << std::setprecision(8) << gbHalfThickness << "\n";
        file << "variable        xhigb equal ${xmid}+" << std::setprecision(8) << gbHalfThickness << "\n";
        // The reference slabs sit outside the boundary region, one in each grain: two rather than
        // one so that the cohesive energy is not taken from whichever grain happens to lie on the
        // low-x side, and set back by a gap so that neither samples strained material.
        file << "variable        xlobulklo equal ${xmid}-" << std::setprecision(8)
             << gbHalfThickness+bulkSlabGap+bulkSlabThickness << "\n";
        file << "variable        xhibulklo equal ${xmid}-" << std::setprecision(8)
             << gbHalfThickness+bulkSlabGap << "\n";
        file << "variable        xlobulkhi equal ${xmid}+" << std::setprecision(8)
             << gbHalfThickness+bulkSlabGap << "\n";
        file << "variable        xhibulkhi equal ${xmid}+" << std::setprecision(8)
             << gbHalfThickness+bulkSlabGap+bulkSlabThickness << "\n";
        file << "region         GB      block ${xlogb} ${xhigb} INF INF INF INF side in units box\n";
        file << "region         BULKLO  block ${xlobulklo} ${xhibulklo} INF INF INF INF side in units box\n";
        file << "region         BULKHI  block ${xlobulkhi} ${xhibulkhi} INF INF INF INF side in units box\n";
        file << "region         BULK    union 2 BULKLO BULKHI\n";
        file << "group           GB region GB\n";
        file << "group           BULK region BULK\n";
        file << "compute         peratom GB pe/atom\n";
        file << "compute         peratombulk BULK pe/atom\n";
        file << "compute         pe GB reduce sum c_peratom\n";
        // Reduced over BULK, not over GB.  It summed over GB when the reference slab lay inside
        // the boundary region, which worked only because every BULK atom was also a GB atom; the
        // slab now sits outside that region, and summing over GB would collect nothing at all --
        // giving a cohesive energy of zero and a "boundary energy" that is just the raw
        // potential energy.
        file << "compute         pebulk BULK reduce sum c_peratombulk\n";
        file << "variable        peGB equal c_pe\n";
        file << "variable        peBULK equal c_pebulk\n";
        file << "variable        atomsGB equal count(GB)\n";
        file << "variable        atomsBULK equal count(BULK)\n";
        file << "thermo 1000\n";
        file << "compute 1 all ke/atom\n";
        file << "compute cna all cna/atom 3.08133\n";
        file << "compute csys all centro/atom  fcc\n";
        file << "compute 3 all pe/atom\n";
        file << "compute 4 all stress/atom NULL pair\n";
        file << "timestep        0.001\n";
        file << "thermo_style custom step temp ke pe etotal press pxx pyy pzz pxy pxz pyz ly lx lz xy xz yz c_pe v_atomsGB v_peBULK v_atomsBULK\n";
        file << "dump                    OUT0 all custom 10 " << output_dump_file << " id type x y z fx fy fz c_3 c_1 vx vy vz c_4[1] c_4[2] c_4[3] c_4[4] c_4[5] c_4[6]\n";
        // The energy of the configuration as constructed, before anything is relaxed.  Captured
        // with $(...) so the value is frozen rather than re-evaluated later, and taken before the
        // tether is applied so it describes the state the enumeration produced.
        file << "run                     0\n";
        file << "variable        peGBunmin equal $(c_pe)\n";
        file << "variable        peBULKunmin equal $(c_pebulk)\n";
        file << "variable        atomsGBunmin equal $(count(GB))\n";
        file << "variable        atomsBULKunmin equal $(count(BULK))\n";
        file << "variable        GBeneUnmin equal (${peGBunmin}-(${peBULKunmin}/${atomsBULKunmin})"
                "*${atomsGBunmin})\n";
    };

    // Restrain the atoms near the boundary to the positions the construction gave them, so that
    // the relaxation cannot carry the mesostate away from the state it represents.  spring/self
    // remembers each atom's position at the moment the fix is defined, which is why it is defined
    // after the setup rather than inside it.  The boundary lies at the middle of the cell.
    const auto emitTether = [&]()
    {
        file << "variable        xlotether equal (xlo+xhi)/2-" << std::setprecision(8)
             << tether_half_width << "\n";
        file << "variable        xhitether equal (xlo+xhi)/2+" << std::setprecision(8)
             << tether_half_width << "\n";
        file << "region          TETHERREG block ${xlotether} ${xhitether} INF INF INF INF"
                " side in units box\n";
        file << "group           TETHERGRP region TETHERREG\n";
        file << "fix             tether TETHERGRP spring/self " << std::setprecision(8)
             << tether_stiffness << "\n";
    };
    const auto emitMinimize = [&]()
    {
        // The groups are fixed at the moment they are defined, so relaxing here does not change
        // which atoms the GB and BULK sums run over -- only where those atoms sit.
        file << "min_style       cg\n";
        file << "minimize        1e-12 1e-12 100000 100000\n";
        file << "run                     0\n";
    };
    // One line of results, in the layout read_python_outfile() parses.  `mode` is "file" for the
    // first line written and "append" for any that follow, so a run reporting two relaxations
    // leaves two lines and the caller reads one row per relaxation.
    const auto emitReport = [&](const char* mode, const char* spring)
    {
        file << "variable        coh equal (${peBULK}/${atomsBULK})\n";
        file << "variable        GBene equal (${peGB}-${coh}*${atomsGB})\n";
        file << "print \"coh = ${coh} energy = ${peGB} numAtoms = ${atomsGB} GBene = ${GBene}"
                " springE = " << spring << " GBeneU = ${GBeneUnmin} area = ${area}\" "
             << mode << " " << outfile << "\n";
    };

    if (both_relaxations && tether_half_width > 0.0 && minimize)
    {
        // Both relaxations in one invocation, but not one after the other: the second begins with
        // its own clear and read_data, so it starts from the configuration the construction
        // produced rather than from wherever the tether left the atoms.  The physics is exactly
        // what two separate invocations give -- what is saved is a process launch, a potential
        // parse and a neighbour-list build, which together are most of what an invocation costs.
        emitSetup();
        emitTether();
        emitMinimize();
        file << "variable        springEnergy equal $(f_tether)\n";
        if (!minimized_dump_file.empty())
            file << "write_dump all custom " << minimized_dump_file << " id type x y z\n";
        emitReport("file", "${springEnergy}");

        emitSetup();
        emitMinimize();
        if (!free_dump_file.empty())
            file << "write_dump all custom " << free_dump_file << " id type x y z\n";
        emitReport("append", "0.0");
        file << "\n";
        return;
    }

    emitSetup();
    if (tether_half_width > 0.0) emitTether();
    if (minimize) emitMinimize();
    else          file << "run                     0\n";

    if (chain_free_minimization && tether_half_width > 0.0 && minimize)
    {
        // The tethered figures are read off, the restraint is released, and the minimiser is run
        // again from where the tethered structure stands.  That is a different calculation, not a
        // faster spelling of the same one: the free minimum reached from the tethered structure
        // need not be the one reached from the as-constructed structure.  Measured at 15% of
        // states differing by up to 1.6 eV, for a 3% saving, so it is off by default.
        file << "variable        springEnergy equal $(f_tether)\n";
        file << "variable        peGBteth equal $(c_pe)\n";
        file << "variable        peBULKteth equal $(c_pebulk)\n";
        file << "variable        atomsGBteth equal $(count(GB))\n";
        file << "variable        atomsBULKteth equal $(count(BULK))\n";
        file << "variable        coh equal (${peBULKteth}/${atomsBULKteth})\n";
        file << "variable        GBene equal (${peGBteth}-${coh}*${atomsGBteth})\n";
        if (!minimized_dump_file.empty())
            file << "write_dump all custom " << minimized_dump_file << " id type x y z\n";
        file << "unfix           tether\n";
        file << "minimize        1e-12 1e-12 100000 100000\n";
        file << "run                     0\n";
        file << "variable        peGBfree equal $(c_pe)\n";
        file << "variable        peBULKfree equal $(c_pebulk)\n";
        file << "variable        atomsGBfree equal $(count(GB))\n";
        file << "variable        atomsBULKfree equal $(count(BULK))\n";
        file << "variable        GBeneF equal (${peGBfree}-(${peBULKfree}/${atomsBULKfree})"
                "*${atomsGBfree})\n";
        if (!free_dump_file.empty())
            file << "write_dump all custom " << free_dump_file << " id type x y z\n";
        file << "print \"coh = ${coh} energy = ${peGBteth} numAtoms = ${atomsGBteth} "
                "GBene = ${GBene} springE = ${springEnergy} GBeneU = ${GBeneUnmin} "
                "GBeneF = ${GBeneF} area = ${area}\" file " << outfile << "\n\n";
        return;
    }

    // How hard the tether had to work.  A large value says the mesostate is not a minimum of the
    // potential on its own.  It is reported but kept out of GBene, which sums pe/atom and so
    // carries no fix contribution -- a restraint is a constraint, not a physical term.
    if (tether_half_width > 0.0)
        file << "variable        springEnergy equal $(f_tether)\n";
    else
        file << "variable        springEnergy equal 0.0\n";
    if (!minimized_dump_file.empty())
        file << "write_dump all custom " << minimized_dump_file << " id type x y z\n";
    emitReport("file", "${springEnergy}");
    file << "\n";
}

/*-----------------------*/
std::vector<std::vector<double>> read_python_outfile(const std::string &path) {
    std::vector<std::vector<double>> data;
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error opening file for reading lammps output: " << path << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string field;
        double state_id, area, total_energy, gb_energy, gb_density;
        double spring_energy= 0.0, unminimized_energy= 0.0;
        // Present only when both relaxations were run in one invocation; a run that did one
        // leaves it absent and the caller sees the zero.
        double free_energy= 0.0;

        while (ss >> field)
        {
            if (field == "coh") {
                ss >> field;
                ss >> state_id;
            }
            else if (field == "energy") {
                ss >> field;
                ss >> total_energy;
            }
            else if (field == "GBene") {
                ss >> field;
                ss >> gb_energy;
            }
            else if (field == "numAtoms") {
                ss >> field;
                ss >> gb_density;
            }
            else if (field == "springE") {
                ss >> field;
                ss >> spring_energy;
            }
            else if (field == "GBeneU") {
                ss >> field;
                ss >> unminimized_energy;
            }
            else if (field == "GBeneF") {
                ss >> field;
                ss >> free_energy;
            }
            else if (field == "area")
            {
                ss >> field;
                ss >> area;
                break;
            }
        }

        //gb_density = gb_density / area;
        data.push_back({state_id, area, gb_energy, gb_density,
                        spring_energy, unminimized_energy, free_energy});
    }

    return data;
}


/*! @param minimize when true the configuration is relaxed in LAMMPS before its energy is read,
 *         so the returned energy is the minimized one.  False (the default) evaluates the
 *         configuration as it stands.
 *  @param minimizedDumpFile if non-empty, the configuration as LAMMPS leaves it -- relaxed, when
 *         \p minimize is set -- is written there as a single dump snapshot.
 *  @param tetherHalfWidth when positive, atoms within this distance of the boundary plane are
 *         restrained to their as-constructed positions by a harmonic spring, so the relaxation
 *         cannot carry the mesostate away from the state it represents.  Zero relaxes freely.
 *  @param tetherStiffness spring constant of that restraint, in eV/Angstrom^2.
 *  @param springEnergy if non-null, receives the energy stored in the restraint -- how hard it
 *         had to work.  Zero when no tether was applied.
 *  @param unminimizedEnergy if non-null, receives the boundary energy of the configuration as
 *         constructed, before any relaxation.
 */
#ifdef OILAB_HAVE_LAMMPS_LIBRARY
/*! What one state's LAMMPS run reports back. */
struct LammpsResult
{
    double density= 0.0, gbEnergy= 0.0, spring= 0.0, unrelaxed= 0.0, freeEnergy= 0.0;
};

/*! Run one state through the LAMMPS library rather than the `lmp` executable.
 *
 *  The calculation is the one write_lammps_input_script() emits: the same commands in the same
 *  order, the same two independent relaxations, the same regions.  What differs is that no
 *  process is spawned, the atoms are handed over in memory instead of through a data file that
 *  one side writes and the other parses, and the results come back as variables instead of a
 *  printed line that has to be parsed.  The per-state dump of forces and per-atom stresses is
 *  not written -- nothing downstream reads it, and it is the one part of the script whose only
 *  product is a file.
 *
 *  \p atoms is one row per atom laid out as id, type, x, y, z -- the layout
 *  write_lammps_datafile() is given -- and \p box the three pairs of bounds.
 */
inline LammpsResult energyThroughLibrary(const Eigen::MatrixXd& atoms,
                                         const std::vector<std::vector<double>>& box,
                                         const std::string& potentialFile,
                                         const bool minimize,
                                         const std::string& tetheredDumpFile,
                                         const double tetherHalfWidth,
                                         const double tetherStiffness,
                                         const std::string& freeDumpFile,
                                         const bool bothRelaxations)
{
    // read_data wraps atoms into the box along the periodic directions; create_atoms drops
    // whatever it is not given ownership of, so the same wrapping is done here.  y and z only:
    // x is open, and an atom the deformation pushed into the vacuum has to be kept, which is
    // what the shrink-wrap flag passed to create_atoms below arranges.
    const int count= static_cast<int>(atoms.rows());
    std::vector<int> id(count), type(count);
    std::vector<double> position(3*count);
    for (int i= 0; i < count; ++i)
    {
        id[i]  = static_cast<int>(atoms(i,0));
        type[i]= static_cast<int>(atoms(i,1));
        for (int k= 0; k < 3; ++k)
        {
            double x= atoms(i,2+k);
            if (k > 0)
            {
                const double length= box[k][1]-box[k][0];
                while (x <  box[k][0]) x+= length;
                while (x >= box[k][1]) x-= length;
            }
            position[3*i+k]= x;
        }
    }

    void* lmp= nullptr;
    {
        // One instance per state, opened and closed here.  Keeping an instance alive across
        // states was measured too: identical energies, but only a further half a percent, which
        // does not pay for having to reason about what the previous state left behind.  The open
        // is serialised because the first one through initialises library-wide state.
        const char* startup[]= {"lmp", "-log", "none", "-screen", "none"};
        #pragma omp critical(oilabLammpsOpen)
        lmp= lammps_open_no_mpi(5, const_cast<char**>(startup), nullptr);
    }
    if (lmp == nullptr) throw std::runtime_error("could not open a LAMMPS instance");

    const auto fail= [&](const std::string& what)
    {
        std::string message= what;
        if (lammps_has_error(lmp))
        {
            char buffer[1024]= {0};
            lammps_get_last_error_message(lmp, buffer, sizeof buffer);
            message+= ": ";
            message+= buffer;
        }
        lammps_close(lmp);
        throw std::runtime_error(message);
    };
    const auto run= [&](const std::string& commands)
    {
        lammps_commands_string(lmp, commands.c_str());
        if (lammps_has_error(lmp)) fail("LAMMPS command failed");
    };
    const auto value= [&](const char* name)
    {
        void* raw= lammps_extract_variable(lmp, name, nullptr);
        if (raw == nullptr) fail(std::string("LAMMPS variable missing: ") + name);
        const double v= *static_cast<double*>(raw);
        lammps_free(raw);
        return v;
    };

    // Everything from the reset to the energy of the as-constructed configuration.  A run that
    // wants both relaxations calls it twice: `clear` and a fresh set of atoms put the second
    // relaxation back at the configuration the construction produced, rather than at wherever
    // the first one left the atoms.
    const auto setup= [&]()
    {
        std::ostringstream o;
        o << std::setprecision(15)
          << "clear\n"
             "units metal\n"
             // Periodic in the boundary plane, open along the normal.  `m` rather than `f`: the
             // box shrink-wraps to the atoms but never inside these bounds, so the vacuum is kept
             // and an atom the deformation pushes outward is followed rather than lost.
             "boundary m p p\n"
             "atom_style atomic\n"
             "neighbor 1.0 bin\n"
             "neigh_modify every 1 delay 2 check yes\n"
             "region cell block " << box[0][0] << " " << box[0][1] << " "
                                  << box[1][0] << " " << box[1][1] << " "
                                  << box[2][0] << " " << box[2][1] << " units box\n"
             "create_box 3 cell\n";
        run(o.str());

        // Shrink-wrap flag set, matching what read_data allows: x is open, and an atom sitting
        // out in the vacuum belongs to the state as much as any other.  Nothing may go missing,
        // so the count is checked rather than trusted -- a silently dropped atom would show up
        // only as an energy slightly and inexplicably off.
        const int created= lammps_create_atoms(lmp, count, id.data(), type.data(),
                                               position.data(), nullptr, nullptr, 1);
        if (created != count)
            fail("LAMMPS took " + std::to_string(created) + " of " + std::to_string(count)
                 + " atoms");

        o.str(""); o.clear();
        o << std::setprecision(15)
          << "pair_style eam/alloy\n"
             // Three species: grain A, grain B, and the atoms the mesostate brings into
             // coincidence.  All are the same element, so every type maps to the same entry of
             // the potential -- which is also where the masses come from, the atoms having been
             // handed over without any.  The third species exists only so the boundary the
             // construction built can be picked out of the output; the overlap removal below
             // fuses each coincident pair and the survivor keeps the species.
             "pair_coeff * * " << potentialFile << " Cu Cu Cu\n"
             "delete_atoms overlap " << lammpsOverlapCutoff << " all all\n"
             // Everything is measured from the middle of the box, which is where the construction
             // puts the boundary, so a structure sitting a little to one side is measured alike.
             "variable xmid equal (xlo+xhi)/2\n"
             "variable xlogb equal ${xmid}-" << gbHalfThickness << "\n"
             "variable xhigb equal ${xmid}+" << gbHalfThickness << "\n"
             // The reference slabs sit outside the boundary region, one in each grain: two rather
             // than one so the cohesive energy is not taken from whichever grain happens to lie
             // on the low-x side, and set back by a gap so neither samples strained material.
             "variable xlobulklo equal ${xmid}-" << gbHalfThickness+bulkSlabGap+bulkSlabThickness << "\n"
             "variable xhibulklo equal ${xmid}-" << gbHalfThickness+bulkSlabGap << "\n"
             "variable xlobulkhi equal ${xmid}+" << gbHalfThickness+bulkSlabGap << "\n"
             "variable xhibulkhi equal ${xmid}+" << gbHalfThickness+bulkSlabGap+bulkSlabThickness << "\n"
             "region GB     block ${xlogb} ${xhigb} INF INF INF INF side in units box\n"
             "region BULKLO block ${xlobulklo} ${xhibulklo} INF INF INF INF side in units box\n"
             "region BULKHI block ${xlobulkhi} ${xhibulkhi} INF INF INF INF side in units box\n"
             "region BULK   union 2 BULKLO BULKHI\n"
             "group GB region GB\n"
             "group BULK region BULK\n"
             "compute peratom GB pe/atom\n"
             "compute peratombulk BULK pe/atom\n"
             "compute pe GB reduce sum c_peratom\n"
             // Reduced over BULK, not over GB: the slab sits outside the boundary region, and
             // summing over GB would collect nothing at all -- a cohesive energy of zero and a
             // "boundary energy" that is just the raw potential energy.
             "compute pebulk BULK reduce sum c_peratombulk\n"
             "variable peGB equal c_pe\n"
             "variable peBULK equal c_pebulk\n"
             "variable atomsGB equal count(GB)\n"
             "variable atomsBULK equal count(BULK)\n"
             // The energy of the configuration as constructed, before anything is relaxed.
             // Captured with $(...) so the value is frozen rather than re-evaluated later, and
             // taken before the tether is applied so it describes the state as enumerated.
             "run 0\n"
             "variable peGBunmin equal $(c_pe)\n"
             "variable peBULKunmin equal $(c_pebulk)\n"
             "variable atomsGBunmin equal $(count(GB))\n"
             "variable atomsBULKunmin equal $(count(BULK))\n"
             "variable GBeneUnmin equal (${peGBunmin}-(${peBULKunmin}/${atomsBULKunmin})"
             "*${atomsGBunmin})\n";
        run(o.str());
    };
    // Restrain the atoms near the boundary to the positions the construction gave them, so the
    // relaxation cannot carry the mesostate away from the state it represents.  spring/self
    // remembers each atom's position at the moment the fix is defined, which is why this comes
    // after the setup rather than inside it.
    const auto tether= [&]()
    {
        std::ostringstream o;
        o << std::setprecision(15)
          << "variable xlotether equal (xlo+xhi)/2-" << tetherHalfWidth << "\n"
             "variable xhitether equal (xlo+xhi)/2+" << tetherHalfWidth << "\n"
             "region TETHERREG block ${xlotether} ${xhitether} INF INF INF INF side in units box\n"
             "group TETHERGRP region TETHERREG\n"
             "fix tether TETHERGRP spring/self " << tetherStiffness << "\n";
        run(o.str());
    };
    // The groups were fixed when they were defined, so relaxing does not change which atoms the
    // GB and BULK sums run over -- only where those atoms sit.
    const auto relax= [&]()
    { run("min_style cg\nminimize 1e-12 1e-12 100000 100000\nrun 0\n"); };
    const auto snapshot= [&](const std::string& path)
    {
        if (!path.empty()) run("write_dump all custom " + path + " id type x y z\n");
    };
    const auto boundaryEnergy= [&]()
    {
        run("variable coh equal (${peBULK}/${atomsBULK})\n"
            "variable GBene equal (${peGB}-${coh}*${atomsGB})\n");
        return value("GBene");
    };
    // How hard the tether had to work.  A large value says the mesostate is not a minimum of the
    // potential on its own.  Reported, but kept out of GBene, which sums pe/atom and so carries
    // no fix contribution -- a restraint is a constraint, not a physical term.
    const auto springWork= [&]()
    { run("variable springEnergy equal $(f_tether)\n"); return value("springEnergy"); };

    LammpsResult result;
    if (bothRelaxations && tetherHalfWidth > 0.0 && minimize)
    {
        // Both relaxations in one instance, but not one after the other: the second begins with
        // its own clear and its own atoms, so it starts from the configuration the construction
        // produced rather than from wherever the tether left things.  The physics is exactly what
        // two separate runs give; what is saved is a potential parse and a neighbour-list build.
        setup();
        tether();
        relax();
        result.spring   = springWork();
        result.unrelaxed= value("GBeneUnmin");
        result.gbEnergy = boundaryEnergy();
        result.density  = value("atomsGB");
        snapshot(tetheredDumpFile);

        setup();
        relax();
        result.freeEnergy= boundaryEnergy();
        snapshot(freeDumpFile);
    }
    else
    {
        setup();
        if (tetherHalfWidth > 0.0) tether();
        if (minimize) relax();
        else          run("run 0\n");
        if (tetherHalfWidth > 0.0) result.spring= springWork();
        result.unrelaxed= value("GBeneUnmin");
        result.gbEnergy = boundaryEnergy();
        result.density  = value("atomsGB");
        snapshot(tetheredDumpFile);
    }
    lammps_close(lmp);
    return result;
}
#endif

std::pair<double, double> energy(const std::string& lammpsLocation,
                                 const std::string& oilabConfigFile,
                                 const std::string& potentialFile,
                                 bool minimize = false,
                                 const std::string& minimizedDumpFile = "",
                                 double tetherHalfWidth = 0.0,
                                 double tetherStiffness = 1.0,
                                 double* springEnergy = nullptr,
                                 double* unminimizedEnergy = nullptr,
                                 bool chainFreeMinimization = false,
                                 const std::string& freeDumpFile = "",
                                 double* freeEnergy = nullptr,
                                 bool bothRelaxations = false)
{
    // Write data
    std::string threadNumber= std::to_string(omp_get_thread_num());
    std::string lammpsInputFile= "in"+  threadNumber +".find_energy";
    std::string lammpsDataFile= "data" + threadNumber + ".lammps_input";
    std::string lammpsDumpFile= "dump" + threadNumber + ".lammpsConfigs";
    std::string outfile = "lmp_mesostate_energies" + threadNumber + ".txt";

    // Read data
    auto [atoms, box, origin] = read_oILAB_output(oilabConfigFile);

    // Find rotation and new box
    Eigen::Matrix3d R= box.transpose();
    for(int i=0; i<3; ++i)
        R.row(i).normalize();
    Eigen::Matrix3d new_box= R*box;

    // New atoms
    Eigen::MatrixXd new_atoms(atoms.rows(),5);
    for (size_t i = 0; i < atoms.rows(); ++i) {
        new_atoms(i,Eigen::seq(2,4))= R*(atoms(i,Eigen::seq(1,3)).transpose());
        /*
        new_atoms[i][2] = new_atoms[i][2] + R[0][0]*atoms[i][1] + R[0][1]*atoms[i][2] + R[0][2]*atoms[i][3];
        new_atoms[i][3] = new_atoms[i][3] + R[1][0]*atoms[i][1] + R[1][1]*atoms[i][2] + R[1][2]*atoms[i][3];
        new_atoms[i][4] = new_atoms[i][4] + R[2][0]*atoms[i][1] + R[2][1]*atoms[i][2] + R[2][2]*atoms[i][3];
         */
        new_atoms(i,0)= i + 1;
        new_atoms(i,1) = atoms(i,0);
    }

    // Determine nbox values
    std::vector<std::vector<double>> nbox(3, std::vector<double>(2, 0.0));
    for (size_t i = 0; i < 3; ++i) {
        if (i == 0) {
            // The crystal keeps its own extent; the vacuum is added outside it, symmetrically,
            // so the middle of the box is still the middle of the crystal and the boundary and
            // reference regions -- which are placed relative to that middle -- do not move.
            nbox[i][0] = -new_box(i,i) / 2 - vacuumThickness;
            nbox[i][1] =  new_box(i,i) / 2 + vacuumThickness;
        } else if (i == 1) {
            //nbox[i][0] = origin(i);
            nbox[i][0] = 0.0;
            nbox[i][1] = new_box(i,i);
        } else {
            nbox[i][0] = 0.0;
            nbox[i][1] = new_box(i,i);
        }
    }

#ifdef OILAB_HAVE_LAMMPS_LIBRARY
    // The library path, when it is turned on and the run is one it covers.  chainFreeMinimization
    // is not: it is off by default, differs in physics rather than in plumbing, and keeping one
    // spelling of it is worth more than making the faster path handle it.
    if (useLammpsLibrary && !chainFreeMinimization)
    {
        const auto result= energyThroughLibrary(new_atoms, nbox, potentialFile, minimize,
                                                minimizedDumpFile, tetherHalfWidth,
                                                tetherStiffness, freeDumpFile, bothRelaxations);
        if (springEnergy)      *springEnergy=      result.spring;
        if (unminimizedEnergy) *unminimizedEnergy= result.unrelaxed;
        // Only a run that did both relaxations has a free energy to report; one that did not
        // leaves the zero, which is what parsing the executable's single line also gives.
        if (freeEnergy)        *freeEnergy=        bothRelaxations ? result.freeEnergy : 0.0;
        return {result.density, result.gbEnergy};
    }
#endif

    // Write files.  The boundary and reference regions are placed relative to the middle of the
    // box by the script itself -- see gbHalfThickness and bulkSlabThickness -- so nothing about
    // where they sit is decided here.
    write_lammps_datafile(lammpsDataFile, nbox, new_atoms, 3);
    write_lammps_input_script(lammpsInputFile, lammpsDataFile, outfile,
                              potentialFile, lammpsDumpFile, minimize, minimizedDumpFile,
                              tetherHalfWidth, tetherStiffness,
                              chainFreeMinimization, freeDumpFile, bothRelaxations);

    // Run the LAMMPS script
    std::string command = lammpsLocation +" -in " + lammpsInputFile + " > /dev/null 2>&1";
    std::system(command.c_str());

    // Read energy
    auto data_energy = read_python_outfile(outfile);

    if (springEnergy)      *springEnergy=      data_energy[0].size()>4 ? data_energy[0][4] : 0.0;
    if (unminimizedEnergy) *unminimizedEnergy= data_energy[0].size()>5 ? data_energy[0][5] : 0.0;
    // Two relaxations in one invocation leave one line each: the first is the tethered run, the
    // second the free one, and the free boundary energy is that row's GBene.  A chained run
    // instead reports both on one line, with the free energy under its own key.
    if (freeEnergy) {
        if (bothRelaxations && data_energy.size()>1) *freeEnergy= data_energy[1][2];
        else *freeEnergy= data_energy[0].size()>6 ? data_energy[0][6] : 0.0;
    }
    return {data_energy[0][3], data_energy[0][2]};
}

#endif //OILAB_LAMMPS_H
