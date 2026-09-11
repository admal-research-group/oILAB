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
#include <atomic>
#ifdef _WIN32
    #include <io.h>
    #include <windows.h>
#else
    #include <unistd.h>
    #include <sys/stat.h>
#endif
#include <Eigen/Eigen>
#include <iomanip>
#include <atomic>
#include <chrono>
#include <vector>
#include <stdexcept>
#include <omp.h>   // the log and the failure message name the thread the instance belongs to
#ifndef OILAB_HAVE_LAMMPS_LIBRARY
#error "oILAB needs the LAMMPS library: install a LAMMPS development package so that CMake \
finds liblammps and library.h.  Energies are computed in this process now -- there is no longer \
a path that spawns an `lmp` executable to fall back to."
#endif
#include "library.h"

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

/*! Whether the tethered and the free relaxation share one LAMMPS setup.
 *
 *  They remain two independent runs either way: with this set, the second begins with its own
 *  fresh set of atoms, so it starts from the configuration the construction produced and not from
 *  where the tether left them.  What is shared is the box, the pair style -- hence the parse of
 *  the potential file -- and the neighbour-list build, together most of what a relaxation costs
 *  outside the minimisation itself, and none of it physics.  Setting it false runs them as two
 *  separate calls, which is the same calculation the slow way, and is how to check that claim. */
constexpr bool bothRelaxationsInOneInvocation = true;

/*! \brief Whether the free (untethered) relaxation is run alongside the tethered one.
 *
 *  The two relaxations are independent runs from the same configuration, so dropping the free
 *  one removes very close to half the LAMMPS work.  A production sweep that only wants the
 *  tethered structure -- the one that still represents the mesostate it was built from -- has
 *  no use for it.
 *
 *  With this off, Relaxations::full is left at zero and the sortedByFull table is not written,
 *  rather than quietly reporting the tethered figure twice.  OILAB_FREE_RELAXATION=0 turns it
 *  off from the environment.
 *
 *  Not const: an application may set it directly, which is how a sweep keeps the choice beside
 *  its other settings rather than in the environment of whoever launches it.  Assign it before
 *  any thread starts -- every LAMMPS call reads it. */
inline bool runFreeRelaxation =
    std::getenv("OILAB_FREE_RELAXATION") == nullptr
    || std::string(std::getenv("OILAB_FREE_RELAXATION")) != "0";

/*! \brief Minimiser iterations, summed over every relaxation, and the number of relaxations.
 *
 *  Here to answer one question: a faceted sweep slows down as it runs, and the phase timers put
 *  95% of a state in LAMMPS, so either the minimiser is doing more work per state as the run
 *  proceeds or something is accumulating across instances.  Those look identical from outside
 *  and completely different from here -- iterations rising says the states are genuinely harder,
 *  iterations flat while the time grows says the cost is not the minimisation at all.
 *
 *  Each relaxation starts from a fresh instance whose timestep is zero, so LAMMPS's `step` after
 *  minimize IS the iteration count.  Costs one variable evaluation per relaxation against a
 *  minimisation of order a hundred milliseconds. */
/*! \brief Whether the per-relaxation diagnostics below are collected.  OFF by default.
 *
 *  They cost a LAMMPS command and three library calls per relaxation -- of order 15-35 us
 *  against a minimisation of about 40 ms, so under 0.1% -- but they sit in the hottest path in
 *  the program and both questions they were written to answer are settled: the tethered
 *  minimisation converges in about nine iterations from the first state to the last, and LAMMPS's
 *  own allocation is constant at 4.236 MB for a whole run.  Neither moves while the sweep decays,
 *  which is how the decay was traced out of LAMMPS altogether and into the allocator.
 *
 *  Turn them back on with OILAB_LAMMPS_STATS=1 if a sweep on another boundary, potential or cell
 *  behaves differently -- iterations rising would mean the states really are getting harder, and
 *  that is worth knowing before anything else is blamed. */
inline const bool lammpsReportStatistics =
    std::getenv("OILAB_LAMMPS_STATS") != nullptr
    && std::string(std::getenv("OILAB_LAMMPS_STATS")) != "0";

/*! LAMMPS's own accounting of what it has allocated, sampled once per state.  The point is to
 *  separate two superimposed decays: if this climbs with session.uses and falls back at every
 *  recycle, the accumulation is inside the instance; if it is flat while the sweep still slows,
 *  the cost is in the host process -- heap, allocator, locality -- and no LAMMPS setting reaches
 *  it.  Kept in kB so the sum stays exact in an integer. */
inline std::atomic<long long> lammpsMemoryKB{0};
inline std::atomic<long long> lammpsMemorySamples{0};
inline std::atomic<long long> lammpsMemoryPeakKB{0};

inline std::atomic<long long> minimizeIterations{0};
inline std::atomic<long long> minimizeRelaxations{0};
inline std::atomic<long long> minimizeCapped{0};   //!< relaxations that hit minimizeMaxIterations

/*! Conjugate gradient, and how hard it is asked to work.
 *
 *  Conjugate gradient rather than fire or hftn: the structures start close to a minimum, which
 *  is the case cg is good at, and it is what every energy quoted from this code so far was
 *  produced with.
 *
 *  The tolerances were 1e-12.  Measured over a 1119-state sweep, 1e-10 costs 16% less wall clock
 *  -- the mean relaxation drops from 144 steps to 113 -- and leaves the 100 lowest-energy states
 *  within 5.8e-5 eV of where 1e-12 put them, which is far below anything physical.  Fifteen of
 *  the 1119 move by more than 1e-4 eV and two by more than 1e-2: those are states whose free
 *  relaxation is soft enough to fall into a neighbouring basin, and no tolerance short of exact
 *  arithmetic settles which one is right.  1e-8 is a different matter and was rejected -- there
 *  72 of the lowest 100 move, one of them by 1.37 eV.
 *
 *  The iteration cap was 100000 and is now 5000.  It binds on nothing: the longest relaxation
 *  measured over that sweep took 756 steps at 1e-12, and fewer at 1e-10.  A cap that cannot be
 *  reached is not a limit on the answer, only on how long a pathological state may run. */
constexpr const char* minimizeStyle = "cg";
constexpr const char* minimizeEnergyTolerance = "1e-10";
constexpr const char* minimizeForceTolerance  = "1e-10";
constexpr int minimizeMaxIterations = 5000;

/*! The energies go through the LAMMPS library, in this process.  There is no executable path.
 *
 *  There was one, and the two agreed: over a 143-state sweep the same states survived, with the
 *  same atom counts and the same fused counts, and every energy within 6e-6 eV -- roughly 2e-5
 *  J/m^2, and smaller than the 1.3e-4 eV that writing the data file at eight digits was costing.
 *  The residue was summation order: create_atoms and read_data hand LAMMPS the same atoms in a
 *  different internal arrangement.  It reordered states in the sorted tables, but only ones tied
 *  to within 5e-6 eV of each other.
 *
 *  What the library buys is a quarter of the LAMMPS time -- 574 against 762 thread-seconds on
 *  that sweep, and LAMMPS is over ninety percent of a sweep.  No process is launched, the atoms
 *  are handed over in memory rather than through a data file one side writes and the other
 *  parses, and the potential is parsed once per thread rather than once per relaxation.
 *
 *  The executable path is gone rather than kept as a fallback: it had to be given the location of
 *  an `lmp` binary, which meant every caller from the application down carried a path to a
 *  program that was never run, and a machine without that binary silently computed no energies at
 *  all even with the library linked in.  A build without the LAMMPS development package now fails
 *  at the #error above instead of at run time. */

inline std::tuple<Eigen::MatrixXd,
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
/*! \brief Whether each thread keeps one LAMMPS instance alive and reuses it across states.
 *
 *  On, and this is what keeps the potential file from being parsed more than once per thread.
 *  Closing the instance destroys the pair style along with everything else, so a fresh instance
 *  per relaxation re-reads the EAM file every time.  That file is 706 kB of ASCII and parsing it
 *  costs 33.8 ms -- measured -- against about 1000 ms of LAMMPS work per state, and without this
 *  it happens twice per state, because the two relaxations are independent.  Opening the instance
 *  is a further 7.7 ms each.
 *
 *  Half of that is recovered whatever this is set to, because the instance survives the two
 *  relaxations OF ONE STATE: the second reuses the box and the pair style instead of clearing
 *  them, so the file is parsed once per state rather than twice.  That much is bit-for-bit
 *  deterministic.
 *
 *  Keeping the instance ACROSS states -- what this flag does -- takes it from once per state to
 *  once per thread, and is worth a further 2%.  It costs a little reproducibility: delete_atoms
 *  and create_atoms leave the surviving atoms in an order that depends on what that thread did
 *  before, and with 100 threads that depends on scheduling, so two runs of the same sweep can
 *  disagree on a couple of states in the last digits of the freely relaxed energy -- measured at
 *  2.2e-7 eV over 1119 states, which is physically nothing.  OILAB_LAMMPS_REUSE=0 turns it off
 *  and buys back exact reproducibility at the price of re-parsing the potential once per state. */
inline const bool reuseLammpsInstance =
    std::getenv("OILAB_LAMMPS_REUSE") == nullptr
    || std::string(std::getenv("OILAB_LAMMPS_REUSE")) != "0";

/*! \brief How many states one instance serves before it is closed and reopened.  0 never recycles.
 *
 *  Reuse is not free in the way it looks.  Measured on the faceted sigma5 (310) sweep, tethered
 *  only, with the per-interval phase timers:
 *
 *      progress     reuse ON     reuse OFF      (LAMMPS core-ms per state)
 *         0%          35.9          88.4
 *        17%          90.7         121.6
 *      slope         3.22/%        1.95/%
 *
 *  The minimiser itself does not change -- about nine iterations a relaxation throughout -- so
 *  none of that rise is physics.  A reused instance grows at two thirds again the rate of a
 *  fresh one, so something accumulates inside it; but a fresh one starts 42 ms per state worse,
 *  because closing the instance destroys the pair style and the 706 kB EAM file is parsed again.
 *  Neither setting is the answer: ON grows faster, OFF is slower everywhere.
 *
 *  So keep the instance and retire it on a schedule.  At 500 states the amortised cost of the
 *  reopen is 42/500 = 0.08 ms a state, which is nothing beside the tens of milliseconds the
 *  accumulation was adding, and the instance never lives long enough to accumulate much.
 *
 *  It does not affect reproducibility beyond what reuse already costs -- see reuseLammpsInstance
 *  on atom ordering.  OILAB_LAMMPS_REUSE_STATES sets it; 0 restores unbounded reuse. */
inline const int lammpsRecycleAfter = []{
    const char* v= std::getenv("OILAB_LAMMPS_REUSE_STATES");
    if (v == nullptr) return 500;
    const int n= std::atoi(v);
    return n >= 0 ? n : 500;
}();

/*! \brief Where LAMMPS writes its log, when it is asked to write one at all.
 *
 *  Empty -- the default -- runs every instance with `-log none -screen none`, which is what a
 *  sweep wants: a hundred threads each narrating their own minimisation is not output anyone
 *  reads, and it was the reason the old executable path was run with its output sent to
 *  /dev/null.
 *
 *  Setting OILAB_LAMMPS_LOG turns the log back on.  The value is a filename prefix, and each
 *  thread gets its own file, `<prefix>.<thread>.log`, since threads share nothing else; the
 *  value `1` means the prefix `lammps`.  Every command this file issues is echoed into it, so
 *  the log is a script that reproduces the run, followed by what LAMMPS made of it.  With
 *  instance reuse on -- the default -- one file holds every state that thread visited, in order,
 *  which is what makes it possible to see what the state before the failing one left behind.
 *
 *  OILAB_LAMMPS_SCREEN sends the same thing to stderr instead of a file.  Only useful with
 *  OILAB_THREADS=1, and then only for watching a single state go by. */
inline const std::string lammpsLogPrefix = []() -> std::string
{
    const char* v= std::getenv("OILAB_LAMMPS_LOG");
    if (v == nullptr || *v == '\0') return "";
    return std::string(v) == "1" ? "lammps" : v;
}();

inline const bool lammpsToScreen = std::getenv("OILAB_LAMMPS_SCREEN") != nullptr;

/*! One thread's LAMMPS instance, with the potential already read.  Closed when the thread ends. */
struct LammpsSession
{
    void* handle= nullptr;
    std::string potential;
    std::vector<double> cell;      //!< the bounds create_box was given, to know when to move them
    long long uses= 0;             //!< states this instance has served, for lammpsRecycleAfter
    bool stateBuilt= false;        //!< whether a previous state's objects are still defined
    bool tetherDefined= false;     //!< whether that state left a spring/self fix behind
    /*! The commands most recently handed to this instance.
     *
     *  LAMMPS reports which command failed but not what it was told to do, and with `-log none`
     *  -- the default -- there is nothing else left behind to look at.  Keeping the last block
     *  means the failure can say it, which is usually enough to place the fault without turning
     *  the log on and running the sweep again. */
    std::string lastCommands;
    ~LammpsSession() { if (handle) lammps_close(handle); }
};

/*! What one state's LAMMPS run reports back. */
struct LammpsResult
{
    double density= 0.0, gbEnergy= 0.0, spring= 0.0, unrelaxed= 0.0, freeEnergy= 0.0;
};

/*! Run one state through the LAMMPS library, in this process.
 *
 *  The calculation is the one the old input script emitted: the same commands in the same order,
 *  the same two independent relaxations, the same regions.  What differs is that no process is
 *  spawned, the atoms are handed over in memory instead of through a data file that one side
 *  writes and the other parses, and the results come back as variables instead of a printed line
 *  that has to be parsed.  The per-state dump of forces and per-atom stresses is not written --
 *  nothing downstream read it, and it was the one part of the script whose only product was a
 *  file.
 *
 *  \p atoms is one row per atom laid out as id, type, x, y, z, and \p box the three pairs of
 *  bounds.
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

    // The instance is this thread's, kept between states so the potential is parsed once -- see
    // reuseLammpsInstance.  The open is serialised because the first one through initialises
    // library-wide state.
    static thread_local LammpsSession session;
    const bool reuse= reuseLammpsInstance;
    // Retired when it has served its quota, as well as when reuse is off or the potential has
    // changed.  Recycling here rather than after the run keeps the decision in one place, and
    // the next block reopens whatever this closed.
    const bool spent= reuse && lammpsRecycleAfter > 0 && session.uses >= lammpsRecycleAfter;
    if (session.handle != nullptr
        && (!reuse || spent || session.potential != potentialFile)) {
        lammps_close(session.handle);
        session= LammpsSession{};
    }
    // Silent unless asked otherwise -- see lammpsLogPrefix.  `-echo log` puts every command into
    // the log alongside what LAMMPS made of it, so a log that was asked for is a script that
    // replays the run rather than a bare record of its output.
    const std::string logFile=
        lammpsLogPrefix.empty()
            ? std::string("none")
            : lammpsLogPrefix + "." + std::to_string(omp_get_thread_num()) + ".log";
    if (session.handle == nullptr) {
        const char* startup[]= {"lmp", "-log", logFile.c_str(),
                                "-screen", lammpsToScreen ? "/dev/stderr" : "none",
                                "-echo", lammpsLogPrefix.empty() ? "none" : "log"};
        #pragma omp critical(oilabLammpsOpen)
        session.handle= lammps_open_no_mpi(7, const_cast<char**>(startup), nullptr);
    }
    void* lmp= session.handle;
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
        // What the run was, so the failure can be placed without turning the log on and repeating
        // the sweep.  The commands come first because they are what usually identifies the fault;
        // the configuration behind them is what says whether the state itself was the problem.
        message+= "\n  thread " + std::to_string(omp_get_thread_num())
                + ",  " + std::to_string(count) + " atom(s)"
                + ",  box x [" + std::to_string(box[0][0]) + ", " + std::to_string(box[0][1])
                + "]  y [" + std::to_string(box[1][0]) + ", " + std::to_string(box[1][1])
                + "]  z [" + std::to_string(box[2][0]) + ", " + std::to_string(box[2][1]) + "]"
                + "\n  potential " + potentialFile
                + "\n  while running:\n" + session.lastCommands;
        if (lammpsLogPrefix.empty())
            message+= "\n  Set OILAB_LAMMPS_LOG=1 (and OILAB_THREADS=1) for the full LAMMPS log "
                      "of a repeat run.";
        else
            message+= "\n  The LAMMPS log of this run is " + logFile + ".";
        // The instance is this thread's and is about to be abandoned; clear the session too, so
        // the destructor does not close it a second time and the next state opens a fresh one.
        lammps_close(lmp);
        session= LammpsSession{};
        throw std::runtime_error(message);
    };
    const auto run= [&](const std::string& commands)
    {
        session.lastCommands= commands;
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
    // Undo what the previous state defined on this instance.  Fixes and computes go before the
    // groups they name, groups before the regions that made them, and the atoms last; nothing
    // survives but the box and the pair style.  Redefining any of these while the old one still
    // exists is an error in LAMMPS, so this is not optional housekeeping.
    const auto teardown= [&]()
    {
        if (session.tetherDefined) { run("unfix tether\ngroup TETHERGRP delete\n"
                                        "region TETHERREG delete\n");
                                     session.tetherDefined= false; }
        run("uncompute pe\nuncompute pebulk\nuncompute peratom\nuncompute peratombulk\n"
            "group GB delete\ngroup BULK delete\n"
            "region GB delete\nregion BULKLO delete\nregion BULKHI delete\n"
            "region BULK delete\n"
            "delete_atoms group all compress no\n");
    };

    const auto setup= [&]()
    {
        std::ostringstream o;
        if (session.handle != nullptr && session.stateBuilt) teardown();
        if (!session.cell.empty())
        {
            // Same instance, and the cell only moves when a state overshoots it; change_box when
            // it has.
            const std::vector<double> wanted{box[0][0],box[0][1],box[1][0],
                                             box[1][1],box[2][0],box[2][1]};
            if (session.cell != wanted) {
                o << std::setprecision(15)
                  << "change_box all x final " << box[0][0] << " " << box[0][1]
                  << " y final " << box[1][0] << " " << box[1][1]
                  << " z final " << box[2][0] << " " << box[2][1] << " units box\n";
                run(o.str());
                o.str(""); o.clear();
                session.cell= wanted;
            }
        }
        else
        {
            o << std::setprecision(15)
              << "units metal\n"
                 // Periodic in the boundary plane, open along the normal.  `m` rather than `f`:
                 // the box shrink-wraps to the atoms but never inside these bounds, so the vacuum
                 // is kept and an atom the deformation pushes outward is followed, not lost.
                 "boundary m p p\n"
                 "atom_style atomic\n"
                 "neighbor 1.0 bin\n"
                 "neigh_modify every 1 delay 2 check yes\n"
                 "region cell block " << box[0][0] << " " << box[0][1] << " "
                                      << box[1][0] << " " << box[1][1] << " "
                                      << box[2][0] << " " << box[2][1] << " units box\n"
                 "create_box 3 cell\n"
                 "pair_style eam/alloy\n"
                 // Three species: grain A, grain B, and the atoms the mesostate brings into
                 // coincidence.  All are the same element, so every type maps to the same entry
                 // of the potential -- which is also where the masses come from, the atoms having
                 // been handed over without any.  The third species exists only so the boundary
                 // can be picked out of the output; the overlap removal fuses each coincident
                 // pair and the survivor keeps the species.
                 "pair_coeff * * " << potentialFile << " Cu Cu Cu\n";
            run(o.str());
            o.str(""); o.clear();
            session.cell= {box[0][0],box[0][1],box[1][0],box[1][1],box[2][0],box[2][1]};
            session.potential= potentialFile;
        }

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
        o.str(""); o.clear();
        o << std::setprecision(15)
          << "delete_atoms overlap " << lammpsOverlapCutoff << " all all\n"
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
        session.stateBuilt= true;
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
        session.tetherDefined= true;
    };
    // The groups were fixed when they were defined, so relaxing does not change which atoms the
    // GB and BULK sums run over -- only where those atoms sit.
    // OILAB_NO_MINIMIZE evaluates the configuration as it stands instead of relaxing it, for
    // separating the cost of the minimisation from everything else LAMMPS does per state.
    static const bool skipMinimize= std::getenv("OILAB_NO_MINIMIZE") != nullptr;
    const auto relax= [&]()
    {
        if (skipMinimize) { run("run 0\n"); return; }
        std::ostringstream m;
        m << "min_style " << minimizeStyle << "\nminimize "
          << minimizeEnergyTolerance << " " << minimizeForceTolerance << " "
          << minimizeMaxIterations << " " << minimizeMaxIterations << "\nrun 0\n";

        // BRACKET the minimisation rather than read `step` after it.  Reading it after only
        // gives the iteration count if the instance started at zero, and it does not: the first
        // attempt at this reported a mean that climbed by a constant ~910 every progress
        // interval, which is a running total being divided by a constant number of relaxations,
        // not a minimiser doing more work.  The difference is right however the timestep is
        // carried, and costs one extra variable evaluation.
        long long before= 0;
        if (lammpsReportStatistics) {
            run("variable oilabStep equal step\n");
            before= static_cast<long long>(value("oilabStep"));
        }

        run(m.str());

        if (lammpsReportStatistics) {
            const long long steps= static_cast<long long>(value("oilabStep")) - before;
            minimizeIterations.fetch_add(steps, std::memory_order_relaxed);
            minimizeRelaxations.fetch_add(1, std::memory_order_relaxed);
            if (steps >= minimizeMaxIterations)
                minimizeCapped.fetch_add(1, std::memory_order_relaxed);
        }
    };
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

        if (runFreeRelaxation) {
            setup();
            relax();
            result.freeEnergy= boundaryEnergy();
            snapshot(freeDumpFile);
        }
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
    // Sampled before any close, so it describes the instance as this state left it.
    if (lammpsReportStatistics)
    {
        double meminfo[3]= {0.0,0.0,0.0};
        lammps_memory_usage(lmp, meminfo);
        const long long kb= static_cast<long long>(meminfo[0]*1024.0);   // meminfo[0] is MB
        lammpsMemoryKB.fetch_add(kb, std::memory_order_relaxed);
        lammpsMemorySamples.fetch_add(1, std::memory_order_relaxed);
        long long peak= lammpsMemoryPeakKB.load(std::memory_order_relaxed);
        while (kb > peak
               && !lammpsMemoryPeakKB.compare_exchange_weak(peak, kb,
                                                            std::memory_order_relaxed)) {}
    }

    ++session.uses;
    if (!reuse) { lammps_close(lmp); session= LammpsSession{}; }
    return result;
}

/*! As energy() below, but taking the configuration already in memory rather than a path.
 *
 *  This is the form the sweep uses.  The atoms were built in this process a moment earlier, so
 *  writing them out and parsing them back was pure overhead once LAMMPS stopped needing a file;
 *  the path-taking overload reads the file and delegates here, for callers that only have a path.
 *
 *  \p atoms is one row per atom -- species, x, y, z, radius -- \p configBox holds the three cell
 *  vectors as columns and \p configOrigin the cell origin, which is what read_oILAB_output()
 *  returns and what GbMesoState::box() now fills in directly.
 */
inline std::pair<double, double> energy(const Eigen::MatrixXd& atoms,
                                 const Eigen::Matrix3d& configBox,
                                 const Eigen::Vector3d& configOrigin,
                                 const std::string& potentialFile,
                                 bool minimize = false,
                                 const std::string& minimizedDumpFile = "",
                                 double tetherHalfWidth = 0.0,
                                 double tetherStiffness = 1.0,
                                 double* springEnergy = nullptr,
                                 double* unminimizedEnergy = nullptr,
                                 const std::string& freeDumpFile = "",
                                 double* freeEnergy = nullptr,
                                 bool bothRelaxations = false)
{
    const Eigen::Matrix3d& box= configBox;
    const Eigen::Vector3d& origin= configOrigin;
    (void)origin;

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

    // The boundary and reference regions are placed relative to the middle of the box by the
    // commands energyThroughLibrary() issues -- see gbHalfThickness and bulkSlabThickness -- so
    // nothing about where they sit is decided here.
    const auto result= energyThroughLibrary(new_atoms, nbox, potentialFile, minimize,
                                            minimizedDumpFile, tetherHalfWidth,
                                            tetherStiffness, freeDumpFile, bothRelaxations);
    if (springEnergy)      *springEnergy=      result.spring;
    if (unminimizedEnergy) *unminimizedEnergy= result.unrelaxed;
    // Only a run that did both relaxations has a free energy to report; one that did not leaves
    // the zero.
    if (freeEnergy)        *freeEnergy=        bothRelaxations ? result.freeEnergy : 0.0;
    return {result.density, result.gbEnergy};
}

/*! Reads an oILAB configuration file and hands it to the overload above.
 *
 *  Kept for callers that have a path rather than the atoms.  The sweep does not go through here
 *  any more: it passes the configuration GbMesoState::box() built, which saves writing the file
 *  and parsing it back. */
inline std::pair<double, double> energy(const std::string& oilabConfigFile,
                                        const std::string& potentialFile,
                                        bool minimize = false,
                                        const std::string& minimizedDumpFile = "",
                                        double tetherHalfWidth = 0.0,
                                        double tetherStiffness = 1.0,
                                        double* springEnergy = nullptr,
                                        double* unminimizedEnergy = nullptr,
                                        const std::string& freeDumpFile = "",
                                        double* freeEnergy = nullptr,
                                        bool bothRelaxations = false)
{
    const auto [atoms, box, origin] = read_oILAB_output(oilabConfigFile);
    return energy(atoms, box, origin, potentialFile, minimize, minimizedDumpFile,
                  tetherHalfWidth, tetherStiffness, springEnergy, unminimizedEnergy,
                  freeDumpFile, freeEnergy, bothRelaxations);
}

#endif //OILAB_LAMMPS_H
