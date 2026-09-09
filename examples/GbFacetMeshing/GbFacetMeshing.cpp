//
// Created by himanshu on 9/1/26.
//

#include <cassert>
#include <TextFileParser.h>
#include <GbMesoStateEnsemble.h>
#include <omp.h>
#include <numbers>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <tuple>
#include <map>
#include <iomanip>
#include <filesystem>
#include <cfloat>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cstring>
#include <chrono>

// Where sweep output is written.  CMake supplies <project>/runs; the fallback only matters if
// the example is compiled outside it.
#ifndef GBFACETMESHING_RUNS_DIR
#define GBFACETMESHING_RUNS_DIR "runs"
#endif

std::vector<std::vector<int>> readPostProcessingOutput(const std::string& filename)
{
	std::ifstream inputFile(filename);
	if (!inputFile.is_open())
        	throw std::runtime_error("Error: could not open file " + filename);
	std::vector<std::vector<int>> data;
	std::string line;

	while(std::getline(inputFile,line))
	{
		if (line.empty()) continue;

		std::istringstream iss(line);
        	std::vector<int> row;
        	int value;
        	while (iss >> value)
            		row.push_back(value);
        	if (!row.empty())
            		data.push_back(std::move(row));
    	}

    	return data;
}

struct LAMMPSData
{
	int natoms = 0;
	int ntypes = 0;
	std::vector<std::vector<double>> box;
	std::vector<std::vector<double>> atoms;
};

LAMMPSData readXYZFile(const std::string&  filename)
{
	LAMMPSData results;
	results.box.resize(3, std::vector<double>(2, 0.0));
        int b_count = 0;
	int a_count = 0;

        std::ifstream inputFile(filename);
        if (!inputFile.is_open())
                throw std::runtime_error("Error: could not open file " + filename);

        std::string line;

        while(std::getline(inputFile,line))
        {
                if(line.empty()) continue;

		std::istringstream iss(line);
                std::string word;
                std::vector<std::string> fields;
                while(iss >> word) fields.push_back(word);

		if (fields.size()==1)
		{
			results.natoms = std::stoi(fields[0]);
			continue;
		}

		if (fields.size()==5)
		{
			std::vector<double> row;
			a_count += 1;
			row.push_back(a_count);
                        for (size_t i = 0; i < fields.size(); ++i)
                                row.push_back(std::stod(fields[i]));
                        if (!row.empty())
                                results.atoms.push_back(std::move(row));
			continue;
		}

		if (fields.size()>5)
		{
			results.box[0][1] = std::stod(fields[3]);
			results.box[1][1] = std::stod(fields[4]);
			results.box[2][1] = std::stod(fields[8]);
			continue;
		}

	}
	inputFile.close();

        // Optional sanity check
        if (results.natoms != 0 && static_cast<int>(results.atoms.size()) != results.natoms)
                std::cerr << "Warning: expected " << results.natoms << " atoms but found " << results.atoms.size() << "\n";


	return results;

}

LAMMPSData readLAMMPSdatafile(const std::string& filename, int mode = 1)
{
	LAMMPSData results;
	results.box.resize(3, std::vector<double>(2, 0.0));
	int count = 0;

	std::ifstream inputFile(filename);
        if (!inputFile.is_open())
                throw std::runtime_error("Error: could not open file " + filename);

	std::string line;

	while(std::getline(inputFile,line))
	{
		if(line.empty()) continue;
		if(line.find("Velocities") != std::string::npos) break;

		std::istringstream iss(line);
		std::string word;
		std::vector<std::string> fields;
		while(iss >> word) fields.push_back(word);

		// number of atoms
        	if (fields.size() == 2 && fields[1] == "atoms")
		{
            		results.natoms = std::stoi(fields[0]);
            		continue;
        	}

		// number of atom types
		if (fields.size() == 2 && fields[1] == "atom")
		{
			results.ntypes = std::stoi(fields[0]);
			continue;
		}

		// box dimensions (xlo xhi, ylo yhi, zlo zhi)
        	if (fields.size() == 4 && (fields[2] == "xlo" || fields[2] == "ylo" || fields[2] == "zlo") && count <3)
		{
                	results.box[count][0] = std::stod(fields[0]);
                	results.box[count][1] = std::stod(fields[1]);
                	++count;
			continue;
            	}

		// parse atom positions and ids
		if ((mode == 1 && fields.size()==5) || (mode == 2 && fields.size()==8))
		{
			std::vector<double> row;
			for (size_t i = 0; i < fields.size(); ++i)
				row.push_back(std::stod(fields[i]));
			if (!row.empty())
                		results.atoms.push_back(std::move(row));
		}

	}

	inputFile.close();

	// Optional sanity check
    	if (results.natoms != 0 && static_cast<int>(results.atoms.size()) != results.natoms)
        	std::cerr << "Warning: expected " << results.natoms << " atoms but found " << results.atoms.size() << "\n";

    return results;
}


using namespace oILAB;

/*! Lagrange-Gauss reduction of a 2D basis.  Nearest-lattice-point rounding is only valid on a
 *  reduced basis: for an oblique cell the closest point can be two cells away in coefficient
 *  space, so an unreduced basis gives wrong minima. */
static void gaussReduce(Eigen::Vector3d& q1, Eigen::Vector3d& q2)
{
    for (int guard=0; guard<100; ++guard) {
        if (q2.squaredNorm() < q1.squaredNorm()) std::swap(q1,q2);
        const double mu= std::round(q1.dot(q2)/q1.squaredNorm());
        if (mu == 0.0) break;
        q2-= mu*q1;
    }
}

/*! Distance between two in-plane vectors on the torus spanned by \p q1 and \p q2, which must be
 *  Gauss-reduced. */
static double torusDistance(const Eigen::Vector3d& d,
                            const Eigen::Vector3d& q1, const Eigen::Vector3d& q2)
{
    Eigen::Matrix<double,3,2> Q; Q.col(0)= q1; Q.col(1)= q2;
    const Eigen::Vector2d c= Q.colPivHouseholderQr().solve(d);
    const Eigen::Vector2d centred(c(0)-std::round(c(0)), c(1)-std::round(c(1)));
    double best= 1.0e300;
    for (int i=-1;i<=1;++i)
        for (int j=-1;j<=1;++j)
            best= std::min(best, (Q*(centred - Eigen::Vector2d(i,j))).norm());
    return best;
}

/*! The atom positions of an extended-XYZ configuration.
 *
 *  Only the positions are wanted, so only they are read: an atom line of these files is
 *  `type x y z radius`, and the count on the first line is trusted over counting lines, which is
 *  what tells a truncated file from a complete one.  readXYZFile() above reads more of the format
 *  but recovers the box by counting tokens on the comment line, which the `Lattice="` field's own
 *  leading whitespace shifts; the box is not needed here and is taken from the mesostate. */
static std::vector<Eigen::Vector3d> configurationPositions(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("could not open configuration " + filename);
    int count= 0;
    file >> count;
    std::string rest;
    std::getline(file, rest);   // finish the count line
    std::getline(file, rest);   // the comment line, which is the mesostate's own header
    std::vector<Eigen::Vector3d> positions;
    positions.reserve(count);
    for (int i=0; i<count; ++i) {
        double type=0.0, x=0.0, y=0.0, z=0.0, radius=0.0;
        if (!(file >> type >> x >> y >> z >> radius))
            throw std::runtime_error("configuration " + filename + " ends after "
                                     + std::to_string(i) + " of " + std::to_string(count)
                                     + " atoms");
        positions.emplace_back(x,y,z);
    }
    return positions;
}

/*! What the overlap removal will find in a deformed configuration.
 *
 *  The enumeration engages a chosen set of coincidence points, but it reaches them with a
 *  displacement field, and a field cannot be aimed at the engaged pairs alone: wherever it
 *  happens to take the value that closes the gap between two other atoms of the two grains,
 *  those meet as well.  With one node engaged the field is uniform -- every node of the
 *  triangulation carries the same displacement -- so the grains translate rigidly and whole
 *  periodic sets of atoms arrive together; with several nodes the field varies and the extra
 *  coincidences are scattered instead.  Either way they are indistinguishable from the engaged
 *  ones: they sit at the same zero separation, the minimisation fuses them alike, and the
 *  boundary that results holds all of them.
 *
 *  So the number of coincidences a state realises is a property of the structure, while the
 *  number of nodes it engages is a property of the path the enumeration took to reach it.  They
 *  are not the same number, and the first is the one that describes the boundary. */
struct Coincidences
{
    int sites   = 0;   //!< points where two or more atoms meet
    int fused   = 0;   //!< atoms the overlap removal deletes: one fewer than meet at each site
    int largest = 0;   //!< how many atoms meet at the most crowded site
};

/*! The coincidences of a deformed configuration, counted the way the minimisation will fuse them.
 *
 *  \p period1 and \p period2 are the two periodic box vectors, which need not be reduced.
 *  Atoms meeting at one point are grouped rather than counted pairwise: three atoms at a point
 *  are one coincidence costing two atoms, not three coincidences costing three, and only
 *  grouping gets that right.
 *
 *  Pairs are found by sweeping a window along the boundary normal.  The separation of two atoms
 *  splits into an in-plane part, which the periodic images move, and a part along the normal,
 *  which they do not; since the two are orthogonal the full separation is the hypotenuse of the
 *  torus distance and the difference in height, and no pair further apart in height than the
 *  cutoff can be a coincidence whatever its in-plane position.  That bound is what keeps this
 *  off the critical path: the atoms lie on planes far further apart than the cutoff, so the
 *  window holds one plane at a time instead of the whole cell. */
static Coincidences countCoincidences(const std::vector<Eigen::Vector3d>& atoms,
                                      const Eigen::Vector3d& period1,
                                      const Eigen::Vector3d& period2,
                                      const double cutoff)
{
    Eigen::Vector3d q1= period1, q2= period2;
    gaussReduce(q1,q2);
    const Eigen::Vector3d nHat= q1.cross(q2).normalized();

    const int n= static_cast<int>(atoms.size());
    std::vector<Eigen::Vector3d> position(n);
    std::vector<double> height(n);
    std::vector<int> order(n);
    for (int i=0; i<n; ++i) {
        position[i]= atoms[i];
        height[i]= position[i].dot(nHat);
        order[i]= i;
    }
    std::sort(order.begin(), order.end(),
              [&height](const int a, const int b){ return height[a] < height[b]; });

    // Union-find over the coincident pairs, so that a point where several atoms meet is counted
    // once however many pairs it contributes.
    std::vector<int> parent(n);
    for (int i=0; i<n; ++i) parent[i]= i;
    const auto root= [&parent](int i)
    {
        while (parent[i]!=i) { parent[i]= parent[parent[i]]; i= parent[i]; }
        return i;
    };

    for (int a=0; a<n; ++a)
        for (int b=a+1; b<n && height[order[b]]-height[order[a]] < cutoff; ++b)
        {
            const int i= order[a], j= order[b];
            const Eigen::Vector3d d= position[i]-position[j];
            const double inPlane= torusDistance(d,q1,q2);
            const double normal= height[i]-height[j];
            if (std::sqrt(inPlane*inPlane + normal*normal) >= cutoff) continue;
            const int ri= root(i), rj= root(j);
            if (ri!=rj) parent[ri]= rj;
        }

    std::map<int,int> groupSize;
    for (int i=0; i<n; ++i) ++groupSize[root(i)];

    Coincidences realized;
    for (const auto& [representative,size] : groupSize) {
        (void)representative;
        if (size < 2) continue;
        ++realized.sites;
        realized.fused+= size-1;
        realized.largest= std::max(realized.largest, size);
    }
    return realized;
}

/*! Species written for the enumeration's coincidence sites when they are appended to a
 *  configuration: every site the sweep could have engaged, and the ones this state did.
 *
 *  Grain A is 1, grain B is 2, and the atoms a state brings into coincidence are
 *  GbMesoState<3>::coincidenceType = 3, so the sites continue from there.  They are markers, not
 *  atoms: they go only into the undeformed configuration, which nothing hands to LAMMPS, so they
 *  cannot affect an energy. */
static constexpr int siteType        = 4;   //!< a site the enumeration considered
static constexpr int engagedSiteType = 5;   //!< one this state engages

/*! Append the enumeration's coincidence sites to an extended-XYZ configuration.
 *
 *  A state's own engaged nodes are visible in the configuration already, as the type-3 atoms --
 *  but in the UNDEFORMED configuration those atoms sit where they started, at \f$x_A\f$ and
 *  \f$x_B\f$, not at the point \f$s\f$ where they will meet.  So the reference configuration on
 *  its own does not say which sites a state used, and says nothing at all about the sites it
 *  passed over.  Both are what a check of whether the sweep sampled the possibilities needs, so
 *  every site is written with every state and marked according to whether that state took it.
 *
 *  The atom count on the first line is rewritten rather than the rows simply appended: extended
 *  XYZ takes that count as authoritative and a reader stops there, so appending without it would
 *  leave the new rows in the file and invisible. */
static void appendSites(const std::string& path,
                        const std::vector<Eigen::Vector3d>& sites,
                        const std::vector<char>& engaged,
                        const double& radius)
{
    std::vector<std::string> lines;
    {
        std::ifstream configuration(path);
        if (!configuration.is_open())
            throw std::runtime_error("could not reopen configuration " + path);
        for (std::string line; std::getline(configuration, line); )
            lines.push_back(line);
    }
    if (lines.size() < 2)
        throw std::runtime_error("configuration " + path + " is too short to hold a frame");

    std::ofstream configuration(path);
    configuration << (std::stoi(lines[0]) + (int)sites.size()) << "\n";
    configuration << lines[1] << "\n";
    for (std::size_t i=2; i<lines.size(); ++i)
        if (!lines[i].empty()) configuration << lines[i] << "\n";
    configuration << std::fixed << std::setprecision(8);
    for (std::size_t i=0; i<sites.size(); ++i)
        configuration << (engaged[i] ? engagedSiteType : siteType)
                      << " " << sites[i](0) << " " << sites[i](1) << " " << sites[i](2)
                      << " " << radius << "\n";
}

/*! A directory name identifying the boundary a run is studying.
 *
 *  Sigma and the misorientation are what name a boundary in conversation, and they are what the
 *  folder is asked for.  The axis and the plane come along because they are what separate two
 *  boundaries that share those: sigma 5 alone does not say whether the plane is (210) or (310).
 *
 *  The Miller indices are joined with dots rather than run together.  Concatenation reads better
 *  for single digits -- "001" over "0.0.1" -- but this project already studies a boundary on
 *  (-5 5 14), and "-5514" could be read three ways. */
static std::string boundaryName(const int sigma, const double thetaDegrees,
                                const Eigen::Vector3d& axis, const Eigen::Vector3d& plane)
{
    const auto miller= [](const Eigen::Vector3d& v)
    {
        std::ostringstream o;
        for (int i=0; i<3; ++i) { if (i) o << '.'; o << (long long)std::llround(v(i)); }
        return o.str();
    };
    std::ostringstream o;
    o << "sigma" << sigma
      << "_theta" << std::fixed << std::setprecision(2) << thetaDegrees
      << "_axis" << miller(axis)
      << "_plane" << miller(plane);
    return o.str();
}

/*! Which of the three energies the shortlist is ranked on.
 *
 *  `Unrelaxed` is the configuration as constructed; `Tethered` is the relaxation that holds the
 *  boundary atoms where the construction put them, so it costs the state the enumeration built;
 *  `Full` is the free relaxation, so it costs the boundary that state leads to.  They rank
 *  differently and the choice is a question about what the sweep is for, not a detail. */
enum class RankBy { Unrelaxed, Tethered, Full };

static const char* rankName(const RankBy by)
{
    switch (by) {
        case RankBy::Unrelaxed: return "unrelaxed";
        case RankBy::Tethered:  return "tethered";
        case RankBy::Full:      return "full";
    }
    return "unknown";
}

/*! Everything pass 1 keeps about a state: its numbers, and the signature needed to build it
 *  again.  The configurations themselves are not kept -- that is the point of the two passes --
 *  so this is the whole of what the shortlist has to choose on. */
struct Surveyed
{
    //! The signature in sparse form: the ensemble indices this state engages.  The 0/1 vector is
    //! the indicator of this set, and runs to one entry per ensemble member, which over a sweep
    //! of this size is gigabytes of mostly zeros.
    std::vector<int> engaged;
    int nodes= 0;              //!< coincidences the deformed structure holds
    int fused= 0;              //!< atoms the overlap removal deletes
    int expelled= 0;           //!< atoms the deformation carried out of their own grain
    double corrugation= 0.0;
    double density= 0.0;
    double unrelaxed= 0.0;
    double tethered= 0.0;
    double spring= 0.0;
    double full= 0.0;

    double rankValue(const RankBy by) const
    {
        switch (by) {
            case RankBy::Unrelaxed: return unrelaxed;
            case RankBy::Tethered:  return tethered;
            case RankBy::Full:      return full;
        }
        return full;
    }
};

/*! The numbers that have to agree for two states to be one and the same.
 *
 *  Distinct placements of the engaged nodes that are translations or symmetry images of one
 *  another relax to the same structure, and nothing recorded tells them apart.  Every recorded
 *  number is included rather than the ranked one alone: two states really are the same only if
 *  nothing separates them, and a bare energy match would merge an accidental degeneracy between
 *  genuinely different structures.
 *
 *  `engaged` is deliberately left out.  It records the route the enumeration took rather than the
 *  structure it arrived at, and the sweep reaches one boundary by engaging different numbers of
 *  its coincidences -- the rest forming on their own -- so including it would preserve exactly
 *  the duplicates this exists to remove. */
static std::vector<double> stateFingerprint(const Surveyed& state)
{
    return { (double)state.nodes, (double)state.fused, state.corrugation, state.density,
             state.unrelaxed, state.tethered, state.spring, state.full };
}

/*! Group states whose fingerprints agree to within \p tolerance, returning a group index for
 *  each entry of \p members.
 *
 *  The comparison is a tolerance, not a rounded key.  A rounded key is cheaper but splits two
 *  values that fall either side of a bucket edge, and that is exactly what a minimiser's
 *  last-digit noise does: two copies of one structure whose free relaxations came back as
 *  5.363889534 and 5.363889495 differ by 4e-8 and land in different buckets at any tolerance
 *  whose edge lies between them.  Sorting on the fingerprint puts such states next to one
 *  another, so one pass over the sorted order groups them, and each is compared against its
 *  group's first member rather than its predecessor so that a long run of near-equal values
 *  cannot drift a group away from where it started. */
static std::vector<int> groupByFingerprint(const std::vector<const Surveyed*>& members,
                                           const double tolerance)
{
    std::vector<std::vector<double>> fingerprints(members.size());
    for (std::size_t i=0; i<members.size(); ++i)
        fingerprints[i]= stateFingerprint(*members[i]);

    std::vector<int> order(members.size());
    for (std::size_t i=0; i<order.size(); ++i) order[i]= (int)i;
    std::sort(order.begin(), order.end(),
              [&fingerprints](const int a, const int b)
              { return fingerprints[a] < fingerprints[b]; });

    std::vector<int> group(members.size(), -1);
    int groups= 0, representative= -1;
    for (const int i : order)
    {
        bool matches= representative >= 0;
        if (matches)
            for (std::size_t field=0; field<fingerprints[i].size(); ++field)
                if (std::abs(fingerprints[i][field]-fingerprints[representative][field])
                    > tolerance) { matches= false; break; }
        if (!matches) { representative= i; ++groups; }
        group[i]= groups-1;
    }
    return group;
}

int main()
{
    /*! [Types] */
    using VectorDimI = LatticeCore<3>::VectorDimI;
    using VectorDimD = LatticeCore<3>::VectorDimD;
    using Vector2d= Eigen::Vector2d;
    using IntScalarType = LatticeCore<3>::IntScalarType;



	VectorDimD axis(0, 0, 1);
	//double theta= 53.130102354155994249*std::numbers::pi/180;       // misorientation angle
	double theta = 36.869897645844041278*std::numbers::pi/180;       // misorientation angle
    //double theta = 43.602818972703637712*std::numbers::pi/180;       // misorientation angle
    //double theta = 28.072486935852960954*std::numbers::pi/180;       // misorientation angle
    //VectorDimD gbNormal(-2,1,0);                        // Miller indices
    VectorDimD gbNormal(-3,1,0);                        // Miller indices
	//VectorDimD gbNormal(-5,2,0);
    //VectorDimD gbNormal(-4,1,0);
    int heightScaling= 4;

    /*! The least the crystal may measure along the boundary normal, in Angstrom.
     *
     *  heightScaling counts CSL repeats, and a repeat is as long as the boundary makes it: four
     *  of them give 64.7 A on sigma5 (210) and only 45.7 A on sigma5 (310).  That is not a
     *  cosmetic difference.  The cohesive energy is read off a slab sitting
     *  gbHalfThickness+bulkSlabGap out from the middle -- 16 A as those are currently set -- and
     *  a 45.7 A crystal leaves that slab 1.9 A short of the free surface, so it samples
     *  under-coordinated atoms instead of bulk.  Measured on sigma5 (310) at axisScaling 3, the
     *  lowest relaxed boundary energy came out 0.311 J/m^2 that way against 0.905 J/m^2 once the
     *  crystal was thick enough: wrong by a factor of three, and quietly.
     *
     *  So heightScaling is a request rather than the last word.  Whatever it asks for, enough
     *  repeats are taken to reach this thickness, which decouples the answer from how long the
     *  boundary happens to make its repeat.  Note this is the crystal, not the LAMMPS box -- the
     *  box is this plus vacuumThickness at each end. */
    const double minimumCrystalThickness= 50.0;
	int periodScaling= 1;
	int axisScaling= 3;


    // enumerateStates: build every mesostate the ensemble's (t,s) pairs admit.  Set it false to
    // read ready-made signatures from `fin` instead.
    const bool enumerateStates            = true;
    // Flat = the original enumeration (CSL shifts along the boundary); Full = ball + slab;
    // Sites = coincidence points first, with the two grains displaced independently.
    const GbShiftSearch searchMode        = GbShiftSearch::Sites;
    // FLAT RUN.  Zero keeps the coincidence points exactly on the boundary plane, |s.n| = 0, so
    // every state this run builds has a planar boundary.
    const double slabHalfThickness        = 0.0;   // Angstrom
    // Choose how far the atoms in grain A and B can move to get to the selected CSL point
    const double dMax                     = 1.95;  // Angstrom
    // Keep only the nodes whose two grains move by equal and opposite amounts, u_A = -u_B, so
    // that the coincidence point is the midpoint of the pair of atoms that meets there.
    // False sweeps assymetric displacments as well.
    const bool symmetricDisplacementsOnly = true;
    // Run both relaxations in one LAMMPS invocation: relax against the tether, read the tethered
    // figures off, release the restraint, and relax again from where that left the structure.
    const bool chainRelaxations           = false;
    // Remove every coincidence the state did not engage -- both of the atoms that met there --
    // so that the configuration holds exactly the coincidences its signature names.
    // False keeps every atom and lets the relaxation fuse whatever the field brought together,
    // which is the honest setting for asking how many coincidences a displacement field really
    // creates, and the one whose energies describe the state the enumeration nominally built.
    const bool dropUnengagedCoincidences  = true;
    // Where results are written.  GBFACETMESHING_RUNS_DIR is <project>/runs, supplied by CMake:
    // outside the build tree, so that regenerating or clearing a build directory cannot reach
    // them.  "." writes them into the working directory instead, which for a run started from
    // the IDE is the directory holding the binary.
    const std::string outputRoot          = GBFACETMESHING_RUNS_DIR;
    // How large the coincidence-site markers are drawn in state_<index>_0.txt.  The atoms carry
    // 0.05, so a larger value picks the sites out of the structure rather than hiding them in it
    // -- but only just larger: at four times the atomic radius the markers sat over the boundary
    // they were meant to annotate and made the structure harder to read, not easier.
    const double siteMarkerRadius         = 0.10;
    // Exclude the nodes that do nothing but slide the grains along the tilt axis.
    //
    // A node whose two displacements both run along the axis, and both carry their atom a whole
    // atomic layer or more, brings nothing new into coincidence: it slides the atoms along the
    // axis until they meet at the coincidence point next door.  The state built on it is a state
    // the sweep already holds, with its atoms moved along the axis -- and those copies cost a
    // construction and two minimisations each while adding nothing to the sampling.
    //
    // They cannot be folded away afterwards on the numbers.  Translating the cell reorders the
    // neighbour lists LAMMPS builds, so the minimiser takes a different path and the copies come
    // back with energies differing in the fourth decimal, far above the tolerance that folds
    // genuine duplicates.  Excluded here, they are never built at all.
    const bool excludeAxialSlides         = false;
    // Two engaged coincidence points must be at least this far apart in projection onto the
    // boundary plane.  Zero admits any pair that is not exactly coincident.  The points lie on a
    // lattice with spacings 0.404, 0.808, 1.617 A here, so the criterion bites in steps: 0.45
    // removes the 0.404 A pairs -- the ones that drive ordinary, non-node atoms of the two grains
    // to within half an Angstrom of each other -- and keeps everything at 0.808 A and beyond.
    const double minSiteSeparation        = 0.45;  // Angstrom
    // Keep only the shifts whose coincidence point lies exactly in the boundary plane.  The
    // deformed surface is then planar and every state built over them is a flat STGB -- the
    // baseline the faceted states depart from.
    //
    // Redundant for a Sites run with slabHalfThickness at zero, which already confines the
    // points to the plane.  It is how to get the flat suite out of a run whose slab is not zero,
    // and it is the only way to get it out of a Flat or Full run: those take their points from a
    // CSL cell of half-thickness sPerpMax*b/2, which is 1 A here, not from the plane.
    const bool restrictToFlatStates       = false;
    // How many (t,s) pairs one state may engage at once.  0 means no limit, in which case the
    // only ceiling is the clash rule itself: no state can engage two pairs that share a lattice
    // site, so no state can be larger than the smaller of the two site counts.  Every state the
    // walk produces is clash-free, so this is a count of real states, not of candidates.
    //
    // No limit here: the whole suite is wanted.  The clash rule stops the walk at 18 engaged
    // nodes -- the 24 sites draw on only 18 atoms of each grain -- and the sizes past 8 are more
    // than half the sweep (581031 of 921599), so capping at 8 would leave most of it unsampled.
    const int maxEngaged                  = 0;
    // ---- shortlist -------------------------------------------------------------
    // The sweep runs in two passes.  The first visits every state and keeps only its numbers, in
    // output_thread_<id>.txt; the second rebuilds the shortlisted states alone and writes their
    // configurations and relaxed structures.  Writing a configuration evaluates the displacement
    // field at every atom, so it has to happen in pass 1 anyway for LAMMPS to be handed a file --
    // but keeping the result is what made the sweep unaffordable, at tens of gigabytes of files
    // that nobody opens.  Only the shortlist is kept.
    //
    // How many states to keep per level of engagement, where the level is the number of
    // coincidences the structure realised, not the number the enumeration engaged.
    const int shortlistTop                = 40;
    // How many states to build from the ranking over the whole sweep, for each of the three
    // energies, on top of the per-level shortlist.
    //
    // The per-level shortlist deliberately refuses to compare across levels, which is right for
    // asking what the best boundary at each complexity is -- but it means the globally cheapest
    // states are not necessarily built, since a level's 41st state can undercut another level's
    // 1st.  These fill that in, and there is one set per energy because the three disagree about
    // which states are good.  Zero switches them off.
    const int globalTop                   = 100;
    // Which energy the ranking is taken on.  Tethered continues what the earlier runs ranked on;
    // Full ranks by the boundary each state relaxes into.
    const RankBy shortlistRankBy          = RankBy::Tethered;
    // Fold together the states the sweep records identically -- translations and symmetry images
    // of one structure.  Without this a single family fills the whole top `shortlistTop` and
    // hides the next distinct state; at one engaged node the families run to 24 copies.
    const bool shortlistDedup             = true;
    const double shortlistDedupTol        = 1.0e-6;
    // Refuse to start a run longer than this many states.  Each one is a mesostate construction
    // and, with energies on, a LAMMPS minimization.
    //
    // The settings above enumerate 921599 states; this leaves room to widen them a little
    // without the run refusing to start, while still catching a setting that blows the count up.
    const long long maxStates             = 1200000;
    const double tMax                     = 0.99;   // ball radius, in units of b
    const double tPerpMax                 = 0.99;   // slab half-thickness, in units of b
    // The translations that give well spread out flat sites point along the GB normal, with
    // |t.n| = 1.617 A here, so a slab of 0.5 b = 1.278 A cuts them out and the sites collapse onto
    // two lattice points.  1.0 b = 2.556 A clears them.
    const double layerHalfThickness       = 1.0;  // Angstrom either side of the flat GB
    const std::string fin                 = "sigma5.txt";

    // ---- LAMMPS ----------------------------------------------------------------
    const bool computeEnergies            = true;
    const bool minimizeInLammps           = true;  // relax before reading the energy
    const double tetherHalfWidth          = 4.0;   // Angstrom; 0 = full minimisation
    const double tetherStiffness          = 10.0;  // eV/Angstrom^2
    const int numThreads                  = 100;
    const std::string potentialName       = "Cu_mishin1.eam.alloy";
    const std::string lmpLocation         = "/usr/bin/lmp";
    // -------------------------------------------------------------------------------

    std::cout << "states = "
              << (slabHalfThickness > 0.0 ? "general (every shift, faceted boundaries allowed)"
                                          : "flat (coincidence points on the boundary plane)")
              << (symmetricDisplacementsOnly ? ", symmetric displacements only" : "") << std::endl;
    std::cout << "search = " << (searchMode==GbShiftSearch::Flat  ? "Flat (original rule)" :
                                 searchMode==GbShiftSearch::Full  ? "Full (ball + slab)"
                                                                  : "Sites (coincidence points "
                                                                    "first)") << std::endl;

    // A missing executable or potential is worth catching here: energy() would otherwise run a
    // command that does nothing and then read an output file that was never written.
    bool energiesRequested = computeEnergies;
    if (energiesRequested && !std::filesystem::exists(lmpLocation)) {
        std::cout << "LAMMPS executable not found at " << lmpLocation
                  << " -- writing configurations only, no energies." << std::endl;
        energiesRequested = false;
    }
    if (energiesRequested && !std::filesystem::exists(potentialName)) {
        std::cout << "Potential file " << potentialName << " not found in "
                  << std::filesystem::current_path().string()
                  << " -- writing configurations only, no energies." << std::endl;
        energiesRequested = false;
    }
    if (energiesRequested)
        std::cout << "energies = LAMMPS at " << lmpLocation << ", potential " << potentialName
                  << (minimizeInLammps ? ", minimized" : ", unrelaxed") << std::endl;
    if (energiesRequested && minimizeInLammps)
        std::cout << "relaxation = "
                  << (tetherHalfWidth > 0.0
                      ? "constrained (atoms within " + std::to_string(tetherHalfWidth)
                        + " A of the boundary tethered, k = " + std::to_string(tetherStiffness)
                        + " eV/A^2)"
                      : "free (no tether)") << std::endl;

    // The original path expects `fin`; the enumerating path does not read it at all.
    std::vector<std::vector<int>> stateData;
    int rows = 0, stateSize = 0;
    if (!enumerateStates) {
        std::cout << "Reading the file: " << fin << std::endl;
        stateData = readPostProcessingOutput(fin);
        rows = stateData.size();
        stateSize = stateData[0].size()-3;
    }

    /*! [Lattice] */
    double strain=0.0; // 0.01 or 0.02
    double a0= (1.0+strain)*3.615;
    Eigen::Matrix3d A;
    A << 0.0, 0.5, 0.5,
            0.5, 0.0, 0.5,
            0.5, 0.5, 0.0;
    A= a0*A;
    Lattice<3> lattice(A);
    std::cout << "Lattice A = " << std::endl;
    std::cout << lattice.latticeBasis << std::endl;
    /*! [Lattice] */

    ReciprocalLatticeVector<3> rAxisGlobal(lattice.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
    std::cout << "Cartesian coordinates of axis = " << std::endl;
    std::cout << rAxisGlobal.cartesian().transpose() << std::endl;

    try
    {
        // construct bicrystal
        Eigen::AngleAxis<double> halfRotation(theta/2,rAxisGlobal.cartesian().normalized());
        Lattice<3> latticeA(lattice.latticeBasis,halfRotation.matrix());
        Lattice<3> latticeB(lattice.latticeBasis,halfRotation.matrix().transpose());
        BiCrystal<3> bc(latticeA,latticeB,false);
        std::cout << "Sigma = " << bc.sigma << std::endl;

        // The boundary plane as it was asked for.  The line below turns gbNormal into the
        // rotated, Cartesian normal the construction needs, and the Miller indices are gone after
        // that -- but they are what names the boundary, so they are kept first.
        const VectorDimD gbNormalMiller= gbNormal;

        // construct GB
        gbNormal= halfRotation.matrix() * gbNormal;
        ReciprocalLatticeDirection<3> rd= latticeA.reciprocalLatticeDirection(gbNormal);
        Gb<3> gb(bc,rd);


        // Define the CSL box vectors
        ReciprocalLatticeVector<3> rAxisA(latticeA.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
        std::cout << gb.getPeriodVector(rAxisA).cartesian().transpose() << std::endl;
        LatticeVector<3> axisA(gb.bc.A.latticeDirection(axis).latticeVector());
        LatticeVector<3> axisC(gb.bc.getLatticeDirectionInC(axisA).latticeVector());
        std::vector<LatticeVector<3>> cslVectors;
        // Enough CSL repeats along the normal to satisfy both heightScaling and the floor above.
        // Gb::box() lays the cell out from -boxVectors[0] to +boxVectors[0], one grain on each
        // side, so the crystal measures twice this vector -- hence the 2 here.  Getting that
        // factor wrong is not harmless in either direction: too few repeats is the error the
        // floor exists to prevent, and too many multiplies the cost of every state in the sweep.
        const LatticeVector<3> heightRepeat=
            gb.bc.csl.latticeDirection(gb.nA.cartesian()).latticeVector();
        const double crystalPerRepeat= 2.0*heightRepeat.cartesian().norm();
        const int heightRepeatsNeeded=
            static_cast<int>(std::ceil(minimumCrystalThickness/crystalPerRepeat));
        const int effectiveHeightScaling= std::max(heightScaling, heightRepeatsNeeded);
        std::cout << "crystal along the normal: " << effectiveHeightScaling << " CSL repeat(s), "
                  << crystalPerRepeat << " A of crystal each = "
                  << effectiveHeightScaling*crystalPerRepeat << " A";
        if (effectiveHeightScaling > heightScaling)
            std::cout << "  (heightScaling " << heightScaling << " would have given only "
                      << heightScaling*crystalPerRepeat << " A, under the "
                      << minimumCrystalThickness << " A floor)";
        std::cout << std::endl;
        cslVectors.push_back(effectiveHeightScaling*heightRepeat);
        cslVectors.push_back(periodScaling*gb.getPeriodVector(rAxisA));
        cslVectors.push_back(axisScaling*axisC);
        gb.box(cslVectors,1,"gb.txt",false);
        bc.box(cslVectors,1,"bcOriented.txt",true);


        // sPerpMax is in units of b and bounds |s.n| by sPerpMax*b/2, so a layer of
        // +-layerHalfThickness about the flat boundary corresponds to 2*layerHalfThickness/b.
        const double b = gb.bc.A.latticeBasis.col(0).norm();
        const double sPerpMax = 2.0*layerHalfThickness/b;
        std::cout << "b = " << b << ",  shifts confined to +-" << sPerpMax*b/2
                  << " A of the flat GB" << std::endl;

        GbMesoStateEnsemble<3> ensemble(gb, rAxisA, cslVectors,
                                        tMax, sPerpMax,
                                        searchMode, tPerpMax, false,
                                        searchMode==GbShiftSearch::Sites ? "nodes.txt"
                                                                         : "translationsNonFlat.txt",
                                        slabHalfThickness, dMax);
        const int ensembleSize = ensemble.initializeState().size();
        std::cout << "Size of the ensemble = " << ensembleSize << std::endl;

        if (enumerateStates)
        {
            // The boundary normal, kept for reporting how far off the flat plane each engaged
            // shift sits -- that offset is what the general states are about.
            const Eigen::Vector3d nHat=
                gb.bc.getReciprocalLatticeDirectionInC(gb.nB.reciprocalLatticeVector())
                    .cartesian().normalized();

            // The general run takes every (t,s) pair the ensemble holds; symmetricDisplacements-
            // Only keeps just those whose two grains move by equal and opposite amounts.  The
            // test is u_A + u_B = 0, which says the coincidence point sits at the midpoint of
            // the two atoms rather than anywhere on the segment between them.  Only the Sites
            // search can produce anything else -- Flat and Full derive the point as that
            // midpoint -- so the filter is meaningful there and inert elsewhere.
            // The atomic layers normal to the tilt axis stand 1/|rAxisA| apart -- an
            // interplanar spacing -- and that is the shortest distance an atom can be carried
            // along the axis and land on the layer next door.  The shortest lattice vector along
            // the axis is generally longer: in FCC along [001] the layers are a/2 apart while
            // the lattice itself does not repeat until a.  It is the layer spacing that matters
            // here, because landing on the next layer is what makes the move a slide.
            const Eigen::Vector3d axisHat= axisC.cartesian().normalized();
            const double layerSpacing= 1.0/rAxisA.cartesian().norm();
            std::cout << "tilt axis: atomic layers " << layerSpacing
                      << " A apart, lattice A repeats after " << axisA.cartesian().norm()
                      << " A, CSL after " << axisC.cartesian().norm() << " A" << std::endl;

            /*! Does this displacement run along the tilt axis and carry its atom a whole layer? */
            const auto slidesAlongAxis= [&axisHat,layerSpacing](const VectorDimD& u)
            {
                const double along= u.dot(axisHat);
                return (u - along*axisHat).norm() < 1.0e-6           // nothing off the axis
                    && std::abs(along) >= layerSpacing - 1.0e-6;     // at least one layer
            };

            std::vector<int> family;
            std::ostringstream slidesDropped;
            int asymmetric= 0, axialSlides= 0;
            for (int i=0; i<ensembleSize; ++i)
            {
                if (symmetricDisplacementsOnly && searchMode==GbShiftSearch::Sites
                    && (ensemble.nodes[i].uA() + ensemble.nodes[i].uB()).norm() > 1.0e-8) {
                    ++asymmetric;
                    continue;
                }
                // Both displacements have to be slides for the node to be one: a node that
                // slides one grain along the axis while moving the other across it is bringing
                // together atoms that no translation would.
                if (excludeAxialSlides && searchMode==GbShiftSearch::Sites
                    && slidesAlongAxis(ensemble.nodes[i].uA())
                    && slidesAlongAxis(ensemble.nodes[i].uB())) {
                    ++axialSlides;
                    // Named rather than counted: this rule removes candidates, which lowers the
                    // engagement ceiling and so reshapes the whole sweep, and it should be
                    // possible to see that it removed the nodes it was meant to and no others.
                    const Eigen::Vector3d uA= ensemble.nodes[i].uA();
                    slidesDropped << "    node " << std::setw(4) << i
                                  << "   |uA| = " << std::fixed << std::setprecision(4)
                                  << uA.norm()
                                  << "   uA.axis = " << uA.dot(axisHat)
                                  << " A  = " << std::setprecision(2)
                                  << uA.dot(axisHat)/layerSpacing << " layer(s)\n";
                    continue;
                }
                family.push_back(i);
            }
            if (symmetricDisplacementsOnly && searchMode==GbShiftSearch::Sites)
                std::cout << "symmetric displacements only (uA = -uB) : dropped "
                          << asymmetric << " asymmetric node(s)" << std::endl;
            if (excludeAxialSlides && searchMode==GbShiftSearch::Sites) {
                std::cout << "axial slides excluded : dropped " << axialSlides
                          << " node(s) whose uA and uB both run along the tilt axis by at least "
                          << layerSpacing << " A\n" << slidesDropped.str() << std::flush;
            }

            // The flat suite: drop every candidate whose coincidence point sits off the
            // boundary plane.  Where the point is read from depends on the search -- Sites holds
            // coincidence nodes, the other two hold (t,s) pairs -- so both are covered here.
            if (restrictToFlatStates)
            {
                const std::size_t before= family.size();
                std::vector<int> flat;
                for (const int i : family) {
                    const Eigen::Vector3d point= (searchMode==GbShiftSearch::Sites)
                        ? ensemble.nodes[i].s : ensemble.tShiftPairs[i].second;
                    if (std::abs(point.dot(nHat)) < 1.0e-9) flat.push_back(i);
                }
                family= std::move(flat);
                std::cout << "flat states only : kept " << family.size() << " of " << before
                          << " candidate(s) whose coincidence point lies in the boundary plane"
                          << std::endl;
                if (family.empty())
                    throw std::runtime_error("No shift lies in the boundary plane, so there is "
                                             "no flat GB to build. Raise tMax or tPerpMax.");
            }

            std::cout << "\ncandidates considered: "
                      << family.size() << " of " << ensembleSize << std::endl;
            if (family.empty())
                throw std::runtime_error("The ensemble is empty. Widen the slab, raise dMax, or "
                                         "raise tMax/tPerpMax for the other searches.");

            // ---- the sites each candidate occupies ----------------------------------------
            // A candidate puts one node on an atom of A and one on an atom of B, and both are
            // decided by the candidate alone.  GbMesoState rejects any state that engages two
            // candidates sharing either atom, so knowing them up front is enough never to
            // enumerate such a state.  Both matter: two candidates can share a B atom while
            // their A atoms differ.
            //
            // A third key covers the coincidence points themselves.  Two that project onto the
            // same point of the boundary plane are one point of the torus the facets are
            // triangulated on, and the mesh collapses; two that are merely close leave the facet
            // rising and falling steeply over a short in-plane distance, which drives ordinary
            // atoms of the two grains together far more tightly than any lattice spacing.  The
            // walk therefore keeps engaged coincidence points at least minSiteSeparation apart,
            // which subsumes the exact case.
            std::vector<int> siteA(family.size()), siteB(family.size()), projected(family.size());
            std::vector<int> basisPairs;
            std::vector<VectorDimD> classPosition;   // one coincidence point per projected class
            // Where every candidate would bring its two atoms together, indexed by ensemble
            // index so that a state's engaged list can be turned into a set of sites.
            std::vector<VectorDimD> siteOfCandidate(ensembleSize, VectorDimD::Zero());
            {
                // The coincidence point projected onto the flat boundary plane, reduced into one
                // in-plane period cell.  Two nodes landing on the same projected point are the
                // same point of the torus the facets are triangulated on, so they collapse the
                // mesh; excluding them here is cheaper and more informative than discovering it
                // when the triangulation fails.
                Eigen::Matrix<double,3,2> inPlane;
                inPlane.col(0)= cslVectors[1].cartesian();
                inPlane.col(1)= cslVectors[2].cartesian();
                const auto inPlaneSolver= inPlane.colPivHouseholderQr();
                std::map<std::pair<long long,long long>,int> projectedClasses;
                std::vector<VectorDimD> projectedRepresentative;
                const auto projectedKey= [&](const VectorDimD& point)
                {
                    Eigen::Vector2d c= inPlaneSolver.solve(point);
                    for (int j=0; j<2; ++j) c(j)-= std::floor(c(j) + 1.0e-9);
                    return std::make_pair((long long)std::llround(c(0)*1.0e6),
                                          (long long)std::llround(c(1)*1.0e6));
                };

                std::vector<LatticeVector<3>> bicrystalBoxVectors(cslVectors);
                bicrystalBoxVectors[0]= 2*cslVectors[0];
                VectorDimD wrapShift;
                wrapShift << -0.5-FLT_EPSILON, -FLT_EPSILON, -FLT_EPSILON;

                std::map<OrderedTuplet<3>,int> aClasses, bClasses;
                std::vector<int> unplaceable;
                for (const int i : family)
                {
                    try {
                        OrderedTuplet<3> keyA, keyB;
                        if (searchMode==GbShiftSearch::Sites) {
                            // An atom and its periodic images are one atom, so the key is taken
                            // after wrapping into the bicrystal box.
                            VectorDimD xA= ensemble.nodes[i].xA.cartesian();
                            VectorDimD xB= ensemble.nodes[i].xB.cartesian();
                            LatticeVector<3>::modulo(xA, bicrystalBoxVectors, wrapShift);
                            LatticeVector<3>::modulo(xB, bicrystalBoxVectors, wrapShift);
                            keyA << gb.bc.A.latticeVector(xA);
                            keyB << gb.bc.B.latticeVector(xB);
                        }
                        else {
                            const auto placement= GbMesoState<3>::nodePlacement(
                                gb, cslVectors, ensemble.tShiftPairs[i].first,
                                ensemble.tShiftPairs[i].second);
                            keyA= placement.siteA;
                            keyB= placement.siteB;
                        }
                        const VectorDimD coincidence=
                            (searchMode==GbShiftSearch::Sites) ? ensemble.nodes[i].s
                                                               : ensemble.tShiftPairs[i].second;
                        const auto pk= projectedKey(coincidence);
                        const int a= aClasses.emplace(keyA, (int)aClasses.size()).first->second;
                        const int b= bClasses.emplace(keyB, (int)bClasses.size()).first->second;
                        const auto inserted= projectedClasses.emplace(pk, (int)projectedClasses.size());
                        if (inserted.second) projectedRepresentative.push_back(coincidence);
                        const int p= inserted.first->second;
                        siteA[basisPairs.size()]= a;
                        siteB[basisPairs.size()]= b;
                        projected[basisPairs.size()]= p;
                        siteOfCandidate[i]= coincidence;
                        basisPairs.push_back(i);
                    }
                    catch (const std::runtime_error&) {
                        // A node that misses its lattice cannot take part in any state.  Dropping
                        // it here is what keeps the sweep alive: the construction would otherwise
                        // hit the same failure with nothing useful to do about it.
                        unplaceable.push_back(i);
                    }
                }
                siteA.resize(basisPairs.size());
                siteB.resize(basisPairs.size());
                projected.resize(basisPairs.size());
                if (!unplaceable.empty())
                    std::cout << "  dropped " << unplaceable.size()
                              << " candidate(s) whose nodes miss their lattice" << std::endl;
                std::cout << "distinct atoms used : " << aClasses.size() << " of A, "
                          << bClasses.size() << " of B" << std::endl;
                std::cout << "distinct projected coincidence points : "
                          << projectedClasses.size() << std::endl;
                classPosition= projectedRepresentative;
            }

            // The distinct sites, one marker per position rather than one per candidate: two
            // candidates can bring different pairs of atoms to the same point, and drawing that
            // point twice would say there are more sites there than there are.
            // The sites are folded into the same period cell the configuration's atoms are
            // folded into.  box() wraps every atom along the two periodic box vectors before
            // writing it, so a site left where the enumeration found it can land on the far face
            // of the cell -- a point at z = 7.23 when the period is 7.23 belongs at z = 0 -- or
            // outside it altogether, and would then be drawn away from the atoms it describes.
            // The same snap-before-floor is used, so a site that belongs at the origin does not
            // end up a rounding error short of the opposite face.
            Eigen::Matrix<double,3,2> sitePeriods;
            sitePeriods.col(0)= cslVectors[1].cartesian();
            sitePeriods.col(1)= cslVectors[2].cartesian();
            const auto siteSolver= sitePeriods.colPivHouseholderQr();
            const auto wrapSite= [&sitePeriods,&siteSolver](const VectorDimD& point)
            {
                Eigen::Vector2d c= siteSolver.solve(point);
                Eigen::Vector2d whole;
                for (int j=0; j<2; ++j) whole(j)= std::floor(c(j) + FLT_EPSILON);
                return VectorDimD(point - sitePeriods*whole);
            };

            std::vector<VectorDimD> distinctSites;
            std::vector<int> siteOfIndex(ensembleSize, -1);
            {
                // Deduplication happens after wrapping, so two candidates whose coincidence
                // points differ by a period count as the one site they are.
                std::map<std::array<long long,3>,int> seen;
                for (const int i : basisPairs) {
                    const VectorDimD point= wrapSite(siteOfCandidate[i]);
                    const std::array<long long,3> key{
                        (long long)std::llround(point(0)*1.0e6),
                        (long long)std::llround(point(1)*1.0e6),
                        (long long)std::llround(point(2)*1.0e6)};
                    const auto inserted= seen.emplace(key, (int)distinctSites.size());
                    if (inserted.second) distinctSites.push_back(point);
                    siteOfIndex[i]= inserted.first->second;
                }
                std::cout << "coincidence sites considered : " << distinctSites.size()
                          << " distinct position(s) from " << basisPairs.size()
                          << " candidate(s); written into every state_<index>_0.txt as species "
                          << siteType << " (passed over) and " << engagedSiteType
                          << " (engaged)" << std::endl;
            }

            const int n= static_cast<int>(basisPairs.size());
            if (n==0)
                throw std::runtime_error("No (t,s) pair places both of its nodes on a lattice.");

            // No state can engage two candidates sharing an atom of either grain, nor two whose
            // coincidence points project onto the same point of the boundary plane, so no state
            // can be larger than the smallest of the three class counts -- asking for more is
            // asking for nothing.
            int distinctA= 0, distinctB= 0, distinctProjected= 0;
            for (int i=0; i<n; ++i) {
                distinctA= std::max(distinctA, siteA[i]+1);
                distinctB= std::max(distinctB, siteB[i]+1);
                distinctProjected= std::max(distinctProjected, projected[i]+1);
            }
            const int ceiling= std::min({distinctA, distinctB, distinctProjected});

            // Which projected classes a chosen one rules out.  The closed neighbourhood of a
            // class is itself together with every class within minSiteSeparation of it on the
            // torus, so blocking the neighbourhood is exactly the pairwise separation condition.
            // With the separation at zero this degenerates to the class itself, which is the
            // exact-coincidence key.
            Eigen::Vector3d reduced1= cslVectors[1].cartesian(), reduced2= cslVectors[2].cartesian();
            gaussReduce(reduced1, reduced2);
            std::vector<std::vector<int>> blocks(distinctProjected);
            {
                long long conflicting= 0;
                for (int i=0; i<distinctProjected; ++i) {
                    blocks[i].push_back(i);
                    for (int j=0; j<distinctProjected; ++j)
                        if (j!=i && torusDistance(classPosition[i]-classPosition[j],
                                                  reduced1, reduced2) < minSiteSeparation) {
                            blocks[i].push_back(j);
                            ++conflicting;
                        }
                }
                std::cout << "minimum separation of coincidence points : " << minSiteSeparation
                          << " A  (" << conflicting/2 << " conflicting class pair(s))" << std::endl;
            }
            int engagedLimit= (maxEngaged > 0 && maxEngaged < n) ? maxEngaged : n;
            if (engagedLimit > ceiling) {
                std::cout << "engaging at most " << ceiling
                          << " candidates at a time (no state can hold more without a collision; "
                             "maxEngaged asked for "
                          << (maxEngaged > 0 ? std::to_string(maxEngaged) : std::string("all"))
                          << ")" << std::endl;
                engagedLimit= ceiling;
            }

            // ---- enumerate the clash-free states ------------------------------------------
            // Walk the subsets as before, but carry the sites already spoken for and refuse to
            // extend with a pair that would reuse one.  Every subset of a clash-free set is
            // itself clash-free, so every node of this walk is a state worth building: nothing is
            // generated and thrown away, and the walk visits exactly the clash-free states.
            std::vector<std::vector<int>> subsets;
            {
                std::vector<char> usedA(distinctA,0), usedB(distinctB,0);
                // counts, not a flag: a class stays blocked while any chosen node still blocks it
                std::vector<int> blockedProjected(distinctProjected,0);
                std::vector<int> combination;
                bool capped= false;
                std::function<void(int)> walk = [&](int start)
                {
                    if (!combination.empty()) subsets.push_back(combination);
                    if (subsets.size() >= (std::size_t)maxStates) { capped= true; return; }
                    if (static_cast<int>(combination.size()) == engagedLimit) return;
                    for (int k=start; k<n && !capped; ++k)
                    {
                        if (usedA[siteA[k]] || usedB[siteB[k]] || blockedProjected[projected[k]])
                            continue;
                        usedA[siteA[k]]= 1; usedB[siteB[k]]= 1;
                        for (const int c : blocks[projected[k]]) ++blockedProjected[c];
                        combination.push_back(basisPairs[k]);
                        walk(k+1);
                        combination.pop_back();
                        for (const int c : blocks[projected[k]]) --blockedProjected[c];
                        usedA[siteA[k]]= 0; usedB[siteB[k]]= 0;
                    }
                };
                walk(0);
                if (capped)
                    throw std::runtime_error("More than maxStates (" + std::to_string(maxStates) +
                        ") clash-free states. Lower maxEngaged or raise maxStates.");
            }
            std::cout << "states to build : " << subsets.size()
                      << "   (clash-free by construction)" << std::endl;

            // Keep each run's output separate.  The undeformed and deformed configurations go in
            // their own folders, named by the state index alone: a viewer picks a set of files up
            // as one sequence by the number in their names, so the index has to be the only one
            // there.  Anything else in the name -- the site count, or box()'s own "reference0" /
            // "reference1" suffix -- gives it a second number to latch onto and the sequence is
            // not recognised.  The site count is recorded in states.txt instead.
            // The engaged-node limit is part of the directory name: runs at different limits are
            // nested rather than disjoint -- a limit of 4 contains everything a limit of 3 does --
            // so they are worth keeping side by side, and they must not overwrite one another.
            // Output goes under a folder named for the boundary, with one folder per run inside
            // it: several sweeps of one boundary at different settings then sit side by side
            // instead of overwriting each other or being told apart by a suffix.
            const std::string boundaryDirectory=
                boundaryName(gb.bc.sigma, theta*180.0/std::numbers::pi, axis, gbNormalMiller);
            // Results live outside the build tree.  A build directory is something a tool may
            // regenerate or clear, and results are not; keeping them apart means no cleanup of
            // one can reach the other.  Only LAMMPS's own per-thread scratch stays in the
            // working directory, because LAMMPS writes it there.
            const std::string outputDirectory=
                outputRoot + "/" + boundaryDirectory + "/"
                + (searchMode==GbShiftSearch::Flat ? "flat" : "full")
                + "_maxEngaged"
                + (maxEngaged > 0 ? std::to_string(maxEngaged) : std::string("all"));
            // Where pass 1 puts the configurations it needs but does not keep.  One pair of files
            // per thread, overwritten by every state that thread visits: LAMMPS has to be handed
            // a file, but nothing downstream wants 921599 of them.
            const std::string scratchDirectory = outputDirectory + "/scratch";
            for (const auto& directory : {outputDirectory, scratchDirectory})
                std::filesystem::create_directories(directory);

            std::cout << "writing to " << std::filesystem::absolute(outputDirectory).string()
                      << std::endl;
            std::cout << "threads = " << numThreads << std::endl;

            auto indexName= [](int i){
                std::ostringstream o; o << std::setw(3) << std::setfill('0') << i; return o.str(); };

            // The signature written out in full: one entry per ensemble member, 0 or 1.  The
            // sparse form beside it -- just the engaged indices -- is the same information and a
            // fraction of the bytes, but only the full vector can be read straight back: it is
            // what constructMesoState() takes, what the screen reports, and what the
            // enumerateStates=false path parses out of `fin`.  Built here from the engaged set
            // rather than streamed from the XTuplet so that the spacing is fixed and parseable
            // whatever that operator does.
            const auto fullSignature= [ensembleSize](const std::vector<int>& engaged)
            {
                std::vector<char> bits(ensembleSize,'0');
                for (const int i : engaged) bits[i]= '1';
                std::ostringstream o;
                for (int i=0; i<ensembleSize; ++i) { if (i) o << ' '; o << bits[i]; }
                return o.str();
            };

            // One line describing an engaged node of a state that has been built.  The values
            // are read off the mesostate rather than off the candidate list, because the two
            // differ: the construction subtracts the common-mode translation, which moves the
            // coincidence point and both displacements.  Reporting the candidate would describe
            // something that was never built.  Reading the constructed pairs also makes the line
            // identical for every search -- a (t,s) pair arrives here as its own u_A = +t/2,
            // u_B = -t/2.
            const auto describe= [&nHat](const auto& mesostate, const std::size_t position)
            {
                const auto& [xA,uA]= mesostate.xuPairsOfFacetedSurfaces.first[position];
                const auto& [xB,uB]= mesostate.xuPairsOfFacetedSurfaces.second[position];
                const Eigen::Vector3d coincidence= xA+uA;
                std::ostringstream o;
                // The displacements go out as vectors, not just magnitudes.  Two nodes can
                // agree on the coincidence point and on |uA| and still be different nodes,
                // pairing different atoms because the displacement points elsewhere -- and then
                // they build different boundaries.  Only the vectors tell them apart, which is
                // what makes two enumerations comparable node by node.
                o << std::fixed << std::setprecision(4)
                  << "s=(" << coincidence(0) << "," << coincidence(1) << "," << coincidence(2)
                  << ")  s.n=" << coincidence.dot(nHat)
                  << "  uA=(" << uA(0) << "," << uA(1) << "," << uA(2)
                  << ")  uB=(" << uB(0) << "," << uB(1) << "," << uB(2)
                  << ")  |uA|=" << uA.norm() << "  |uB|=" << uB.norm()
                  << "  |t|=" << (uA-uB).norm();
                return o.str();
            };

            // How far the deformed surface departs from the flat boundary.  This has to be asked
            // of the deformed surface: the reference facets are stepped, because t varies between
            // sites.  Reported, not enforced -- a non-zero value is the faceting itself.  The
            // out-of-plane offset is kept alongside the corrugation because a state engaging one
            // node has a flat boundary wherever it sits, and reporting only the offset would let
            // a displaced plane read as a faceted one.
            const auto surfaceExtent= [&nHat](const auto& mesostate)
            {
                const Eigen::MatrixXd deformed= mesostate.facetA.deformedVertices();
                double lowest= 1.0e300, highest= -1.0e300;
                for (int r=0; r<deformed.rows(); ++r) {
                    const double height= Eigen::Vector3d(deformed.row(r)).dot(nHat);
                    lowest= std::min(lowest, height);
                    highest= std::max(highest, height);
                }
                return std::make_pair(lowest,highest);
            };

            int rejected=0;
            std::map<std::string,int> reasons;

            // ================================================================ PASS 1
            // Visit every state, keep its numbers, throw its configuration away.
            // Pass 2 builds the shortlist for real.
            std::cout << "\npass 1 of 2 : surveying " << subsets.size()
                      << " state(s), keeping numbers only" << std::endl;

            // Where pass 1 spends its time.  One accumulator per thread, so the additions need
            // no synchronisation and the totals are core-seconds rather than wall-clock: that is
            // what says which phase to attack, since wall-clock hides everything behind whichever
            // phase happens to be running when a thread stalls.
            struct Phases { double construct=0, box=0, coincidences=0, energy=0; long long states=0; };
            std::vector<Phases> perThread(numThreads);
            const auto tick= []{ return std::chrono::steady_clock::now(); };
            const auto secondsSince= [](const std::chrono::steady_clock::time_point& t)
            { return std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count(); };
            const auto passOneBegan= tick();

            std::vector<Surveyed> surveyed;
            surveyed.reserve(subsets.size());
            long long visited= 0;
            const long long progressEvery=
                std::max<long long>(1, (long long)subsets.size()/200);
            std::ofstream out_file;

            #pragma omp parallel for num_threads(numThreads) schedule(dynamic) private(out_file)
            for (long long subsetIndex=0; subsetIndex < (long long)subsets.size(); ++subsetIndex)
            {
                const std::vector<int>& engaged= subsets[subsetIndex];

                if (!out_file.is_open())
                {
                    const std::string energyFileName= outputDirectory + "/output_thread_" +
                                                      std::to_string(omp_get_thread_num()) + ".txt";
                    out_file.open(energyFileName);
                    if (out_file.is_open())
                        out_file << "# one line per mesostate surveyed by this thread\n"
                                    "# nodes   = coincidences the deformed structure holds\n"
                                    "# engaged = nodes the enumeration engaged to reach it\n"
                                    "# fused   = atoms the overlap removal deletes\n"
                                    "# expelled = atoms the deformation carried out of their own\n"
                                    "#   grain, left out of the state entirely\n"
                                    "# unrelaxed / tethered / full = GB energy as constructed,\n"
                                    "#   after the tethered relaxation, and after the free one;\n"
                                    "#   the two relaxations are separate runs from the same start\n"
                                    "# engages = the ensemble indices this state engages\n"
                                    "# signature = the same thing in full, one 0/1 per ensemble\n"
                                    "#   member; this is what constructMesoState() takes, so a\n"
                                    "#   state can be rebuilt or inspected from this line alone\n"
                                    "# nodes  engaged  fused  expelled  corrugation  density"
                                    "  unrelaxed  tethered  spring  full  |  engages"
                                    "  |  signature\n";
                    else
                    {
#pragma omp critical (report)
                        std::cerr << "Failed to open file " << energyFileName << std::endl;
                    }
                }

                XTuplet state(ensembleSize);
                state.setZero();
                for (const int i : engaged) state(i)= 1;

                Phases& phase= perThread[omp_get_thread_num()];
                try {
                    auto mark= tick();
                    const auto& mesostate= ensemble.constructMesoState(state);
                    const auto [lowest,highest]= surfaceExtent(mesostate);
                    phase.construct+= secondsSince(mark);

                    // The scratch configuration: written because LAMMPS needs a file and the
                    // coincidence count needs the deformed positions, and then left to be
                    // overwritten by this thread's next state.
                    const std::string scratchBase= scratchDirectory + "/thread"
                                                 + std::to_string(omp_get_thread_num());
                    int expelled= 0, droppedCoincidences= 0;
                    mark= tick();
                    mesostate.box(scratchBase, &expelled, dropUnengagedCoincidences,
                                  &droppedCoincidences);
                    phase.box+= secondsSince(mark);
                    const std::string deformedFile= scratchBase + "_reference1.txt";

                    Surveyed record;
                    record.engaged    = engaged;
                    record.corrugation= highest-lowest;

                    mark= tick();
                    const Coincidences realized= countCoincidences(
                        configurationPositions(deformedFile),
                        mesostate.mesoStateCslVectors[1].cartesian(),
                        mesostate.mesoStateCslVectors[2].cartesian(),
                        lammpsOverlapCutoff);
                    phase.coincidences+= secondsSince(mark);
                    record.nodes= realized.sites;
                    record.fused= realized.fused;
                    record.expelled= expelled;

                    if (energiesRequested) {
                        // Both relaxations, each from this configuration; no dumps, since the
                        // structures themselves are not being kept in this pass.
                        mark= tick();
                        const auto relaxed= mesostate.relaxations(
                            lmpLocation, potentialName, deformedFile,
                            tetherHalfWidth, tetherStiffness, "", "", chainRelaxations);
                        phase.energy+= secondsSince(mark);
                        record.density  = relaxed.density;
                        record.unrelaxed= relaxed.unrelaxed;
                        record.tethered = relaxed.tethered;
                        record.spring   = relaxed.spring;
                        record.full     = relaxed.full;
                    }

                    std::ostringstream line;
                    line << record.nodes << "  " << record.engaged.size() << "  " << record.fused
                         << "  " << record.expelled
                         << "  " << std::fixed << std::setprecision(4) << record.corrugation
                         << "  " << std::setprecision(8) << record.density
                         << "  " << record.unrelaxed << "  " << record.tethered
                         << "  " << record.spring << "  " << record.full << "  | ";
                    for (const int i : record.engaged) line << " " << i;
                    line << "  |  " << fullSignature(record.engaged);

                    ++phase.states;
#pragma omp critical (report)
                    {
                        surveyed.push_back(std::move(record));
                        if (out_file.is_open()) out_file << line.str() << "\n";
                        // A line per state would be a million lines of screen for a sweep whose
                        // point is that nobody reads it state by state; pass 2 reports in full.
                        if (++visited % progressEvery == 0 || visited == (long long)subsets.size())
                            std::cout << "  surveyed " << visited << " / " << subsets.size()
                                      << "  (" << (100*visited)/(long long)subsets.size()
                                      << "%)" << std::endl;
                    }
                }
                // An exception that escapes a parallel region terminates the program, so
                // everything the construction can throw is caught here, not just runtime_error.
                catch (const std::exception& e) {
#pragma omp critical (report)
                    {
                        ++rejected;
                        reasons[std::string(e.what()).substr(0,70)]++;
                    }
                }
            }
            // The private copies of out_file are destroyed -- and so flushed and closed -- when
            // the parallel region ends.

            // ---- where pass 1 spent itself ------------------------------------------------
            {
                Phases total;
                for (const auto& p : perThread) {
                    total.construct+= p.construct; total.box+= p.box;
                    total.coincidences+= p.coincidences; total.energy+= p.energy;
                    total.states+= p.states;
                }
                const double accounted= total.construct+total.box+total.coincidences+total.energy;
                const double wall= secondsSince(passOneBegan);
                const std::string path= outputDirectory + "/profile.txt";
                std::ofstream report(path);
                const auto row= [&](const char* name, const double core, const char* note)
                {
                    report << "  " << std::left << std::setw(26) << name << std::right
                           << std::setw(12) << std::fixed << std::setprecision(1) << core
                           << std::setw(9) << std::setprecision(1) << (100*core/accounted) << "%"
                           << std::setw(11) << std::setprecision(2)
                           << (total.states ? 1000*core/total.states : 0.0)
                           << "   " << note << "\n";
                };
                report << "profile of pass 1\n"
                       << "=================\n\n"
                       << "  boundary            : " << boundaryDirectory << "\n"
                       << "  states surveyed     : " << total.states << "\n"
                       << "  threads             : " << numThreads << "\n"
                       << "  displacement field  : " << GbContinuum<3>::facetImageShells
                       << " image shell(s), refinement " << GbContinuum<3>::facetRefinement << "\n"
                       << "  energies            : " << (energiesRequested ? "on" : "off") << "\n\n"
                       << "  wall clock          : " << std::fixed << std::setprecision(1)
                       << wall << " s\n"
                       << "  core-seconds        : " << accounted << " s   ("
                       << std::setprecision(1) << (wall>0 ? accounted/wall : 0.0)
                       << "x the wall clock, against " << numThreads << " threads)\n\n"
                       << "  phase                     core-seconds     share    ms/state   what it is\n"
                       << "  ---------------------------------------------------------------------------\n";
                row("constructMesoState",  total.construct,
                    "triangulate the facets, build the continuum");
                row("box",                 total.box,
                    "evaluate the displacement field at every atom and write the file");
                row("countCoincidences",   total.coincidences,
                    "re-read that file and group atoms that met");
                row("relaxations (LAMMPS)",total.energy,
                    "both minimisations, in one invocation");
                report << "\n  Anything outside these four -- the walk, the manifest, the screen --\n"
                          "  is the difference between the core-seconds above and threads x wall.\n";
                report.close();
                std::cout << "  profile written to " << path << std::endl;
                std::cout << "    construct " << std::fixed << std::setprecision(1)
                          << (100*total.construct/accounted) << "%,  box "
                          << (100*total.box/accounted) << "%,  coincidences "
                          << (100*total.coincidences/accounted) << "%,  LAMMPS "
                          << (100*total.energy/accounted) << "%" << std::endl;
            }

            std::cout << "\nstates surveyed : " << surveyed.size() << std::endl;
            std::cout << "  rejected (clash / not triangulable) : " << rejected << std::endl;
            std::cout << "  candidate states examined : " << subsets.size() << std::endl;
            for (const auto& [message,count] : reasons)
                std::cout << "      " << count << " x  " << message << std::endl;
            std::cout << "  numbers in " << outputDirectory << "/output_thread_<id>.txt"
                      << std::endl;
            if (surveyed.empty())
                throw std::runtime_error("no state survived pass 1, so there is nothing to build.");

            // ============================================================== SHORTLIST
            // States engaging different numbers of coincidence points are not comparable -- a
            // boundary holding eight costs more than one holding a single point whatever it does
            // -- so the ranking is taken within each level of engagement separately, and the
            // level is the count the structure realised, not the count the enumeration asked for.
            //
            // Within a level the enumeration produces whole families of states that differ only
            // by where along the periodic directions their coincidences sit, and these are
            // degenerate: every number recorded for them agrees.  Left alone one family fills the
            // whole top `shortlistTop` with copies of a single structure and hides the next
            // distinct one, so they are folded together and the survivor carries the count.
            // Families are found once over the whole sweep rather than per level, so that the
            // per-level and the global selections agree about what counts as a duplicate.
            std::vector<const Surveyed*> everyState;
            everyState.reserve(surveyed.size());
            for (const auto& record : surveyed) everyState.push_back(&record);
            std::vector<int> familyOf;
            std::map<const Surveyed*,int> familyIndex;
            if (shortlistDedup) {
                familyOf= groupByFingerprint(everyState, shortlistDedupTol);
                for (std::size_t i=0; i<everyState.size(); ++i)
                    familyIndex[everyState[i]]= familyOf[i];
            }
            else {
                // No dedup: every state is its own family, so nothing is ever folded.
                for (std::size_t i=0; i<everyState.size(); ++i)
                    familyIndex[everyState[i]]= (int)i;
            }

            std::map<int,std::vector<const Surveyed*>> byLevel;
            for (const auto& record : surveyed) byLevel[record.nodes].push_back(&record);

            struct Shortlisted
            {
                const Surveyed* record= nullptr;
                int level= 0;
                int rank= 0;
                int copies= 1;
                //! Why this state is being built: its place in the per-level shortlist, and its
                //! place in whichever global rankings reached it.  A state can be picked by
                //! several of them, and it is built once.
                std::string why;
            };
            std::vector<Shortlisted> shortlist;
            // Which shortlist entry a family already occupies, so that a state reached by a
            // second criterion is annotated rather than built twice.
            std::map<int,int> slotOfFamily;

            std::cout << "\nshortlisting the lowest " << shortlistTop << " state(s) per level by "
                      << rankName(shortlistRankBy)
                      << (shortlistDedup ? ", duplicates folded" : ", no dedup") << std::endl;
            for (auto& [level,members] : byLevel)
            {
                std::stable_sort(members.begin(), members.end(),
                                 [&](const Surveyed* a, const Surveyed* b)
                                 { return a->rankValue(shortlistRankBy)
                                        < b->rankValue(shortlistRankBy); });

                // The whole level is walked even once `shortlistTop` distinct states are in hand:
                // the copies of a kept state are spread through the ranking, and counting them is
                // what turns "20 of 258" into a statement about how large each family is.  A
                // family that begins past the cut is given a slot of -1, so that its own copies
                // are not counted against a state they are not copies of.
                std::map<int,int> seenHere;         // family -> its slot, or -1 past the cut
                const std::size_t levelBegin= shortlist.size();
                long long folded= 0;
                for (const Surveyed* member : members)
                {
                    const int family= familyIndex.at(member);
                    const auto found= seenHere.find(family);
                    if (found != seenHere.end()) {
                        if (found->second >= 0) ++shortlist[found->second].copies;
                        ++folded;
                        continue;
                    }
                    if ((int)(shortlist.size()-levelBegin) >= shortlistTop) {
                        seenHere.emplace(family, -1);
                        continue;
                    }
                    seenHere.emplace(family, (int)shortlist.size());
                    const int rank= (int)(shortlist.size()-levelBegin)+1;
                    std::ostringstream why;
                    why << "L" << level << "#" << rank;
                    slotOfFamily[family]= (int)shortlist.size();
                    shortlist.push_back({member, level, rank, 1, why.str()});
                }
                const std::size_t kept= shortlist.size()-levelBegin;
                std::cout << "  " << std::setw(2) << level << " node(s): "
                          << std::setw(8) << members.size() << " state(s) -> kept "
                          << std::setw(2) << kept;
                if (kept > 0)
                    std::cout << ",  " << rankName(shortlistRankBy) << " "
                              << std::fixed << std::setprecision(6)
                              << shortlist[levelBegin].record->rankValue(shortlistRankBy)
                              << " .. "
                              << shortlist.back().record->rankValue(shortlistRankBy);
                if (shortlistDedup) std::cout << ",  " << folded << " copy(s) folded in";
                std::cout << std::endl;
            }

            // ---- and the globally cheapest, by each energy in turn ------------------------
            // The per-level shortlist above never compares one level with another, so a state
            // that is cheap for the whole sweep can go unbuilt because its own level held forty
            // cheaper ones.  These selections close that gap.  A state already chosen -- by a
            // level, or by an earlier energy -- is annotated rather than built again, so the
            // three sets overlap without costing anything.
            if (globalTop > 0)
            {
                std::cout << "\nadding the globally cheapest " << globalTop
                          << " state(s) by each energy" << std::endl;
                for (const auto& [by,tag] : {std::pair<RankBy,const char*>{RankBy::Unrelaxed,"U"},
                                             {RankBy::Tethered,"T"},
                                             {RankBy::Full,"F"}})
                {
                    std::vector<const Surveyed*> ordered= everyState;
                    std::stable_sort(ordered.begin(), ordered.end(),
                                     [by](const Surveyed* a, const Surveyed* b)
                                     { return a->rankValue(by) < b->rankValue(by); });
                    std::set<int> takenHere;
                    int kept= 0, added= 0;
                    double lowest= 0.0, highest= 0.0;
                    for (const Surveyed* record : ordered)
                    {
                        if (kept >= globalTop) break;
                        const int family= familyIndex.at(record);
                        if (!takenHere.insert(family).second) continue;   // a copy of one already
                        if (kept==0) lowest= record->rankValue(by);
                        highest= record->rankValue(by);
                        ++kept;
                        std::ostringstream why;
                        why << tag << "#" << kept;
                        const auto found= slotOfFamily.find(family);
                        if (found != slotOfFamily.end()) {
                            shortlist[found->second].why += "," + why.str();
                            continue;                                    // already being built
                        }
                        slotOfFamily[family]= (int)shortlist.size();
                        shortlist.push_back({record, record->nodes, 0, 1, why.str()});
                        ++added;
                    }
                    std::cout << "  " << std::setw(9) << rankName(by) << " : " << kept
                              << " distinct state(s), " << added << " newly built,  "
                              << std::fixed << std::setprecision(6)
                              << lowest << " .. " << highest << " eV" << std::endl;
                }
                std::cout << "  shortlist now holds " << shortlist.size() << " state(s)" << std::endl;
            }

            // ---- every state, ordered by each energy in turn ------------------------------
            // Pass 2 builds only the shortlist, so without these the other twelve thousand
            // states exist only as unordered lines spread across the per-thread files.  The
            // three orderings are written separately rather than one file sorted three ways
            // because the three energies disagree about which states are good -- a state can
            // be cheap as constructed and unremarkable once relaxed, or the reverse -- and
            // which ordering is the right one is a question about what is being asked.
            //
            // Each row carries the signature twice over: as the indices it engages, which is
            // what can be read, and in full, which is what constructMesoState() takes.  The
            // length is written alongside because a signature means nothing without it -- the
            // same list of engaged indices describes a different state in an ensemble of a
            // different size.
            {
                const auto writeRanking = [&](const std::string& name, const RankBy by)
                {
                    std::vector<const Surveyed*> ordered;
                    ordered.reserve(surveyed.size());
                    for (const auto& record : surveyed) ordered.push_back(&record);
                    std::stable_sort(ordered.begin(), ordered.end(),
                                     [by](const Surveyed* a, const Surveyed* b)
                                     { return a->rankValue(by) < b->rankValue(by); });

                    const std::string path= outputDirectory + "/sortedBy" + name + ".txt";
                    std::ofstream file(path);
                    file << "# every state of the sweep, ordered by " << rankName(by) << "\n"
                            "# nodes   = coincidences the deformed structure holds\n"
                            "# engaged = nodes the enumeration engaged to reach it\n"
                            "# fused   = atoms the overlap removal deletes\n"
                            "# unrelaxed / tethered / full = GB energy as constructed, after the\n"
                            "#   tethered relaxation, and after the free one, in eV\n"
                            "# engages = the ensemble indices this state engages\n"
                            "# signature = the same thing in full, one 0/1 per ensemble member;\n"
                            "#   this is what constructMesoState() takes\n"
                            "# rank  engaged  signatureLength  nodes  fused  corrugation  density"
                            "  unrelaxed  tethered  spring  full  |  engages  |  signature\n";
                    int rank= 0;
                    for (const Surveyed* record : ordered) {
                        file << ++rank
                             << "  " << record->engaged.size()
                             << "  " << ensembleSize
                             << "  " << record->nodes
                             << "  " << record->fused
                             << "  " << std::fixed << std::setprecision(4) << record->corrugation
                             << "  " << std::setprecision(8) << record->density
                             << "  " << record->unrelaxed
                             << "  " << record->tethered
                             << "  " << record->spring
                             << "  " << record->full
                             << "  | ";
                        for (const int i : record->engaged) file << " " << i;
                        file << "  |  " << fullSignature(record->engaged) << "\n";
                    }
                    std::cout << "  " << path << std::endl;
                };
                std::cout << "\nordering all " << surveyed.size() << " state(s) by each energy:"
                          << std::endl;
                writeRanking("Unrelaxed", RankBy::Unrelaxed);
                writeRanking("Tethered",  RankBy::Tethered);
                writeRanking("Full",      RankBy::Full);
            }

            // ================================================================ PASS 2
            // Build the shortlisted states for real: their configurations, both relaxed
            // structures, and a manifest.  This is the only output anyone opens, and it is a few
            // hundred states rather than a million.
            std::cout << "\npass 2 of 2 : building " << shortlist.size()
                      << " shortlisted state(s)" << std::endl;

            std::ofstream manifest(outputDirectory + "/states.txt");
            manifest << "# state_<index>_0.txt = undeformed, state_<index>_1.txt = deformed,\n"
                        "# dump.state_<index>_2 = tethered relaxation, _3 = free relaxation\n"
                        "# nodes   = coincidences the deformed structure holds\n"
                        "# engaged = nodes the enumeration engaged to reach it\n"
                        "# fused   = atoms the relaxation deletes (one per atom past the first\n"
                        "#           at each coincidence)\n"
                        "# expelled= atoms the deformation carried out of their own grain, and\n"
                        "#           which the state therefore does not contain at all\n"
                        "# copies  = states in this one's family, all recorded identically\n"
                        "# The displacement field cannot be aimed at the engaged nodes alone, so\n"
                        "# it closes other pairs of atoms with them and nodes >= engaged.  Those\n"
                        "# extra coincidences are as real as the engaged ones, so `nodes` is what\n"
                        "# describes the boundary and `engaged` only how it was reached.\n"
                        "# unrelaxed / tethered / full = GB energy as constructed, after the\n"
                        "#   tethered relaxation, and after the free one; the two relaxations are\n"
                        "#   separate runs from the same starting configuration.\n"
                        "# selectedBy = why this state was built: L<level>#<rank> for its place\n"
                        "#   in that level's shortlist, and U#/T#/F# for its place in the ranking\n"
                        "#   over the whole sweep by unrelaxed / tethered / full energy\n"
                     << "# index  nodes  engaged  fused  expelled  copies  selectedBy  corrugation"
                        "  density  unrelaxed  tethered  spring  full  engaged nodes\n";

            // The shortlisted signatures on their own, in the layout readPostProcessingOutput()
            // parses: three leading columns and then the signature, which is where that reader
            // starts.  Pointing `fin` at this file with enumerateStates=false rebuilds exactly
            // these states, so the shortlist can be revisited without surveying the sweep again.
            // They are kept apart from states.txt because 540 columns in the middle of it would
            // sit between the named columns and the free-text node descriptions that follow, and
            // nothing that reads states.txt expects that.
            std::ofstream signatures(outputDirectory + "/signatures.txt");
            signatures << "# the shortlisted states' signatures, in the layout `fin` is read in:\n"
                          "# index  density  nodes  then one 0/1 per ensemble member.\n"
                          "# Every column is an integer, since the reader stops at the first token\n"
                          "# that is not one; the index is the state_<index>_*.txt number without\n"
                          "# its leading zeros, and the density is the atom count.\n"
                          "# Set enumerateStates=false and fin to this file to rebuild them.\n";

            int built=0, unbuilt=0;
            std::map<std::string,int> buildFailures;
            // Pass 2 reports every state in full, but a few hundred blocks scroll past without
            // saying how far along they are, so each carries its own place in the run.
            long long constructed= 0;

#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
            for (long long entry=0; entry < (long long)shortlist.size(); ++entry)
            {
                const Shortlisted& chosen= shortlist[entry];
                const std::string index= indexName((int)entry);

                XTuplet state(ensembleSize);
                state.setZero();
                for (const int i : chosen.record->engaged) state(i)= 1;

                try {
                    const auto& mesostate= ensemble.constructMesoState(state);
                    const auto [lowest,highest]= surfaceExtent(mesostate);

                    // box() writes <name>_reference0.txt and _reference1.txt; rename them to
                    // state_<index>_<config>.txt, config 0 undeformed and 1 deformed.
                    const std::string base= scratchDirectory + "/build" + index;
                    int expelled= 0, droppedCoincidences= 0;
                    mesostate.box(base, &expelled, dropUnengagedCoincidences,
                                  &droppedCoincidences);
                    for (const int configuration : {0,1})
                        std::filesystem::rename(
                            base + "_reference" + std::to_string(configuration) + ".txt",
                            outputDirectory + "/state_" + index + "_"
                                            + std::to_string(configuration) + ".txt");

                    // Every site the enumeration considered, marked with whether this state
                    // took it.  Only the undeformed configuration carries them, so the deformed
                    // one that LAMMPS reads is untouched.
                    {
                        std::vector<char> engagedSite(distinctSites.size(), 0);
                        for (const int i : chosen.record->engaged)
                            if (siteOfIndex[i] >= 0) engagedSite[siteOfIndex[i]]= 1;
                        appendSites(outputDirectory + "/state_" + index + "_0.txt",
                                    distinctSites, engagedSite, siteMarkerRadius);
                    }

                    const std::string deformedFile= outputDirectory + "/state_" + index + "_1.txt";
                    const Coincidences realized= countCoincidences(
                        configurationPositions(deformedFile),
                        mesostate.mesoStateCslVectors[1].cartesian(),
                        mesostate.mesoStateCslVectors[2].cartesian(),
                        lammpsOverlapCutoff);

                    // Both relaxed structures are kept this time, each in its own dump: _2 is
                    // the tethered one, _3 the free one.
                    GbMesoState<3>::Relaxations relaxed;
                    if (energiesRequested)
                        relaxed= mesostate.relaxations(
                            lmpLocation, potentialName, deformedFile,
                            tetherHalfWidth, tetherStiffness,
                            std::filesystem::absolute(outputDirectory + "/dump.state_" + index
                                                      + "_2").string(),
                            std::filesystem::absolute(outputDirectory + "/dump.state_" + index
                                                      + "_3").string(),
                            chainRelaxations);

                    std::ostringstream report;
                    report << "  [" << index << "] " << chosen.level << " node(s), "
                           << chosen.record->engaged.size() << " engaged, "
                           << chosen.copies << " copy(s), chosen as " << chosen.why
                           << "   -- built %PROGRESS%"
                           << "\n           GB signature: " << state;
                    for (std::size_t e=0; e<chosen.record->engaged.size(); ++e)
                        report << "\n           " << describe(mesostate,e);
                    report << "\n           surface x.n in [" << std::fixed << std::setprecision(4)
                           << lowest << ", " << highest << "] A,  corrugation = "
                           << highest-lowest << " A"
                           << "\n           coincidences: " << realized.sites << " site(s), "
                           << realized.fused << " atom(s) fused, "
                           << realized.largest << " at the most crowded";
                    if (expelled > 0)
                        report << "\n           " << expelled
                               << " atom(s) expelled: the deformation carried them out of their "
                                  "own grain";
                    if (droppedCoincidences > 0)
                        report << "\n           " << droppedCoincidences
                               << " atom(s) dropped: they had drifted into coincidences this "
                                  "state did not engage";
                    if (realized.sites < (int)chosen.record->engaged.size())
                        report << "\n           WARNING: " << chosen.record->engaged.size()
                               << " node(s) engaged but only " << realized.sites << " realised";
                    if (energiesRequested)
                        report << "\n           density = " << std::setprecision(6)
                               << relaxed.density
                               << "   unrelaxed = " << relaxed.unrelaxed
                               << "   tethered = " << relaxed.tethered
                               << " (spring " << relaxed.spring << ")"
                               << "   full = " << relaxed.full;

                    std::ostringstream manifestLine;
                    manifestLine << index
                                 << "  " << realized.sites
                                 << "  " << chosen.record->engaged.size()
                                 << "  " << realized.fused
                                 << "  " << expelled
                                 << "  " << chosen.copies
                                 << "  " << (chosen.why.empty() ? std::string("-") : chosen.why)
                                 << "  " << std::fixed << std::setprecision(4) << highest-lowest
                                 << "  " << std::setprecision(8) << relaxed.density
                                 << "  " << relaxed.unrelaxed
                                 << "  " << relaxed.tethered
                                 << "  " << relaxed.spring
                                 << "  " << relaxed.full;
                    for (std::size_t e=0; e<chosen.record->engaged.size(); ++e)
                        manifestLine << "   " << describe(mesostate,e);

                    // Every column here is an integer, because readPostProcessingOutput() parses
                    // integers and stops at the first token that is not one: a density written as
                    // 270.00000000 ends the row at the decimal point and the signature never gets
                    // read.  The density is an atom count, so nothing is lost by rounding it.
                    std::ostringstream signatureLine;
                    signatureLine << (int)entry << "  "
                                  << (long long)std::llround(relaxed.density) << "  "
                                  << realized.sites << "  "
                                  << fullSignature(chosen.record->engaged);

#pragma omp critical (report)
                    {
                        ++built;
                        // Filled in here rather than above, because the count is only meaningful
                        // once the critical section has claimed this state's place in the order.
                        std::ostringstream progress;
                        progress << ++constructed << " / " << shortlist.size() << "  ("
                                 << (100*constructed)/(long long)shortlist.size() << "%)";
                        std::string line= report.str();
                        const std::size_t marker= line.find("%PROGRESS%");
                        if (marker != std::string::npos)
                            line.replace(marker, std::strlen("%PROGRESS%"), progress.str());
                        std::cout << line << std::endl;
                        manifest << manifestLine.str() << "\n";
                        signatures << signatureLine.str() << "\n";
                    }
                }
                catch (const std::exception& e) {
#pragma omp critical (report)
                    {
                        ++unbuilt;
                        buildFailures[std::string(e.what()).substr(0,70)]++;
                    }
                }
            }
            manifest.close();
            signatures.close();

            // The scratch configurations have served their purpose; leaving a directory of
            // per-thread leftovers behind would be the small version of the problem this pass
            // exists to solve.
            std::error_code removalFailed;
            std::filesystem::remove_all(scratchDirectory, removalFailed);

            std::cout << "\nshortlisted states built : " << built << " of " << shortlist.size()
                      << std::endl;
            if (unbuilt > 0) {
                std::cout << "  failed to rebuild : " << unbuilt << std::endl;
                for (const auto& [message,count] : buildFailures)
                    std::cout << "      " << count << " x  " << message << std::endl;
            }
            std::cout << "  configurations and manifest in "
                      << std::filesystem::absolute(outputDirectory).string() << std::endl;
            std::cout << "  every state's numbers and full signature in " << outputDirectory
                      << "/output_thread_<id>.txt" << std::endl;
            std::cout << "  the shortlisted signatures, ready to be read back as `fin`, in "
                      << outputDirectory << "/signatures.txt" << std::endl;
            return 0;
        }

		std::cout << "Size of input mesostate = " << stateSize << std::endl;

        // A state vector of the wrong length must not be allowed through.  getEngagedTsPairs()
        // only asserts this, and asserts are compiled out of a Release build, so a longer state
        // is silently truncated to the ensemble's length: the run then quietly describes a
        // different mesostate than the input file asked for.
        if (stateSize != ensembleSize)
            throw std::runtime_error("The mesostate signatures in " + fin + " have " +
                std::to_string(stateSize) + " entries but this ensemble has " +
                std::to_string(ensembleSize) + " constraints. Each row of the file must hold "
                "3 leading columns (count, rho, and one unused) followed by exactly " +
                std::to_string(ensembleSize) + " zeros/ones.");

        XTuplet state(stateSize);

        // Same layout as tests/testGbMesoState: one line per mesostate, holding the GB signature
        // and the (density, energy) LAMMPS reports for it.
        std::ofstream energyFile;
        if (energiesRequested) {
            energyFile.open("energies.txt");
            energyFile << "# GB signature  density  energy"
                       << (minimizeInLammps ? "  (minimized)" : "  (unrelaxed)") << "\n";
        }

        for(int state_count = 0; state_count<rows;state_count++)
        {
            // Define state
            int start = 3;
            int length = stateData[state_count].size()-start;
            if (state.size() != length)
                throw std::runtime_error("Row " + std::to_string(state_count+1) + " of " + fin +
                    " has " + std::to_string(length) + " state entries; expected " +
                    std::to_string(state.size()) + ".");
            for (int k=0; k<length; ++k)
                state(k)= stateData[state_count][k+3];

            std::cout <<"rho = "<< stateData[state_count][1] << " count = "<< stateData[state_count][0]
                      << "\n           GB signature: " << state << std::endl;
            const auto& mesostate= ensemble.constructMesoState(state);

            // What this signature engages, printed next to the signature itself.  Which
            // container holds it depends on the search: Sites builds coincidence nodes and
            // leaves tShiftPairs empty, while Flat and Full do the reverse.  Reading the wrong
            // one is not a wrong number but an out-of-range access -- this path used to assume
            // (t,s) pairs and abort on a Sites run, which is exactly the run that produces the
            // signatures.txt this branch is most useful for reading back.
            for (int i=0; i<stateSize; ++i) {
                if (state(i) != 1) continue;
                if (searchMode==GbShiftSearch::Sites) {
                    const auto& node= ensemble.nodes[i];
                    std::cout << "           s=(" << std::fixed << std::setprecision(4)
                              << node.s(0) << "," << node.s(1) << "," << node.s(2)
                              << ")  |uA|=" << node.uA().norm()
                              << "  |uB|=" << node.uB().norm()
                              << "  |t|=" << (node.uA()-node.uB()).norm() << std::endl;
                }
                else {
                    const Eigen::Vector3d t= ensemble.tShiftPairs[i].first.cartesian();
                    const Eigen::Vector3d sh= ensemble.tShiftPairs[i].second;
                    std::cout << "           t=(" << std::fixed << std::setprecision(4)
                              << t(0) << "," << t(1) << "," << t(2) << ")  s=("
                              << sh(0) << "," << sh(1) << "," << sh(2) << ")" << std::endl;
                }
            }

            // Output the file for mesostate configuration
            std::string preamble = "mesostate_data_";
            std::string spacer = "_";
            std::string rho = std::to_string(stateData[state_count][1]);
            std::string count = std::to_string(stateData[state_count][0]);
            std::string outfilename = preamble + rho + spacer + count;
            std::cout << outfilename<< std::endl;
            mesostate.box(outfilename);
            std::cout << "Done writing mesostate file: " << outfilename << std::endl;
            mesostate.exportMesh(outfilename);

            if (energiesRequested) {
                const auto [density,gbEnergy]=
                    mesostate.densityEnergy(lmpLocation, potentialName, minimizeInLammps);
                std::cout << "           density = " << std::fixed << std::setprecision(6)
                          << density << "   energy = " << gbEnergy
                          << (minimizeInLammps ? "  (minimized)" : "  (unrelaxed)") << std::endl;
                energyFile << state << "  " << std::setprecision(8) << density
                           << "  " << gbEnergy << std::endl;
            }
        }
  }
  catch(std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}