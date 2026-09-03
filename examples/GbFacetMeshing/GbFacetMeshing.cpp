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

int main()
{
    /*! [Types] */
    using VectorDimI = LatticeCore<3>::VectorDimI;
    using VectorDimD = LatticeCore<3>::VectorDimD;
    using Vector2d= Eigen::Vector2d;
    using IntScalarType = LatticeCore<3>::IntScalarType;


    /*/ Sigma 123 [110](-5 5 14)
    VectorDimD axis(1,1,0);
    double theta= 53.594515175286005615*std::numbers::pi/180;       // misorientation angle
    VectorDimD gbNormal(-5,5,14);                        // Miller indices
    int heightScaling= 2;
    int periodScaling= 1;
    int axisScaling= 1;
    double bScaling= 0.55;
	*/
	// Sigma 5[001]()
	VectorDimD axis(0, 0, 1);
	double theta= 53.130102354155994249*std::numbers::pi/180;       // misorientation angle
	VectorDimD gbNormal(-2,1,0);                        // Miller indices
	int heightScaling= 3;
	int periodScaling= 1;
	int axisScaling= 1;

    // ------------------------------------------------------------------ USER CHOICE
    // Every setting of the run lives in this block.  Edit it and rebuild; nothing is read from
    // the command line.
    //
    // enumerateStates: build every mesostate the ensemble's (t,s) pairs admit.  Set it false to
    // read ready-made signatures from `fin` instead.
    const bool enumerateStates            = true;
    // Flat = the original enumeration (CSL shifts along the boundary); Full = ball + slab;
    // Sites = coincidence points first, with the two grains displaced independently.
    const GbShiftSearch searchMode        = GbShiftSearch::Sites;
    // Sites only.  The slab about the flat boundary that coincidence points are taken from, and
    // how far one atom may move to reach one.  The points fall on planes a fixed distance apart
    // (0.404 A for this boundary), so a slab thinner than that spacing admits only the boundary
    // plane itself and the run is flat whatever the search.
    const double slabHalfThickness        = 0.5;   // Angstrom
    const double dMax                     = 1.5;   // Angstrom
    // Two engaged coincidence points must be at least this far apart in projection onto the
    // boundary plane.  Zero admits any pair that is not exactly coincident.  The points lie on a
    // lattice with spacings 0.404, 0.808, 1.617 A here, so the criterion bites in steps: 0.45
    // removes the 0.404 A pairs -- the ones that drive ordinary, non-node atoms of the two grains
    // to within half an Angstrom of each other -- and keeps everything at 0.808 A and beyond.
    const double minSiteSeparation        = 0.45;  // Angstrom
    // FLAT STGB (commented out).  Restricting to the shifts that lie in the boundary plane
    // leaves the deformed surface planar, which is the flat symmetric-tilt suite.  This run is
    // the general one instead: every (t,s) pair the ensemble holds is available and the boundary
    // is free to facet.  To get the flat suite back, restore this line and the two blocks marked
    // "FLAT STGB" further down.
    //const bool restrictToFlatStates     = true;
    // How many (t,s) pairs one state may engage at once.  0 means no limit, in which case the
    // only ceiling is the clash rule itself: no state can engage two pairs that share a lattice
    // site, so no state can be larger than the smaller of the two site counts.  Every state the
    // walk produces is clash-free, so this is a count of real states, not of candidates.
    const int maxEngaged                  = 8;
    // Refuse to start a run longer than this many states.  Each one is a mesostate construction
    // and, with energies on, a LAMMPS minimization.
    const long long maxStates             = 150000;
    const double tMax                     = 0.99;   // ball radius, in units of b
    const double tPerpMax                 = 0.99;   // slab half-thickness, in units of b
    // The translations that give well spread out flat sites point along the GB normal, with
    // |t.n| = 1.617 A here, so a slab of 0.5 b = 1.278 A cuts them out and the sites collapse onto
    // two lattice points.  1.0 b = 2.556 A clears them.
    const double layerHalfThickness       = 1.0;  // Angstrom either side of the flat GB
    const std::string fin                 = "sigma5.txt";

    // ---- LAMMPS ----------------------------------------------------------------
    // Same setup as tests/testGbMesoState: each mesostate is written out, handed to LAMMPS and
    // its (density, GB energy) read back.  Set computeEnergies false to only write the
    // configurations, and lmpLocation to the serial LAMMPS executable on this machine.
    const bool computeEnergies            = true;
    const bool minimizeInLammps           = true;  // relax before reading the energy
    // Threads that build and evaluate mesostates in parallel, one LAMMPS process each.  Same
    // shape as tests/testGbMesoState, which runs with num_threads(1); raise it once a serial
    // pass has been seen to work.
    const int numThreads                  = 80;
    const std::string potentialName       = "Cu_mishin1.eam.alloy";
    const std::string lmpLocation         = "/usr/bin/lmp";
    // -------------------------------------------------------------------------------

    std::cout << "states = general (every shift, faceted boundaries allowed)" << std::endl;
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
        cslVectors.push_back(heightScaling*gb.bc.csl.latticeDirection(gb.nA.cartesian()).latticeVector());
        cslVectors.push_back(periodScaling*gb.getPeriodVector(rAxisA));
        cslVectors.push_back(axisScaling*axisC);
        gb.box(cslVectors,1,"gb.txt",false);
        bc.box(cslVectors,1,"bcOriented.txt",true);


        // material parameter
        // source: https://openkim.org/id/1999--Mishin-Y-Farkas-D-Mehl-M-J-Papaconstantopoulos-D-A--Al
        double c11= 113.796/160.2176621;
        double c12= 61.55/160.2176621;
        //GbMaterialTensors::lambda= c12;
        //GbMaterialTensors::mu= (c11-c12)/2;

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

            // The general run takes every (t,s) pair the ensemble holds.
            std::vector<int> family;
            for (int i=0; i<ensembleSize; ++i)
                family.push_back(i);

            // FLAT STGB : keeping only the shifts that lie exactly in the
            // boundary plane leaves the deformed surface planar, so every signature over them is
            // a flat STGB -- the baseline the faceted states depart from.
            //std::vector<int> family;
            //for (int i=0; i<ensembleSize; ++i)
            //    if (std::abs(ensemble.tShiftPairs[i].second.dot(nHat)) < 1.0e-9)
            //        family.push_back(i);
            //if (family.empty())
            //    throw std::runtime_error("No shift lies in the boundary plane, so there is no "
            //                             "flat GB to build. Raise tMax or tPerpMax.");

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
            const std::string outputDirectory=
                std::string("generalGB_")
                + (searchMode==GbShiftSearch::Flat ? "flat" : "full");
            const std::string meshDirectory = outputDirectory + "/mesh";
            for (const auto& directory : {outputDirectory, meshDirectory})
                std::filesystem::create_directories(directory);
            std::ofstream manifest(outputDirectory + "/states.txt");
            manifest << "# state_<index>_0.txt = undeformed, state_<index>_1.txt = deformed\n"
                     << "# index  sites  relief  density  energy  engaged (t,s) pairs\n";

            std::cout << "writing to " << std::filesystem::absolute(outputDirectory).string()
                      << std::endl;
            std::cout << "threads = " << numThreads << std::endl;

            auto indexName= [](int i){
                std::ostringstream o; o << std::setw(3) << std::setfill('0') << i; return o.str(); };

            // One line describing engaged candidate i, whichever search produced it.  A Sites
            // candidate is a coincidence point with the two grains displaced independently, so
            // it is reported as the point and the two displacements; the other searches split
            // one translation evenly, so they are reported as (t,s) as before.
            const auto describe= [&ensemble,&searchMode,&nHat](const int i)
            {
                std::ostringstream o;
                o << std::fixed << std::setprecision(4);
                if (searchMode==GbShiftSearch::Sites) {
                    const auto& node= ensemble.nodes[i];
                    o << "s=(" << node.s(0) << "," << node.s(1) << "," << node.s(2) << ")"
                      << "  s.n=" << node.s.dot(nHat)
                      << "  |uA|=" << node.uA().norm() << "  |uB|=" << node.uB().norm()
                      << "  |t|=" << node.t().norm();
                }
                else {
                    const Eigen::Vector3d t= ensemble.tShiftPairs[i].first.cartesian();
                    const Eigen::Vector3d sh= ensemble.tShiftPairs[i].second;
                    o << "t=(" << t(0) << "," << t(1) << "," << t(2) << ")  s=("
                      << sh(0) << "," << sh(1) << "," << sh(2) << ")  s.n=" << sh.dot(nHat);
                }
                return o.str();
            };

            int built=0, rejected=0;
            std::map<std::string,int> reasons;

            // The energies, written the way tests/testGbMesoState writes them: one line per
            // mesostate holding the GB signature, the density and the energy, in a file named
            // after the thread that produced it.  The stream is private to each thread and
            // opened on that thread's first state, so the threads never share a file handle.
            std::ofstream out_file;

            // Every thread runs its own LAMMPS process.  energy() and box() both name their
            // scratch files after omp_get_thread_num() -- in<id>.find_energy,
            // data<id>.lammps_input, temp<id>_reference1.txt -- so the processes do not read
            // each other's input.  Everything shared between threads (the state counter, the
            // manifest, the reason tally, the screen) is written inside a critical section.
#pragma omp parallel for num_threads(numThreads) schedule(dynamic) private(out_file)
            for (long long subsetIndex=0; subsetIndex < (long long)subsets.size(); ++subsetIndex)
            {
                const std::vector<int>& engaged= subsets[subsetIndex];

                if (energiesRequested && !out_file.is_open())
                {
                    const std::string energyFileName= outputDirectory + "/output_thread_" +
                                                      std::to_string(omp_get_thread_num()) + ".txt";
                    out_file.open(energyFileName);
                    if (out_file.is_open())
                        out_file << "# GB signature  density  energy"
                                 << (minimizeInLammps ? "  (minimized)" : "  (unrelaxed)") << "\n";
                    else
                    {
#pragma omp critical (report)
                        std::cerr << "Failed to open file " << energyFileName << std::endl;
                    }
                }

                XTuplet state(ensembleSize);
                state.setZero();
                for (const int i : engaged) state(i)= 1;

                try {
                    const auto& mesostate= ensemble.constructMesoState(state);

                    // How far the deformed surface departs from the flat boundary.  This has to
                    // be asked of the deformed surface: the reference facets are stepped,
                    // because t varies between sites.  In the general run it is reported, not
                    // enforced -- a non-zero value is the faceting itself.
                    const Eigen::MatrixXd deformed= mesostate.facetA.deformedVertices();
                    double outOfPlane= 0.0;
                    for (int r=0; r<deformed.rows(); ++r)
                        outOfPlane= std::max(outOfPlane,
                            std::abs(Eigen::Vector3d(deformed.row(r)).dot(nHat)));

                    // FLAT STGB (commented out): reject anything that is not planar, so that
                    // only the flat symmetric-tilt states survive.
                    //if (outOfPlane > 1.0e-6)
                    //    throw std::runtime_error("state is not flat, |s.n| up to " +
                    //                             std::to_string(outOfPlane));

                    // Claim this state's index.  With several threads the subsets finish out of
                    // order, so the index records the order they were built in, not the mask --
                    // which pairs an index holds is what states.txt is for.
                    int index;
#pragma omp critical (stateIndex)
                    index= built++;

                    const std::string name= meshDirectory + "/state_" + indexName(index);
                    mesostate.box(name);
                    //mesostate.exportMesh(name);

                    // box() writes <name>_reference0.txt and _reference1.txt; rename them to
                    // state_<state>_<config>.txt, config 0 undeformed and 1 deformed.
                    for (const int configuration : {0,1})
                        std::filesystem::rename(
                            name + "_reference" + std::to_string(configuration) + ".txt",
                            outputDirectory + "/state_" + indexName(index) + "_"
                                            + std::to_string(configuration) + ".txt");

                    // Hand the deformed configuration to LAMMPS and read back what it costs.
                    // The file renamed just above is the one LAMMPS wants, so it is passed over
                    // rather than regenerated: writing a configuration evaluates the displacement
                    // field at every atom and is by far the most expensive step of a state, so
                    // letting densityEnergy() write its own scratch copy doubled the cost of the
                    // whole run.  LAMMPS leaves its relaxed structure in dump.state_<index>_2,
                    // beside the undeformed and deformed configurations.
                    double density= 0.0, gbEnergy= 0.0;
                    if (energiesRequested) {
                        const std::string deformedFile= outputDirectory + "/state_"
                                                      + indexName(index) + "_1.txt";
                        const std::string minimizedDump= std::filesystem::absolute(
                            outputDirectory + "/dump.state_" + indexName(index) + "_2").string();
                        std::tie(density,gbEnergy)=
                            mesostate.densityEnergy(lmpLocation, potentialName, minimizeInLammps,
                                                    deformedFile, minimizedDump);
                        if (out_file.is_open())
                            out_file << state << "  " << std::setprecision(8) << density
                                     << "  " << gbEnergy << std::endl;
                    }

                    // Build the whole screen block first and print it in one go: with several
                    // threads reporting, a line-by-line print interleaves them into nonsense.
                    // The pair ordering differs between the two search rules, so print the (t,s)
                    // values themselves -- those are what can be compared across them -- with
                    // the signature, which is what identifies the mesostate to the rest of the
                    // ensemble and what the energy is recorded against.
                    std::ostringstream report;
                    report << "  [" << index << "] thread " << omp_get_thread_num() << ", "
                           << engaged.size() << " site(s) -> " << name
                           << "\n           GB signature: " << state;
                    for (const int i : engaged)
                        report << "\n           " << describe(i);
                    report << "\n           facet relief |x.n| = " << std::fixed
                           << std::setprecision(4) << outOfPlane << " A";
                    if (energiesRequested)
                        report << "\n           density = " << std::fixed << std::setprecision(6)
                               << density << "   energy = " << gbEnergy
                               << (minimizeInLammps ? "  (minimized)" : "  (unrelaxed)");

                    std::ostringstream manifestLine;
                    manifestLine << indexName(index) << "  " << engaged.size()
                                 << "  " << std::fixed << std::setprecision(4) << outOfPlane;
                    if (energiesRequested)
                        manifestLine << "  " << std::fixed << std::setprecision(8)
                                     << density << "  " << gbEnergy;
                    for (const int i : engaged)
                        manifestLine << "   " << describe(i);

#pragma omp critical (report)
                    {
                        std::cout << report.str() << std::endl;
                        manifest << manifestLine.str() << "\n";
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
            // the parallel region ends; only the manifest is left to close here.
            manifest.close();

            std::cout << "\nstates written : " << built << std::endl;
            std::cout << "  rejected (clash / not triangulable) : " << rejected
                      << std::endl;
            std::cout << "  candidate states examined : " << subsets.size() << std::endl;
            for (const auto& [message,count] : reasons)
                std::cout << "      " << count << " x  " << message << std::endl;
            if (energiesRequested)
                std::cout << "  energies in " << outputDirectory << "/output_thread_<id>.txt"
                          << std::endl;
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

            // The (t,s) pairs this signature engages, printed next to the signature itself.
            for (int i=0; i<stateSize; ++i) {
                if (state(i) != 1) continue;
                const Eigen::Vector3d t= ensemble.tShiftPairs[i].first.cartesian();
                const Eigen::Vector3d sh= ensemble.tShiftPairs[i].second;
                std::cout << "           t=(" << std::fixed << std::setprecision(4)
                          << t(0) << "," << t(1) << "," << t(2) << ")  s=("
                          << sh(0) << "," << sh(1) << "," << sh(2) << ")" << std::endl;
                // this path reads ready-made (t,s) signatures, so Sites does not reach it
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