//
// Enumerates the (t,s) pairs a grain boundary admits, and checks the invariants the rest of the
// code relies on.  Choose between the two searches with the constants marked USER CHOICE below,
// or override them on the command line:
//
//     testGbShiftsNonFlat [flat|full] [tMax] [sPerpMax] [tPerpMax] [dedup 0|1]
//
// Flat is the original enumeration, restricted to translations whose CSL shift lies along the
// boundary.  Full also admits shifts that leave the boundary plane, which is what allows non-flat
// mesostates; it is drawn from a ball |t| <= tMax*b intersected with a slab |t.n| <= tPerpMax*b.
//

#include "../../include/Lattices/GbMesoStateEnsemble.h"

#include <cassert>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using namespace oILAB;

namespace {

// ----------------------------------------------------------------- USER CHOICE
GbShiftSearch searchMode              = GbShiftSearch::Full;  // Flat or Full
double        tMax                = 2.0;   // ball radius, in units of b
double        sPerpMax            = 1.0;   // shift allowed off the boundary, in units of b

// The shift limit is more natural to think about as a thickness: the shifts are confined to a
// layer reaching this far either side of the flat boundary.  Set it positive to have it override
// sPerpMax once b is known (the filter is |s.n| <= sPerpMax*b/2, so sPerpMax = 2*halfThickness/b);
// set it <= 0 to use sPerpMax as given.  Passing sPerpMax on the command line also disables it.
double        layerHalfThickness  = 0.41;  // Angstrom, either side of the flat GB
double        tPerpMax            = 0.5;   // slab half-thickness, in units of b (Full only)
bool          oneTranslationPerSite = false; // keep only the shortest t per site (Full only)

// Dump the atoms behind the pairs as extended XYZ, for viewing in OVITO.
bool          writeVisualization  = true;
std::string   visualizationBase   = "tsPairAtoms";
// -----------------------------------------------------------------------------

int failures = 0;

void check(const bool ok, const std::string& what, const double& value = 0.0)
{
    std::cout << (ok ? "  [ ok ] " : "  [FAIL] ") << what;
    if (!ok || value != 0.0) std::cout << "   (" << value << ")";
    std::cout << std::endl;
    if (!ok) ++failures;
}

bool onLattice(const Eigen::Vector3d& x, const Lattice<3>& lat, const double& tol = 1.0e-6)
{
    const Eigen::Vector3d nd = lat.reciprocalBasis.transpose() * x;
    return (nd - nd.array().round().matrix()).norm() <= tol;
}

/*! One atom of the dump. */
struct Atom
{
    Eigen::Vector3d position;
    Eigen::Vector3d displacement;   // where the atom has to move to close the boundary
    Eigen::Vector3d translation;    // the t of the pair it belongs to
    int             type;           // see the table in collectAtoms()
    int             pairIndex;      // -1 for the bulk context atoms
    int             newCsl;         // 1 = its coincidence is new, 0 = already a CSL point, -1 = bulk
    double          radius;
};

/*! \brief The atoms that give rise to the (t,s) pairs.
 *
 *  Each pair contributes three: the site \f$x_{\mathcal A}=s-t/2\f$ on grain A, the site
 *  \f$x_{\mathcal B}=s+t/2\f$ on grain B, and the point \f$s\f$ where the two arrive once A is
 *  displaced by \f$+t/2\f$ and B by \f$-t/2\f$.  Those displacements are carried on the atoms, so
 *  warping by them in a viewer collapses each A/B pair onto its coincidence site -- which is the
 *  interface closing up.
 *
 *  A shift does not always produce a coincidence that was not there before: if \f$s\f$ is itself a
 *  point of the CSL then the two grains already coincided at that spot, whereas otherwise the
 *  enumerated shift has created a *new* coincidence.  The two cases get separate types so that
 *  colouring by type in a viewer distinguishes them, and the \p newCSL field carries the same
 *  information for filtering.
 *
 *  Types extend the convention BiCrystal::box() already writes (1 = lattice A, 2 = lattice B,
 *  3 = CSL):
 *      1 / 2 / 3   grain A atom, grain B atom, and their meeting point, at an existing CSL point
 *      4 / 5 / 6   the same three, where the shift forms a NEW CSL point
 */
std::vector<Atom> collectAtoms(const Gb<3>& gb,
                               const std::vector<std::pair<LatticeVector<3>,Eigen::Vector3d>>& pairs)
{
    std::vector<Atom> atoms;

    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const Eigen::Vector3d t = pairs[i].first.cartesian();
        const Eigen::Vector3d s = pairs[i].second;
        const Eigen::Vector3d u = t/2;
        const int index = static_cast<int>(i);

        // Was this already a coincidence, or does the shift create one?
        const bool fresh = !onLattice(s, gb.bc.csl);
        const int base = fresh ? 4 : 1;
        const int flag = fresh ? 1 : 0;

        atoms.push_back({s - u,  u, t, base + 0, index, flag, 0.40});   // grain A atom
        atoms.push_back({s + u, -u, t, base + 1, index, flag, 0.40});   // grain B atom
        atoms.push_back({s, Eigen::Vector3d::Zero(), t, base + 2, index, flag, 0.25}); // meeting point
    }

    return atoms;
}

/*! Writes the atoms as extended XYZ, in the same style BiCrystal::box() uses, for OVITO. */
void writeXyz(const std::string& filename, const std::vector<Atom>& atoms,
              const std::vector<LatticeVector<3>>& cslVectors)
{
    std::ofstream os(filename);
    if (!os) { std::cout << "  [warn] could not open " << filename << std::endl; return; }
    os << std::setprecision(17);
    os << atoms.size() << "\n";
    os << "Lattice=\"";
    for (int i = 0; i < 3; ++i) os << cslVectors[i].cartesian().transpose() << " ";
    os << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1:displacement:R:3:"
          "newCSL:I:1:pairIndex:I:1 PBC=\"F T T\"\n";
    for (const auto& a : atoms)
        os << a.type << " " << a.position.transpose() << " " << a.radius << " "
           << a.displacement.transpose() << " " << a.newCsl << " " << a.pairIndex << "\n";
    std::cout << "  [info] wrote " << filename << " (" << atoms.size() << " atoms)" << std::endl;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc > 1) {
        const std::string mode(argv[1]);
        if      (mode == "flat") searchMode = GbShiftSearch::Flat;
        else if (mode == "full") searchMode = GbShiftSearch::Full;
        else { std::cout << "First argument must be \"flat\" or \"full\"." << std::endl; return -1; }
    }
    if (argc > 2) tMax                  = std::stod(argv[2]);
    if (argc > 3) { sPerpMax = std::stod(argv[3]); layerHalfThickness = -1.0; }
    if (argc > 4) tPerpMax              = std::stod(argv[4]);
    if (argc > 5) oneTranslationPerSite = (std::stoi(argv[5]) != 0);

    const bool full = (searchMode == GbShiftSearch::Full);

    try {
        // ---- sigma 5 [001](-2 1 0) in fcc, the same boundary examples/GbFacetMeshing uses
        const Eigen::Vector3d axis(0,0,1);
        const double theta = 53.130102354155994249 * std::numbers::pi/180;
        Eigen::Vector3d gbNormal(-2,1,0);
        const double a0 = 3.615;

        Eigen::Matrix3d A;
        A << 0.0, 0.5, 0.5,
             0.5, 0.0, 0.5,
             0.5, 0.5, 0.0;
        A *= a0;
        Lattice<3> lattice(A);

        ReciprocalLatticeVector<3> rAxisGlobal(
            lattice.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
        const Eigen::AngleAxis<double> halfRotation(theta/2, rAxisGlobal.cartesian().normalized());
        Lattice<3> latticeA(lattice.latticeBasis, halfRotation.matrix());
        Lattice<3> latticeB(lattice.latticeBasis, halfRotation.matrix().transpose());
        BiCrystal<3> bc(latticeA, latticeB, false);

        gbNormal = halfRotation.matrix() * gbNormal;
        Gb<3> gb(bc, latticeA.reciprocalLatticeDirection(gbNormal));

        ReciprocalLatticeVector<3> rAxisA(
            latticeA.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
        LatticeVector<3> axisA(gb.bc.A.latticeDirection(axis).latticeVector());
        LatticeVector<3> axisC(gb.bc.getLatticeDirectionInC(axisA).latticeVector());

        std::vector<LatticeVector<3>> cslVectors;
        cslVectors.push_back(2 * gb.bc.csl.latticeDirection(gb.nA.cartesian()).latticeVector());
        cslVectors.push_back(1 * gb.getPeriodVector(rAxisA));
        cslVectors.push_back(1 * axisC);

        const double b = gb.bc.A.latticeBasis.col(0).norm();
        if (layerHalfThickness > 0.0) sPerpMax = 2.0*layerHalfThickness/b;
        auto nC = gb.bc.getReciprocalLatticeDirectionInC(gb.nB.reciprocalLatticeVector());
        const Eigen::Vector3d nHat = nC.cartesian().normalized();

        std::cout << "\n=============================================================" << std::endl;
        std::cout << "search   = " << (full ? "Full (non-flat allowed)" : "Flat (original)") << std::endl;
        std::cout << "b        = " << b << std::endl;
        std::cout << "tMax     = " << tMax << " b = " << tMax*b << std::endl;
        std::cout << "sPerpMax = " << sPerpMax << " b = " << sPerpMax*b
                  << "   (shifts confined to +-" << sPerpMax*b/2 << " A of the flat GB";
        if (layerHalfThickness > 0.0) std::cout << ", set from layerHalfThickness";
        std::cout << ")" << std::endl;
        if (full) {
            std::cout << "tPerpMax = " << tPerpMax << " b = " << tPerpMax*b
                      << (tPerpMax >= tMax ? "   (slab inactive)" : "") << std::endl;
            std::cout << "one translation per site = " << (oneTranslationPerSite ? "yes":"no")
                      << std::endl;
        }
        std::cout << "=============================================================\n" << std::endl;

        GbMesoStateEnsemble<3> ensemble(gb, rAxisA, cslVectors,
                                        tMax, sPerpMax,
                                        searchMode, tPerpMax, oneTranslationPerSite,
                                        full ? "translationsNonFlat.txt" : "");

        const auto& pairs = ensemble.tShiftPairs;

        // ---- the full list
        std::cout << "\n" << pairs.size() << " (t,s) pairs:" << std::endl;
        std::cout << std::fixed << std::setprecision(6);
        std::cout << std::setw(5) << "idx"
                  << std::setw(31) << "t"
                  << std::setw(31) << "s"
                  << std::setw(11) << "|t|/b"
                  << std::setw(11) << "t.n/b"
                  << std::setw(11) << "s.n/b" << std::endl;
        for (std::size_t i = 0; i < pairs.size(); ++i) {
            const Eigen::Vector3d t = pairs[i].first.cartesian();
            const Eigen::Vector3d s = pairs[i].second;
            std::cout << std::setw(5) << i
                      << std::setw(11) << t(0) << std::setw(10) << t(1) << std::setw(10) << t(2)
                      << std::setw(11) << s(0) << std::setw(10) << s(1) << std::setw(10) << s(2)
                      << std::setw(11) << t.norm()/b
                      << std::setw(11) << t.dot(nHat)/b
                      << std::setw(11) << s.dot(nHat)/b << std::endl;
        }
        std::cout << std::endl;

        // ---- invariants
        check(!pairs.empty(), "the enumeration is non-empty",
              static_cast<double>(pairs.size()));

        // The ensemble always engages pair 0, so it must be the zero translation.
        check(pairs[0].first.cartesian().norm() < 1.0e-9,
              "pair 0 is the zero translation (the ensemble always engages it)",
              pairs[0].first.cartesian().norm());

        // Every pair must respect the limits it was generated under.
        double worstT = 0.0, worstTn = 0.0, worstSn = 0.0;
        for (const auto& [t, s] : pairs) {
            worstT  = std::max(worstT,  t.cartesian().norm());
            worstTn = std::max(worstTn, std::abs(t.cartesian().dot(nHat)));
            worstSn = std::max(worstSn, std::abs(s.dot(nHat)));
        }
        check(worstT <= tMax*b + 1.0e-6, "every |t| is within tMax*b", worstT);
        check(worstSn <= sPerpMax*b/2 + 1.0e-6, "every |s.n| is within sPerpMax*b/2", worstSn);
        if (full && tPerpMax < tMax)
            check(worstTn <= tPerpMax*b + 1.0e-6, "every |t.n| is within the slab", worstTn);

        // The precondition GbMesoState::getFacetedSurfaces enforces with exit(0): the nodes the
        // pair generates have to be lattice vectors of the respective grains.
        int offA = 0, offB = 0;
        for (const auto& [t, s] : pairs) {
            const Eigen::Vector3d u = t.cartesian()/2;
            if (!onLattice(s - u, gb.bc.A)) ++offA;
            if (!onLattice(s + u, gb.bc.B)) ++offB;
        }
        check(offA == 0, "every xA = s - t/2 lies on lattice A", offA);
        check(offB == 0, "every xB = s + t/2 lies on lattice B", offB);

        // Distinct pairs; and, when asked for, distinct sites.
        int dupPairs = 0, dupSites = 0;
        for (std::size_t i = 0; i < pairs.size(); ++i)
            for (std::size_t j = i+1; j < pairs.size(); ++j) {
                const bool sameT = (pairs[i].first.cartesian()-pairs[j].first.cartesian()).norm() < 1.0e-9;
                const bool sameS = (pairs[i].second-pairs[j].second).norm() < 1.0e-9;
                if (sameT && sameS) ++dupPairs;
                if (sameS) ++dupSites;
            }
        check(dupPairs == 0, "no duplicated (t,s) pair", dupPairs);
        if (full && oneTranslationPerSite)
            check(dupSites == 0, "one translation per site: every site is used once", dupSites);
        else
            std::cout << "  [info] pairs sharing a site: " << dupSites
                      << "  (engaging two of them together is rejected downstream as a clash)"
                      << std::endl;

        // What is actually available for building relief.
        std::vector<double> heights;
        for (const auto& [t, s] : pairs) {
            const double h = s.dot(nHat);
            if (std::none_of(heights.begin(), heights.end(),
                             [&](double y){ return std::abs(h-y) < 1.0e-6; }))
                heights.push_back(h);
        }
        std::sort(heights.begin(), heights.end());
        std::cout << "  [info] " << heights.size() << " distinct step heights s.n/b :";
        for (const double h : heights) std::cout << " " << h/b;
        std::cout << std::endl;
        check(heights.size() >= 1, "at least one step height is available",
              static_cast<double>(heights.size()));
        if (full)
            check(heights.size() > 1,
                  "the full search reaches shifts off the boundary plane",
                  static_cast<double>(heights.size()));
        // ---- the atoms behind the pairs
        if (writeVisualization) {
            const auto atoms = collectAtoms(gb, pairs);
            writeXyz(visualizationBase + ".xyz", atoms, cslVectors);

            // Every A/B site must land on its coincidence site once displaced -- the same property
            // the two facets rely on, checked here at the level of the individual atoms.
            // Atoms come in triples: grain A atom, grain B atom, meeting point -- at types
            // 1/2/3 for an existing coincidence and 4/5/6 for a new one, so the pattern has to be
            // recognised modulo that offset rather than against the literal types 1/2/3.
            double worstGap = 0.0;
            int checkedTriples = 0;
            for (std::size_t i = 0; i + 2 < atoms.size(); i += 3) {
                const int base = (atoms[i].type >= 4) ? 4 : 1;
                if (atoms[i].type != base || atoms[i+1].type != base+1 || atoms[i+2].type != base+2)
                    break;
                worstGap = std::max(worstGap,
                    ((atoms[i].position   + atoms[i].displacement)   - atoms[i+2].position).norm());
                worstGap = std::max(worstGap,
                    ((atoms[i+1].position + atoms[i+1].displacement) - atoms[i+2].position).norm());
                ++checkedTriples;
            }
            check(checkedTriples == static_cast<int>(pairs.size()),
                  "every pair contributed a triple to the dump",
                  static_cast<double>(checkedTriples));
            check(worstGap <= 1.0e-9,
                  "displaced A and B atoms meet at their coincidence site", worstGap);

            int freshPairs = 0;
            for (const auto& [t, s] : pairs) if (!onLattice(s, gb.bc.csl)) ++freshPairs;
            std::cout << "  [info] " << freshPairs << " of " << pairs.size()
                      << " pairs form a NEW CSL point; " << (pairs.size()-freshPairs)
                      << " sit on an existing one" << std::endl;
            std::cout << "  [info] atom_type: 1/2/3 = A atom, B atom, meeting point at an EXISTING "
                         "CSL point\n"
                         "                    4/5/6 = the same three where the shift forms a NEW "
                         "CSL point\n"
                         "         newCSL is 1 for a new coincidence and 0 for an existing one.\n"
                         "         Displacing by \"displacement\" brings each A/B atom onto its "
                         "coincidence point." << std::endl;
        }
    }
    catch (const std::exception& e) {
        std::cout << "  [FAIL] threw: " << e.what() << std::endl;
        ++failures;
    }

    std::cout << std::endl
              << (failures ? "testGbShiftsNonFlat FAILED with " : "testGbShiftsNonFlat passed, ")
              << failures << " failure(s)." << std::endl;
    return failures ? -1 : 0;
}
