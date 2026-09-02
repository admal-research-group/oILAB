//
// Lists the coincidence sites available around a grain boundary, and for each one the
// atoms of A and of B that could be brought together there.
//
// This is the site-first enumeration of mesostateEnumerationDesign.tex, taken only as far
// as the listing: no mesostate is constructed, no pair record is committed to.  A site is
// a point of 1/2*DSCL -- the lattice of midpoints of A and B atoms -- inside a slab about
// the flat boundary; a candidate atom is one that can reach the site by moving no further
// than dMax.  The displacements of the two grains are independent, so the site need not be
// the midpoint of the pair that meets there.
//
// Outputs
//   * a summary of the geometry and the level structure of the sites;
//   * sites.txt, one block per site listing its candidate atoms;
//   * sites.xyz, extended XYZ for OVITO: type 1 = site, 2 = candidate A atom,
//     3 = candidate B atom, so the listing can be checked against a real structure.
//

#include "../../include/Lattices/GbMesoStateEnsemble.h"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numbers>
#include <set>
#include <vector>

using namespace oILAB;
using VectorDimD = LatticeCore<3>::VectorDimD;

namespace {

// ----------------------------------------------------------------- USER CHOICE
double slabHalfThickness = 1.0;   // Angstrom, either side of the flat GB: the region sites are drawn from
double dMax              = 1.5;   // Angstrom, how far one atom may move to reach its site
bool   writeVisualization = true;
// -----------------------------------------------------------------------------

/*! Atoms of \p lattice within \p radius of \p target, as (position, integer coordinates). */
std::vector<std::pair<VectorDimD,Eigen::Vector3i>>
atomsNear(const Lattice<3>& lattice, const VectorDimD& target, const double& radius)
{
    std::vector<std::pair<VectorDimD,Eigen::Vector3i>> out;
    const Eigen::Matrix3d M= lattice.latticeBasis, Mi= M.inverse();
    const Eigen::Vector3d centre= Mi*target;
    int range[3];
    for (int i=0; i<3; ++i) range[i]= (int) std::ceil(Mi.row(i).norm()*radius) + 1;
    for (int i=-range[0]; i<=range[0]; ++i)
    for (int j=-range[1]; j<=range[1]; ++j)
    for (int k=-range[2]; k<=range[2]; ++k)
    {
        const Eigen::Vector3i n((int)std::round(centre(0))+i,
                                (int)std::round(centre(1))+j,
                                (int)std::round(centre(2))+k);
        const VectorDimD x= M*n.cast<double>();
        if ((target-x).norm() <= radius + FLT_EPSILON) out.emplace_back(x,n);
    }
    std::sort(out.begin(), out.end(),
              [&target](const auto& l, const auto& r)
              { return (target-l.first).norm() < (target-r.first).norm(); });
    return out;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc > 1) slabHalfThickness = std::stod(argv[1]);
    if (argc > 2) dMax              = std::stod(argv[2]);

    // Sigma 5 [001](-2 1 0), the boundary of the meshing example.
    const VectorDimD axis(0,0,1);
    const double theta= 53.130102354155994249*std::numbers::pi/180;
    VectorDimD gbNormal(-2,1,0);
    const int heightScaling= 2, periodScaling= 1, axisScaling= 1;
    const double a0= 3.615;

    Eigen::Matrix3d A;
    A << 0.0, 0.5, 0.5,  0.5, 0.0, 0.5,  0.5, 0.5, 0.0;
    A= a0*A;
    Lattice<3> lattice(A);

    try
    {
        ReciprocalLatticeVector<3> rAxisGlobal(lattice.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
        Eigen::AngleAxis<double> halfRotation(theta/2, rAxisGlobal.cartesian().normalized());
        Lattice<3> latticeA(lattice.latticeBasis, halfRotation.matrix());
        Lattice<3> latticeB(lattice.latticeBasis, halfRotation.matrix().transpose());
        BiCrystal<3> bc(latticeA, latticeB, false);
        gbNormal= halfRotation.matrix()*gbNormal;
        Gb<3> gb(bc, latticeA.reciprocalLatticeDirection(gbNormal));

        ReciprocalLatticeVector<3> rAxisA(latticeA.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
        LatticeVector<3> axisA(bc.A.latticeDirection(axis).latticeVector());
        LatticeVector<3> axisC(bc.getLatticeDirectionInC(axisA).latticeVector());
        std::vector<LatticeVector<3>> cslVectors;
        cslVectors.push_back(heightScaling*bc.csl.latticeDirection(gb.nA.cartesian()).latticeVector());
        cslVectors.push_back(periodScaling*gb.getPeriodVector(rAxisA));
        cslVectors.push_back(axisScaling*axisC);

        const auto nC= bc.getReciprocalLatticeDirectionInC(gb.nB.reciprocalLatticeVector());
        const VectorDimD nHat= nC.cartesian().normalized();
        const VectorDimD p1= cslVectors[1].cartesian(), p2= cslVectors[2].cartesian();
        const double b= bc.A.latticeBasis.col(0).norm();

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Sigma = " << bc.sigma << ",  b = " << b << " A\n";
        std::cout << "GB normal (unit) = " << nHat.transpose() << "\n";
        std::cout << "in-plane periods |p1| = " << p1.norm() << " A, |p2| = " << p2.norm()
                  << " A,  cell area = " << p1.cross(p2).norm() << " A^2\n";
        std::cout << "CSL plane spacing along n = " << nC.planeSpacing() << " A\n";
        std::cout << "slab half-thickness = " << slabHalfThickness
                  << " A,  dMax = " << dMax << " A\n\n";

        // Sites: points of 1/2*DSCL in the slab, over one in-plane period cell.
        const Eigen::Matrix3d H= bc.dscl.latticeBasis/2.0, Hi= H.inverse();
        Eigen::Matrix<double,3,2> P; P.col(0)= p1; P.col(1)= p2;
        const auto Pqr= P.colPivHouseholderQr();

        Eigen::Vector3d lo= Eigen::Vector3d::Constant(1e300), hi= -lo;
        for (int a=0;a<2;++a) for (int c=0;c<2;++c) for (int g=0;g<2;++g) {
            const Eigen::Vector3d n= Hi*(a*p1 + c*p2 + (g? slabHalfThickness : -slabHalfThickness)*nHat);
            lo= lo.cwiseMin(n); hi= hi.cwiseMax(n);
        }

        std::vector<VectorDimD> sites;
        for (int i=(int)std::floor(lo(0))-1; i<=(int)std::ceil(hi(0))+1; ++i)
        for (int j=(int)std::floor(lo(1))-1; j<=(int)std::ceil(hi(1))+1; ++j)
        for (int k=(int)std::floor(lo(2))-1; k<=(int)std::ceil(hi(2))+1; ++k)
        {
            const VectorDimD s= H*Eigen::Vector3d(i,j,k);
            if (std::abs(s.dot(nHat)) > slabHalfThickness + 1.0e-9) continue;
            const Eigen::Vector2d c= Pqr.solve(s);
            if (c(0) < -1.0e-9 || c(0) >= 1.0-1.0e-9) continue;
            if (c(1) < -1.0e-9 || c(1) >= 1.0-1.0e-9) continue;
            sites.push_back(s);
        }
        std::sort(sites.begin(), sites.end(),
                  [&nHat](const VectorDimD& l, const VectorDimD& r)
                  { if (std::abs(l.dot(nHat)-r.dot(nHat)) > 1.0e-6) return l.dot(nHat) < r.dot(nHat);
                    return l(0) < r(0); });

        // The level structure: sites lie on planes parallel to the boundary.
        std::map<long long,int> levelCount;
        for (const auto& s : sites) levelCount[std::llround(s.dot(nHat)*1.0e4)]++;
        std::cout << "sites found : " << sites.size() << ",  on " << levelCount.size()
                  << " level(s) parallel to the GB\n";
        double previous= 0.0; bool first= true;
        for (const auto& [key,count] : levelCount) {
            const double h= key/1.0e4;
            std::cout << "    s.n = " << std::setw(8) << h << " A : " << count << " sites";
            if (!first) std::cout << "   (+" << h-previous << " A from the previous level)";
            std::cout << (std::abs(h) < 1.0e-6 ? "   <-- the flat GB plane\n" : "\n");
            previous= h; first= false;
        }
        std::cout << std::endl;

        // Candidate atoms per site.
        std::ofstream list("sites.txt");
        list << std::fixed << std::setprecision(6);
        list << "# Coincidence sites and the atoms that can reach them.\n";
        list << "# slabHalfThickness = " << slabHalfThickness << " A, dMax = " << dMax << " A\n";
        list << "# site  <index>  s = (x y z)  s.n = <h>  nA = <count>  nB = <count>\n";
        list << "#   A  <i j k>  x = (x y z)  t_A = (x y z)  |t_A|\n";
        list << "#   B  <i j k>  x = (x y z)  t_B = (x y z)  |t_B|\n";

        std::set<std::array<long long,3>> seenA, seenB;
        std::vector<std::pair<int,VectorDimD>> cloud;   // (type, position) for the XYZ
        std::size_t usable= 0, nodeCount= 0;
        std::size_t minA= 1000, maxA= 0, minB= 1000, maxB= 0;

        for (std::size_t idx=0; idx<sites.size(); ++idx)
        {
            const VectorDimD& s= sites[idx];
            const auto nearA= atomsNear(bc.A, s, dMax);
            const auto nearB= atomsNear(bc.B, s, dMax);
            if (nearA.empty() || nearB.empty()) continue;
            ++usable;
            nodeCount += nearA.size()*nearB.size();
            minA= std::min(minA, nearA.size()); maxA= std::max(maxA, nearA.size());
            minB= std::min(minB, nearB.size()); maxB= std::max(maxB, nearB.size());

            list << "site " << idx << "   s = " << s.transpose()
                 << "   s.n = " << s.dot(nHat)
                 << "   nA = " << nearA.size() << "   nB = " << nearB.size() << "\n";
            for (const auto& [x,n] : nearA) {
                const VectorDimD t= s-x;
                list << "   A  " << n.transpose() << "   x = " << x.transpose()
                     << "   t_A = " << t.transpose() << "   " << t.norm() << "\n";
            }
            for (const auto& [x,n] : nearB) {
                const VectorDimD t= s-x;
                list << "   B  " << n.transpose() << "   x = " << x.transpose()
                     << "   t_B = " << t.transpose() << "   " << t.norm() << "\n";
            }
            cloud.emplace_back(1,s);
            for (const auto& [x,n] : nearA)
                if (seenA.insert({n(0),n(1),n(2)}).second) cloud.emplace_back(2,x);
            for (const auto& [x,n] : nearB)
                if (seenB.insert({n(0),n(1),n(2)}).second) cloud.emplace_back(3,x);
        }
        list.close();

        std::cout << "sites with atoms of both grains in range : " << usable
                  << " of " << sites.size() << "\n";
        if (usable) {
            std::cout << "   candidate A atoms per site : " << minA << " to " << maxA << "\n";
            std::cout << "   candidate B atoms per site : " << minB << " to " << maxB << "\n";
            std::cout << "   distinct A atoms used : " << seenA.size()
                      << ",  distinct B atoms used : " << seenB.size() << "\n";
            std::cout << "   (site, A atom, B atom) combinations : " << nodeCount << "\n";
        }
        std::cout << "\nwrote sites.txt" << std::endl;

        if (writeVisualization && !cloud.empty())
        {
            // One in-plane cell, the slab thickness, and enough out-of-plane room to see the
            // atoms that reach into it.
            std::ofstream xyz("sites.xyz");
            xyz << cloud.size() << "\n";
            xyz << "Lattice=\" " << std::setprecision(12)
                << (4*(slabHalfThickness+dMax)*nHat).transpose() << " "
                << p1.transpose() << " " << p2.transpose()
                << "\" Properties=species:S:1:pos:R:3 PBC=\"F T T\" origin=\" "
                << (-2*(slabHalfThickness+dMax)*nHat).transpose() << "\"\n";
            xyz << std::setprecision(6);
            for (const auto& [type,x] : cloud)
                xyz << (type==1 ? "S" : (type==2 ? "A" : "B")) << " " << x.transpose() << "\n";
            std::cout << "wrote sites.xyz  (species S = coincidence site, A / B = candidate atoms)"
                      << std::endl;
        }
    }
    catch(std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    return 0;
}
