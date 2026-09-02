//
// Created by Nikhil Chandra Admal on 2/4/24.
//
#include <cassert>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "../../include/Lattices/GbShifts.h"
#include "../../include/Utilities/randomInteger.h"

namespace oILAB {
template <int dim>
GbShifts<dim>::GbShifts(const Gb<dim> &gb,
                        const ReciprocalLatticeVector<dim>& axis,
                        const std::vector<LatticeVector<dim>>& gbCslVectors,
                        const double& tMax,
                        const double& sPerpMax,
                        const GbShiftSearch& search,
                        const double& tPerpMax,
                        const bool& oneTranslationPerSite,
                        const std::string& filename) try :
/* init */ gb(gb)
/* init */,axis(axis)
/* init */,gbCslVectors(gbCslVectors)
/* init */,tShiftPairs(search==GbShiftSearch::Flat ?
                       getbShiftPairs(gb, gbCslVectors, tMax, sPerpMax) :
                       getNonFlatShiftPairs(gb, gbCslVectors, tMax, sPerpMax,
                                            tPerpMax, oneTranslationPerSite, filename))
/* init */,search(search)
{
    std::cout << "--------------------GBShifts class construction "
                 "---------------------------"
              << std::endl;
    std::cout << "GB CSL vectors = " << std::endl;
    for (const auto& elem : gbCslVectors)
        std::cout << elem.cartesian().transpose() << std::endl;
    std::cout << std::endl;

    double latticeConstant= gb.bc.A.latticeBasis.col(0).norm();
    std::cout << "Maximum translation = "
              << tMax * latticeConstant << std::endl;
    std::cout << std::endl;
    std::cout << "Maximum shift offset = "
              << sPerpMax * latticeConstant << std::endl;
    std::cout << std::endl;
    std::cout << "Search = "
              << (search==GbShiftSearch::Flat ? "Flat (CSL shifts along the boundary)"
                                             : "Full (non-flat shifts allowed)") << std::endl;
    if (search==GbShiftSearch::Full) {
        std::cout << "Maximum translation off the boundary plane = "
                  << (tPerpMax>1.0e299 ? tMax*latticeConstant : tPerpMax*latticeConstant)
                  << (tPerpMax>=tMax ? "  (slab inactive: the ball is already this thin)" : "")
                  << std::endl;
        std::cout << "One translation per site = "
                  << (oneTranslationPerSite ? "yes" : "no") << std::endl;
    }
    std::cout << std::endl;

    VectorDimD normal;
    if (dim == 3)
        normal = gbCslVectors[0].cross(gbCslVectors[1]).cartesian().normalized();
    else
        normal = gbCslVectors[0].cross().cartesian().normalized();

    Eigen::IOFormat fmt(12, 0, " ", " ", "", "", "", "");
    std::cout << std::fixed << std::setprecision(6);


    std::cout << "Exploring the following translation-shift pairs:" << std::endl;
    for (const auto &[t, s] : tShiftPairs) {
        std::cout << "t = " << t.cartesian().transpose().format(fmt);
        std::cout << "; s = " << s.transpose().format(fmt) << std::endl;
        assert(t.cartesian().norm() <= tMax*latticeConstant && "Translation exceeds the max translation.\n");
        assert(s.dot(normal) <= sPerpMax*latticeConstant/2.0 && "Shift exceeds the max offset along the GB normal.\n");
    }
    std::cout << "----------------------------" << std::endl;
    std::cout << std::endl;
}
catch(std::runtime_error& e)
{
    std::cout << e.what() << std::endl;
    throw(std::runtime_error("GB construction failed. "));
}

template<int dim>
std::vector<std::pair<LatticeVector<dim>, typename GbShifts<dim>::VectorDimD>> GbShifts<dim>::getbShiftPairs(const Gb<dim>& gb,
                                                                                                             const std::vector<LatticeVector<dim>>& gbCslVectors,
                                                                                                             const double& tMax,
                                                                                                             const double& sPerpMax)
{
    std::vector<std::pair<LatticeVector<dim>, VectorDimD>> output;

    // ensure the input CSL vectors describe the GB.
    assert(gbCslVectors.size()==dim-1);

    // form the CSL cell for modulo operations to identify the shift vectors
    // - the cell is spanned by the input GB CSL vectors and an out-of-plane vector scaled
    //   by a factor determined by sPerpMax
    auto nC= gb.bc.getReciprocalLatticeDirectionInC(gb.nB.reciprocalLatticeVector());
    auto gbPlaneParallelCslBasis= gb.bc.csl.planeParallelLatticeBasis(nC,true);
    std::vector<LatticeVector<dim>> cslSubLatticeVectors;
    double latticeConstant= gb.bc.A.latticeBasis.col(0).norm();
    int factor= floor(sPerpMax*latticeConstant/nC.planeSpacing() +FLT_EPSILON);
    factor= (factor>0 ? factor : 1);
    cslSubLatticeVectors.push_back(factor*gbPlaneParallelCslBasis[0].latticeVector());
    cslSubLatticeVectors.push_back(gbCslVectors[0]);
    cslSubLatticeVectors.push_back(gbCslVectors[1]);
    auto cslPoints= gb.bc.csl.box(cslSubLatticeVectors,"cslSubLattice.txt");

    // Collect all DSCL translations within a parallelopiped of size tMax
    //  - orient this parallelopiped along the tilt axis, period vector directions,
    //    and the third direction is close to being parallel to the GB normal
    std::vector<LatticeVector<dim>> orthogonalCslLatticeVectors;
    for (int i=0; i<dim; ++i)
        orthogonalCslLatticeVectors.push_back(cslSubLatticeVectors[i]);
    gb.bc.updateBoxVectors(orthogonalCslLatticeVectors,0.8);
    std::vector<LatticeVector<dim>> latticeVectorsDscl;
    for (int i=0; i<dim; ++i)
        latticeVectorsDscl.push_back(gb.bc.getLatticeDirectionInD(orthogonalCslLatticeVectors[i]).latticeVector());
    for(int i=0; i<dim; ++i)
    {
        int factor= floor(tMax*latticeConstant/latticeVectorsDscl[i].cartesian().norm()+FLT_EPSILON);
        factor= (factor>0 ? factor : 1);
        latticeVectorsDscl[i]= factor*latticeVectorsDscl[i];
    }
    auto allTranslations= gb.bc.dscl.box(latticeVectorsDscl,"translations.txt");

    VectorDimD shiftT, shiftC;
    shiftT << -0.5, -0.5, -0.5;
    shiftC << -0.5, -FLT_EPSILON, -FLT_EPSILON;
    for(auto& translation: allTranslations) {
        LatticeVector<dim>::modulo(translation, latticeVectorsDscl, shiftT);
        if (translation.cartesian().norm() > tMax*gb.bc.A.latticeBasis.col(0).norm())
            continue;
        auto cslShift = LatticeVector<dim>((gb.bc.LambdaA * translation).eval(), gb.bc.dscl);
        for(const auto& cslPoint : cslPoints) {
            VectorDimD cslShiftCentered = cslShift.cartesian() + cslPoint.cartesian() - translation.cartesian() / 2;
            LatticeVector<dim>::modulo(cslShiftCentered, cslSubLatticeVectors, shiftC);
            if (abs(cslShiftCentered.dot(nC.cartesian().normalized())) <= sPerpMax*latticeConstant/2.0)
                output.push_back(std::make_pair(translation, cslShiftCentered));
        }
    }
    return output;
}
template<int dim>
std::vector<std::pair<LatticeVector<dim>, typename GbShifts<dim>::VectorDimD>>
GbShifts<dim>::getNonFlatShiftPairs(const Gb<dim>& gb,
                                   const std::vector<LatticeVector<dim>>& gbCslVectors,
                                   const double& tMax,
                                   const double& sPerpMax,
                                   const double& tPerpMax,
                                   const bool& oneTranslationPerSite,
                                   const std::string& filename)
{
    assert(gbCslVectors.size()==dim-1);
    std::vector<std::pair<LatticeVector<dim>, VectorDimD>> output;

    const double latticeConstant= gb.bc.A.latticeBasis.col(0).norm();

    // The mesostate CSL cell that the shifts are reduced into: the two in-plane GB CSL vectors,
    // plus enough out-of-plane CSL planes to accommodate sPerpMax.  Same cell the flat enumeration
    // uses, so the two searches place their shifts on the same footing.
    auto nC= gb.bc.getReciprocalLatticeDirectionInC(gb.nB.reciprocalLatticeVector());
    auto gbPlaneParallelCslBasis= gb.bc.csl.planeParallelLatticeBasis(nC,true);
    int factor= floor(sPerpMax*latticeConstant/nC.planeSpacing() +FLT_EPSILON);
    factor= (factor>0 ? factor : 1);
    std::vector<LatticeVector<dim>> cslSubLatticeVectors;
    cslSubLatticeVectors.push_back(factor*gbPlaneParallelCslBasis[0].latticeVector());
    cslSubLatticeVectors.push_back(gbCslVectors[0]);
    cslSubLatticeVectors.push_back(gbCslVectors[1]);

    // The CSL sites inside that cell.  Each translation is paired with every one of them, which is
    // what lets a single translation appear at several places along the boundary.  Without this the
    // shift would be a function of the translation alone, and a flat boundary -- many sites sharing
    // one translation -- would be unreachable: its sites would all fold onto a single lattice point
    // and GbMesoState would reject every combination of them as a clash.
    const auto cslPoints= gb.bc.csl.box(cslSubLatticeVectors,"");

    const VectorDimD nHat= nC.cartesian().normalized();
    const double tRadius= tMax*latticeConstant;
    const double tSlab= (tPerpMax>1.0e299 ? tRadius : tPerpMax*latticeConstant);
    const double sLimit= sPerpMax*latticeConstant/2.0;

    // Enumerate the DSCL translations in the ball.  With t = D n the integer coordinates obey
    // |n_i| = |row_i(D^-1).t| <= ||row_i(D^-1)|| |t|, which bounds the scan exactly.
    const Eigen::Matrix<double,dim,dim> D= gb.bc.dscl.latticeBasis;
    const Eigen::Matrix<double,dim,dim> Dinv= D.inverse();
    Eigen::Matrix<int,dim,1> limit;
    for(int i=0; i<dim; ++i)
        limit(i)= (int) std::ceil(Dinv.row(i).norm()*tRadius) + 1;

    VectorDimD shiftC;
    shiftC << -0.5, -FLT_EPSILON, -FLT_EPSILON;

    for(int i=-limit(0); i<=limit(0); ++i)
    for(int j=-limit(1); j<=limit(1); ++j)
    for(int k=-limit(2); k<=limit(2); ++k)
    {
        VectorDimI integerCoordinates;
        integerCoordinates << i,j,k;
        LatticeVector<dim> translation(integerCoordinates, gb.bc.dscl);

        const VectorDimD tCartesian= translation.cartesian();
        if (tCartesian.norm() > tRadius+FLT_EPSILON) continue;              // ball
        if (std::abs(tCartesian.dot(nHat)) > tSlab+FLT_EPSILON) continue;   // slab

        // Lambda_A t is the CSL shift produced by translating lattice A by t.
        const LatticeVector<dim> cslShift((gb.bc.LambdaA*translation).eval(), gb.bc.dscl);

        for(const auto& cslPoint : cslPoints)
        {
        VectorDimD s= cslShift.cartesian() + cslPoint.cartesian() - tCartesian/2;

        // Reduce the component off the boundary plane, so that the shift sits in the layer that
        // sPerpMax admits.
        LatticeVector<dim>::modulo(s, cslSubLatticeVectors, shiftC);

        // Then bring the in-plane part into [0,1) of the box the caller asked for.  The reduction
        // above cannot do that on its own: its out-of-plane vector comes from
        // planeParallelLatticeBasis() and is sheared -- it carries in-plane components -- so a
        // shift reduced in that basis can still sit up to a quarter of a period outside the box.
        // The two in-plane CSL vectors are lattice translations lying in the boundary plane, so
        // subtracting whole multiples of them moves the shift onto an equivalent site and leaves
        // s.n untouched, which is why this cannot disturb the filter above.
        Eigen::Matrix<double,dim,dim-1> inPlane;
        for(int c=0; c<dim-1; ++c) inPlane.col(c)= gbCslVectors[c].cartesian();
        const Eigen::Matrix<double,dim-1,1> coefficients=
            inPlane.colPivHouseholderQr().solve(s);
        // Snap before flooring.  A shift that belongs at coordinate 1 -- the periodic image of 0 --
        // routinely lands a rounding error below it, and a bare floor() leaves it there: the site
        // then sits on the far face of the box and is counted as distinct from its own image.
        Eigen::Matrix<double,dim-1,1> whole;
        for(int c=0; c<dim-1; ++c) whole(c)= std::floor(coefficients(c) + FLT_EPSILON);
        s-= inPlane*whole;

        if (std::abs(s.dot(nHat)) <= sLimit)
            output.push_back(std::make_pair(translation,s));
        }
    }

    // Shortest translation first, so that the zero pair -- which the ensemble always engages --
    // heads the list, and so that the ordering (hence the meaning of a signature) is reproducible.
    std::sort(output.begin(), output.end(),
              [](const std::pair<LatticeVector<dim>,VectorDimD>& l,
                 const std::pair<LatticeVector<dim>,VectorDimD>& r)
              {
                  const double nl= l.first.cartesian().norm();
                  const double nr= r.first.cartesian().norm();
                  if (std::abs(nl-nr) > FLT_EPSILON) return nl < nr;
                  return std::lexicographical_compare(l.first.data(), l.first.data()+dim,
                                                      r.first.data(), r.first.data()+dim);
              });

    if (oneTranslationPerSite)
    {
        // The list is sorted by |t|, so the first pair reaching a site is the shortest translation
        // that lands there.
        std::vector<std::pair<LatticeVector<dim>, VectorDimD>> unique;
        for(const auto& candidate : output)
        {
            bool seen= false;
            for(const auto& kept : unique)
                if ((candidate.second-kept.second).norm() < FLT_EPSILON) { seen= true; break; }
            if (!seen) unique.push_back(candidate);
        }
        output= unique;
    }

    if (!filename.empty())
    {
        std::ofstream file(filename);
        if (!file)
            std::cout << "Warning: could not open " << filename << " for the (t,s) pair list."
                      << std::endl;
        else
        {
            file << std::scientific << std::setprecision(15);
            file << "# non-flat (t,s) pairs\n";
            file << "# b = " << latticeConstant
                 << "  tMax = " << tMax << "b  tPerpMax = " << (tPerpMax>1.0e299 ? tMax : tPerpMax)
                 << "b  sPerpMax = " << sPerpMax << "b\n";
            file << "# index  t_x t_y t_z  s_x s_y s_z  |t|  s.n\n";
            for(std::size_t p=0; p<output.size(); ++p)
            {
                const VectorDimD& t= output[p].first.cartesian();
                const VectorDimD& s= output[p].second;
                file << p;
                for(int c=0;c<dim;++c) file << " " << t(c);
                for(int c=0;c<dim;++c) file << " " << s(c);
                file << " " << t.norm() << " " << s.dot(nHat) << "\n";
            }
        }
    }

    return output;
}

/*
template<int dim>
std::vector<LatticeVector<dim>> GbShifts<dim>::getGbCslVectors(const Gb<dim>& gb, const ReciprocalLatticeVector<dim>& axis)
{
    std::vector<LatticeVector<dim>> output;
    LatticeVector<dim> axisA(gb.bc.A.latticeDirection(axis.cartesian()).latticeVector());
    LatticeVector<dim> axisC(gb.bc.getLatticeDirectionInC(axisA).latticeVector());
    for(int i=0; i<dim-1; ++i) {
        if (i==0) output.push_back(gb.getPeriodVector(axis));
        if (i==1) output.push_back(axisC);
    }
    return output;
}
 */

//template class GbShifts<2>;
template class GbShifts<3>;

} // namespace oILAB
