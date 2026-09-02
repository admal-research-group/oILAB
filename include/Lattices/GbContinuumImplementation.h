//
// Created by Nikhil Chandra Admal on 5/27/24.
//

#ifndef OILAB_GBCONTINUUMIMPLEMENTATION_H
#define OILAB_GBCONTINUUMIMPLEMENTATION_H
#include "GbMesoState.h"

namespace oILAB {

namespace detail {
    template<int dim>
    inline Eigen::Vector3d toVec3(const typename LatticeCore<dim>::VectorDimD& v) {
        Eigen::Vector3d v3 = Eigen::Vector3d::Zero();
        for (int k = 0; k < dim && k < 3; ++k) v3(k) = v(k);
        return v3;
    }

    template<int dim, typename Pairs>
    inline std::vector<std::vector<double>> xuPairsToCloud(const Pairs& pairs) {
        std::vector<std::vector<double>> pc;
        pc.reserve(pairs.size());
        for (const auto& [x, u] : pairs) {
            const Eigen::Vector3d x3 = toVec3<dim>(x);
            const Eigen::Vector3d u3 = toVec3<dim>(u);
            pc.push_back({x3(0), x3(1), x3(2), u3(0), u3(1), u3(2)});
        }
        return pc;
    }

    /*! The deformed nodes \f$\mathbf x_i+\mathbf u_i\f$, i.e. where the surface actually sits.
     *  These are what the shared triangulation is built from. */
    template<int dim, typename Pairs>
    inline std::vector<Eigen::Vector3d> deformedNodes(const Pairs& pairs) {
        std::vector<Eigen::Vector3d> s;
        s.reserve(pairs.size());
        for (const auto& [x, u] : pairs)
            s.push_back(toVec3<dim>(x) + toVec3<dim>(u));
        return s;
    }

    /*! The two in-plane period vectors held in \p boxDim. */
    inline Eigen::Matrix<double,3,2> inPlanePeriods(const std::vector<double>& boxDim) {
        if (boxDim.size() < 6)
            throw std::runtime_error("GbContinuum: boxDim must hold the two in-plane period "
                                     "vectors as {p1x,p1y,p1z,p2x,p2y,p2z}.");
        Eigen::Matrix<double,3,2> P;
        P.col(0) << boxDim[0], boxDim[1], boxDim[2];
        P.col(1) << boxDim[3], boxDim[4], boxDim[5];
        return P;
    }

    /*! Resolves \p d into an integer multiple of the in-plane periods, or returns false if it is
     *  not one (to within \p tol). */
    inline bool asPeriodTranslation(const Eigen::Vector3d& d,
                                    const Eigen::Matrix<double,3,2>& P,
                                    const double& tol,
                                    Eigen::Vector3d& translation) {
        const Eigen::Vector2d c = P.colPivHouseholderQr().solve(d);
        const Eigen::Vector2d rounded(std::round(c(0)), std::round(c(1)));
        translation = P * rounded;
        return (d - translation).norm() <= tol;
    }

    /*! \brief Point cloud of grain B, with each node moved to the periodic image whose deformed
     *  position coincides exactly with that of the corresponding node of grain A.
     *
     *  GbMesoState::getFacetedSurfaces() forms xA = s-t/2 and xB = s+t/2 (so both deform to the
     *  CSL shift s) but then wraps xA and xB into the bicrystal box *independently*.  Nodes that
     *  wrap differently end up in different periodic images, and their deformed positions then
     *  agree only modulo the in-plane periods.  Since the box vectors are CSL vectors, shifting a
     *  node by one of them lands on another lattice site of the same grain, so undoing the
     *  mismatch here is free -- and it is what lets the two facets share a connectivity and still
     *  describe the same embedded surface rather than two copies a period apart.
     */
    template<int dim, typename Pairs>
    inline std::vector<std::vector<double>> gluedCloudB(const Pairs& pairsA,
                                                        const Pairs& pairsB,
                                                        const std::vector<double>& boxDim,
                                                        const double& tol=1.0e-8) {
        if (pairsA.size() != pairsB.size())
            throw std::runtime_error("GbContinuum: the two faceted surfaces have "
                                     + std::to_string(pairsA.size()) + " and "
                                     + std::to_string(pairsB.size()) + " nodes -- they cannot be glued.");

        const Eigen::Matrix<double,3,2> P = inPlanePeriods(boxDim);

        std::vector<std::vector<double>> pc;
        pc.reserve(pairsB.size());
        for (std::size_t i = 0; i < pairsB.size(); ++i) {
            const Eigen::Vector3d dA = toVec3<dim>(pairsA[i].first) + toVec3<dim>(pairsA[i].second);
            const Eigen::Vector3d xB = toVec3<dim>(pairsB[i].first);
            const Eigen::Vector3d uB = toVec3<dim>(pairsB[i].second);

            Eigen::Vector3d translation;
            if (!asPeriodTranslation(xB + uB - dA, P, tol, translation))
                throw std::runtime_error("GbContinuum: node " + std::to_string(i)
                    + " of the two faceted surfaces deforms to points that differ by "
                    + std::to_string((xB + uB - dA).norm())
                    + ", which is not an in-plane period of the mesostate box.  The two surfaces "
                      "cannot be glued into one.");

            const Eigen::Vector3d xBglued = xB - translation;
            pc.push_back({xBglued(0), xBglued(1), xBglued(2), uB(0), uB(1), uB(2)});
        }
        return pc;
    }
}

    template <int dim>
    GbContinuum<dim>::GbContinuum(const Eigen::Matrix<double, dim, dim - 1> &domain,
                                  const std::pair<XuPairs,XuPairs>& xuPairsOfFacetedSurfaces,
                                  const VectorDimD& normalGrainA,
                                  const std::vector<double>& boxDim,
                                  const int& facetRefinement,
                                  const bool &verbosity) :
    /*init*/ gbDomain(domain),
    /*init*/ xuPairsOfFacetedSurfaces(xuPairsOfFacetedSurfaces),
    /*init*/ normalGrainA(normalGrainA),
    /*init*/ boxDim(boxDim),
    /*init*/ sharedTopology(GbFacet::triangulate(detail::deformedNodes<dim>(xuPairsOfFacetedSurfaces.first),
                                                 boxDim,
                                                 detail::toVec3<dim>(normalGrainA))),
    // The two facets share a connectivity but face opposite ways: grain A occupies the +nA side
    // and grain B the -nA side, and each grain has to sit on the POSITIVE side of its own facet
    // for the solid angle -- hence the sign of its displacement -- to come out right.
    /*init*/ facetA(detail::xuPairsToCloud<dim>(xuPairsOfFacetedSurfaces.first),  boxDim, sharedTopology, senseA, 4, facetRefinement),
    /*init*/ facetB(detail::gluedCloudB<dim>(xuPairsOfFacetedSurfaces.first,
                                            xuPairsOfFacetedSurfaces.second,
                                            boxDim), boxDim, sharedTopology, senseB, 4, facetRefinement){

        assertFacetsGlue();

        if (verbosity) {
            std::cout << "-------------------------------------------------------------"
                         "-----------------"
                      << std::endl;
            std::cout << std::endl;
        }
    }

    template<int dim>
    GbContinuum<dim>::VectorDimD GbContinuum<dim>::displacement(const VectorDimD& x,
                                                                const int& i) const
    {
        VectorDimD u;
        u.setZero();

        /*
        We also need to be careful about the following: sometimes coincident points other
        than those engaged are mistakenly engaged, i.e. form coincidence, after deformation.
        In such cases we delete those unwanted coincidences. If done successfully, the number of
        coincidence points/density should exactly match with that expected from the mesostate signature
        */

        Eigen::Vector3d u3;
        Eigen::Vector3d x3 = to3D(x);

        if (i==1 || i==-1) // inside grain A
            u3 = facetA.displacement(x3);
        else if (i==2 || i==-2) // inside grain B
            u3 = facetB.displacement(x3);
        else
            return u;  // outside both grains — zero

        if constexpr (dim == 3)
            u = u3;
        else
            u = u3.head<dim>();

        return u;
    }

    template<int dim>
    bool GbContinuum<dim>::inGrainA(const VectorDimD &x) const {
        // check if x belongs to the grain 1, i.e. the appropriate side of the first
        // faceted surface. You will have to use the variable normalGrainA. It is sometimes
        // possible that x belongs to both grainA region and grainB region
        // normalGrainA points out of grain A, so grain A is the material *below* facet A, and an
        // atom of lattice A lying above it does not belong to the grain and has to go.  facetA was
        // built with the topology's own orientation (+1), i.e. its normals run along normalGrainA,
        // so its positive side is the one above.
        //
        // The facet's own side test is used rather than a signed distance along a fixed normal:
        // for a point down inside a corrugation of the facet the nearest point can lie on a
        // neighbouring ridge, and a distance test then reports the wrong side.
        // Atoms sitting on the facet are kept: they are the sites the two grains are brought
        // together at, so they belong to the grain rather than to neither side of it.
        return facetA.isOnSurface(to3D(x)) || sideAlongNormalGrainA(facetA, senseA, x) < 0;
    }
    template<int dim>
    bool GbContinuum<dim>::inGrainB(const VectorDimD &x) const {
        // check if x belongs to the grain 2, i.e. the appropriate side of the second
        // faceted surface. You will have to use the variable normalGrainA. It is sometimes
        // possible that x belongs to both grainA region and grainB region.
        // Grain B is the material above facet B, so an atom of lattice B below it has to go.
        // facetB carries the reversed orientation (-1), so its own positive side is the one below
        // and the sense has to be undone before comparing against normalGrainA.
        return facetB.isOnSurface(to3D(x)) || sideAlongNormalGrainA(facetB, senseB, x) > 0;
    }

    template<int dim>
    int GbContinuum<dim>::sideAlongNormalGrainA(const GbFacet& facet, const int& facetSense,
                                                const VectorDimD& x) const {
        // facet.sideOf() answers with respect to the facet's own normals, which carry the sense it
        // was built with; multiplying by that sense re-expresses the answer along normalGrainA.
        return facetSense * facet.sideOf(to3D(x));
    }

    template <int dim>
    void GbContinuum<dim>::assertFacetsGlue(const double& tol) const {
        const auto& PA = xuPairsOfFacetedSurfaces.first;
        const auto& PB = xuPairsOfFacetedSurfaces.second;

        // Nodewise: grain A node i and grain B node i must deform onto the same point.  This is
        // exact by construction in GbMesoState::getFacetedSurfaces (xA = s-t/2 with u = +t/2,
        // xB = s+t/2 with u = -t/2, so both deform to the CSL shift s), so any failure here means
        // the two faceted surfaces were not built as a matching pair.
        if (PA.size() != PB.size())
            throw std::runtime_error("GbContinuum: the two faceted surfaces have "
                                     + std::to_string(PA.size()) + " and " + std::to_string(PB.size())
                                     + " nodes -- they cannot be glued.");

        // Corresponding nodes must deform onto the same point *of the periodic surface*, i.e. up
        // to an in-plane period -- getFacetedSurfaces() wraps the two grains independently, and
        // gluedCloudB() is what removes the resulting period offsets.
        const Eigen::Matrix<double,3,2> P = detail::inPlanePeriods(boxDim);
        double worstNode = 0.0;
        std::size_t worstNodeIndex = 0;
        for (std::size_t i = 0; i < PA.size(); ++i) {
            const Eigen::Vector3d d = detail::toVec3<dim>(PB[i].first + PB[i].second)
                                    - detail::toVec3<dim>(PA[i].first + PA[i].second);
            Eigen::Vector3d translation;
            const bool isPeriod = detail::asPeriodTranslation(d, P, tol, translation);
            const double residual = isPeriod ? (d - translation).norm() : d.norm();
            if (residual > worstNode) { worstNode = residual; worstNodeIndex = i; }
        }
        if (worstNode > tol)
            throw std::runtime_error("GbContinuum: node " + std::to_string(worstNodeIndex)
                                     + " of the two faceted surfaces deforms to two different points "
                                       "(gap " + std::to_string(worstNode)
                                     + " after removing in-plane periods).");

        // Meshwise: identical connectivity, and identical deformed vertex positions.  Together
        // these say the two deformed facets are the same surface, not merely surfaces through the
        // same nodes.
        if (facetA.faces().rows() != facetB.faces().rows() || facetA.faces() != facetB.faces())
            throw std::runtime_error("GbContinuum: the two facets do not share connectivity.");

        const double worstVertex = (facetA.deformedVertices() - facetB.deformedVertices())
                                       .rowwise().norm().maxCoeff();
        if (worstVertex > tol)
            throw std::runtime_error("GbContinuum: the deformed facets differ by up to "
                                     + std::to_string(worstVertex) + " -- they do not form one surface.");
    }

    template <int dim>
    Eigen::Vector3d GbContinuum<dim>::to3D(const VectorDimD& x) const {
        Eigen::Vector3d x3 = Eigen::Vector3d::Zero();
        for (int k = 0; k < dim; ++k) x3(k) = x(k);
        return x3;
    }

    template<int dim>
    void GbContinuum<dim>::exportMesh(const std::string& name) const
    {
        facetA.export_to_vtp(name + "_facetA.vtp");
        facetB.export_to_vtp(name + "_facetB.vtp");
    }
} // namespace oILAB

#endif