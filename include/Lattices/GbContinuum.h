//
// Created by Nikhil Chandra Admal on 5/16/24.
//

#ifndef OILAB_GBPLASTICITY_H
#define OILAB_GBPLASTICITY_H

#include "LatticeCore.h"
#include "LatticeVector.h"
#include "GbFacet.h"
#include "Eigen/Dense"
#include <deque>
#include <string>
#include <utility>
#include <vector>

namespace oILAB {

    template<int dim>
    class GbContinuum {
        using VectorDimD= LatticeCore<dim>::VectorDimD;
        using XuPairs = std::deque<std::pair<VectorDimD,VectorDimD>>;
    public:
        /*! Describes the normally flat GB. */
        const Eigen::Matrix<double,dim,dim-1> gbDomain;

        /*! This is a pair of XuPairs, where each XuPairs is a deque of nodal
         * position-displacement pairs that describe the faceted boundaries and
         * their displacements of respective grains.
         */
        const std::pair<XuPairs,XuPairs> xuPairsOfFacetedSurfaces;

        /*! Describes the outward unit normal to grainA in a flat GB. This is
         * used to infer which side of the faceted boundaries grainA occupies.
         */
        const VectorDimD normalGrainA;

        /*! The two in-plane period vectors of the mesostate box, in Cartesian coordinates, laid
         * out as \f$\{p_{1x},p_{1y},p_{1z},p_{2x},p_{2y},p_{2z}\}\f$.  Declared before
         * \p sharedTopology so that it is initialised first.
         */
        const std::vector<double> boxDim;

        /*! Periodic connectivity shared by both facets.
         *
         * Built by triangulating the DEFORMED nodes.  Node \p i of grain A and node \p i of
         * grain B have the same deformed position, so giving both facets the same connectivity
         * makes their deformed surfaces identical triangle-for-triangle -- that is what glues
         * the boundary into one continuous surface.  Declared before \p facetA / \p facetB so
         * that it is initialised first.
         */
        const GbFacet::Topology sharedTopology;

        /*! Orientation each facet is built with, relative to the topology's own.  Every grain has
         * to sit on the POSITIVE side of its own facet: the solid angle is then +2*pi there and
         * the grain picks up +u rather than -u.  Grain A lies below facet A and grain B above
         * facet B, so the two senses are opposite.
         */
        static constexpr int senseA = -1;
        static constexpr int senseB = +1;

        /*! Triangulated (and periodic) representation of the faceted surface of
         * grain A, together with its nodal displacements.
         */
        GbFacet facetA;

        /*! Triangulated (and periodic) representation of the faceted surface of
         * grain B, together with its nodal displacements.
         */
        GbFacet facetB;

        /*! @param boxDim the two in-plane period vectors of the mesostate box in
         * Cartesian coordinates, laid out as
         * \f$\{p_{1x},p_{1y},p_{1z},p_{2x},p_{2y},p_{2z}\}\f$. Used by GbFacet to
         * build the periodic triangulation.
         * @param facetRefinement how finely the facets are subdivided for the displacement
         * quadrature; see GbFacet's constructor.  The triangulation only joins the nodes, so
         * without subdivision the displacement is held constant over triangles as large as the
         * node spacing.  1 restores that behaviour.
         */
        GbContinuum(const Eigen::Matrix<double, dim,dim-1>& domain,
                    const std::pair<XuPairs,XuPairs>& xuPairsOfFacetedSurfaces,
                    const VectorDimD& normalGrainA,
                    const std::vector<double>& boxDim,
                    const int& facetRefinement=4,
                    const bool& verbosity=false);

        /*! Returns the displacement at position \f$\mathbf{x}\f$ in the region occupied
         * by grain \f$i\f$. */
        VectorDimD displacement(const VectorDimD& x, const int& i) const;

        /*! \brief Tells if the point \f$\mathbf{x}\f$ belongs to undeformed Grain A, i.e. lies
         *  below facet A.  \p normalGrainA points out of grain A, so "below" means on the side
         *  \f$-\hat n_{\mathcal A}\f$. */
        bool inGrainA(const VectorDimD& x) const;

        /*! \brief Tells if the point \f$\mathbf{x}\f$ belongs to undeformed Grain B, i.e. lies
         *  above facet B, on the side \f$+\hat n_{\mathcal A}\f$. */
        bool inGrainB(const VectorDimD& x) const;

        /*! Writes the two faceted surfaces to \p name_facetA.vtp and \p name_facetB.vtp. */
        void exportMesh(const std::string& name) const;

        /*! \brief Throws unless the two facets glue into a single continuous surface.
         *
         * Checks that corresponding nodes of the two grains deform onto the same point, that the
         * two facets share connectivity, and that their deformed vertices coincide.
         */
        void assertFacetsGlue(const double& tol=1.0e-8) const;

    private:
        /*! \brief Which side of \p facet the point \p x lies on, measured along \p normalGrainA:
         *  +1 above, -1 below.  \p facetSense is the orientation the facet was built with, and is
         *  what converts the facet's own answer into one along \p normalGrainA. */
        int sideAlongNormalGrainA(const GbFacet& facet, const int& facetSense,
                                  const VectorDimD& x) const;

        /*! Zero-pads (or truncates) a \p dim-dimensional vector to 3D, as required by GbFacet. */
        Eigen::Vector3d to3D(const VectorDimD& x) const;
    };

} // namespace oILAB

#include "GbContinuumImplementation.h"
#endif //OILAB_GBPLASTICITY_H
