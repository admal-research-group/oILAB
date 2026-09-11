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
#include <cstdlib>
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

        /*! How many shells of periodic images the facet quadrature integrates explicitly.
         *
         *  The remainder is closed off analytically, and that closure absorbs the mean nodal
         *  displacement exactly, so this controls only how well the variation about that mean is
         *  resolved far from the facet.  It converges as 1/n while the cost grows as (2n+1)^2. */
        /*! A convergence-study hook: both quadrature parameters may be overridden from the
         *  environment, so that a scan over them needs neither a recompile per point nor an edit
         *  to whatever example is being run.  Unset -- the normal case -- leaves the production
         *  values in place. */
        static int quadratureSetting(const char* variable, const int fallback)
        {
            const char* value= std::getenv(variable);
            if (value == nullptr) return fallback;
            const int parsed= std::atoi(value);
            return parsed > 0 ? parsed : fallback;
        }

        inline static int facetImageShells =
            quadratureSetting("OILAB_FACET_IMAGE_SHELLS", 4);

        /*! How finely each triangle is subdivided for that quadrature.
         *
         *  The triangulation only joins the nodes, so without subdivision the displacement is
         *  held constant over triangles as large as the node spacing.  This is the dominant cost
         *  of building a mesostate, so it is the first thing to try lowering -- but it changes
         *  the displacement field, so lower it only against a measured comparison. */
        inline static int facetRefinement =
            quadratureSetting("OILAB_FACET_REFINEMENT", 4);

        /*! The deformed boundary surface, as a facet in its own right.
         *
         *  facetA and facetB are the two grains' surfaces in the REFERENCE configuration, and
         *  they differ; their deformed images coincide, and that common surface is the boundary
         *  the atoms actually see once the displacement has been applied.  Only that surface can
         *  answer whether a displaced atom is still inside its own grain, so it is built here:
         *  the same shared connectivity over the deformed nodes, carrying no displacement of its
         *  own.
         *
         *  It is subdivided only once.  Refinement exists for the displacement quadrature, and
         *  nothing asks this facet for a displacement -- only which side of it a point is on,
         *  which subdividing a flat triangle cannot change.
         */
        GbFacet deformedSurface;

        /*! @param boxDim the two in-plane period vectors of the mesostate box in
         * Cartesian coordinates, laid out as
         * \f$\{p_{1x},p_{1y},p_{1z},p_{2x},p_{2y},p_{2z}\}\f$. Used by GbFacet to
         * build the periodic triangulation.
         */
        GbContinuum(const Eigen::Matrix<double, dim,dim-1>& domain,
                    const std::pair<XuPairs,XuPairs>& xuPairsOfFacetedSurfaces,
                    const VectorDimD& normalGrainA,
                    const std::vector<double>& boxDim,
                    const bool& verbosity=false);

        /*! Returns the displacement at position \f$\mathbf{x}\f$ in the region occupied
         * by grain \f$i\f$. */
        VectorDimD displacement(const VectorDimD& x, const int& i) const;

        /*! \brief Tells if the point \f$\mathbf{x}\f$ belongs to undeformed Grain A, i.e. lies
         *  below facet A.  \p normalGrainA points out of grain A, so "below" means on the side
         *  \f$-\hat n_{\mathcal A}\f$. */
        bool inGrainA(const VectorDimD& x) const;

        /*! \brief Tells if the DEFORMED position \f$\mathbf{x}\f$ is still inside grain A, i.e.
         *  on grain A's side of the deformed boundary surface.
         *
         *  inGrainA() cuts the reference crystal at the reference facet, which is the right cut
         *  to make there; but the displacement field can carry an atom that was on its own side
         *  of that facet through the boundary and out the other side, leaving it inside a crystal
         *  it does not belong to.  Nothing downstream notices: the atom is at a perfectly ordinary
         *  distance from its new neighbours, so the overlap removal keeps it and the energy simply
         *  comes back too high.  This is the test that catches it. */
        bool inGrainAAfterDeformation(const VectorDimD& x) const;

        /*! As inGrainAAfterDeformation(), for grain B.  The deformed surfaces of the two grains
         *  are one surface, so both questions are asked of it and differ only in which side
         *  counts as inside. */
        bool inGrainBAfterDeformation(const VectorDimD& x) const;

        /*! \brief Whether the DEFORMED position \f$\mathbf x\f$ lies on the boundary surface
         *  itself.
         *
         *  The two side tests above keep such a point, because the points the construction brings
         *  the grains together at sit exactly there and belong to the boundary rather than to
         *  neither side of it.  That reasoning holds for an atom that meets another one there; it
         *  does not hold for an atom the field carries onto the surface alone.  Telling the two
         *  apart needs to know what else is at that point, which is a question about the
         *  configuration rather than about the geometry, so this only answers where the point is
         *  and leaves the decision to the caller.
         *
         *  \p tolerance has to be given, and has to be the distance at which the caller already
         *  treats two atoms as one.  GbFacet::isOnSurface() defaults to 1e-6 A, which is a
         *  question about exact arithmetic, not about atoms: the field routinely leaves an atom a
         *  few thousandths of an angstrom off the surface, which is far outside 1e-6 and far
         *  inside the distance at which the relaxation fuses a pair.  Such an atom stands on the
         *  boundary for every purpose that matters, and a 1e-6 test walks straight past it. */
        bool onDeformedSurface(const VectorDimD& x, const double& tolerance) const;

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
