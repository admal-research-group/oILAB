//
// Created by Nikhil Chandra Admal on 5/17/24.
//

#ifndef OILAB_GBMESOSTATES_H
#define OILAB_GBMESOSTATES_H

#include "GbShifts.h"
#include <deque>
#include "GbMesoState.h"
#include "../MonteCarlo/Ensemble.h"

namespace oILAB {
/*! Class template that aids in the construction of an ensemble of GB
 * mesostates.
 *
 */
template <int dim>
class GbMesoStateEnsemble : public GbShifts<dim>,
                            public Ensemble<XTuplet,GbMesoState<dim>, GbMesoStateEnsemble<dim>> {
    using VectorDimD = typename LatticeCore<dim>::VectorDimD;
    using BicrystalLatticeVectors = std::vector<LatticeVector<dim>>;
    using Constraints = XTuplet;

    /*! The coincidence nodes a signature engages. Used when the ensemble was built with
     * \p GbShiftSearch::Sites, where a state is a set of nodes rather than of (t,s) pairs. */
    static std::deque<GbNode<dim>> getEngagedNodes(
        const std::vector<GbNode<dim>>& nodes,
        const Constraints& constraints);

    static std::deque<std::pair<LatticeVector<dim>, VectorDimD>> getEngagedTsPairs(
        const std::vector<std::pair<LatticeVector<dim>,VectorDimD>> &bShiftPairs,
        const Constraints &constraints);

public:
    /*!
    * CSL vectors that define the ensemble's grain boundary region
    */
    std::vector<LatticeVector<dim>> ensembleCslVectors;

    /*! @param search \p GbShiftSearch::Flat (default) keeps the original enumeration; \p Full
     *         also admits translations whose CSL shift leaves the boundary plane, which is what
     *         makes non-flat mesostates reachable.  The remaining arguments are forwarded to
     *         GbShifts and ignored unless \p search is \p Full.
     */
    GbMesoStateEnsemble(const Gb<dim> &gb,
                      const ReciprocalLatticeVector<dim> &axis,
                      std::vector<LatticeVector<dim>> &ensembleCslVectors,
                      const double &tMax=1,
                      const double& sPerpMax=1,
                      const GbShiftSearch& search= GbShiftSearch::Flat,
                      const double& tPerpMax= 1.0e300,
                      const bool& oneTranslationPerSite= false,
                      const std::string& filename= "translationsNonFlat.txt",
                      const double& slabHalfThickness= 1.0,
                      const double& dMax= 1.5,
                      const bool& dropInvertedNodes= true,
                      const bool& dropZeroJumpNodes= true);

    /*!
    * \brief Constructs an ensemble of mesostates
    * @param filename-
    * @return A deque of mesostates
    */
    std::map<Constraints, GbMesoState<dim>>
    collectMesoStates(const std::string &filename = "") const;

    static std::deque<Constraints> enumerateConstraints(const int& size);

    /*!
    * \brief Evove mesostates using a Monte Carlo algorithm
    * @param filename-
    * @return A deque of mesostates
    */
    GbMesoState<dim> constructMesoState(const Constraints &constraints) const;

    Constraints sampleNewState(const Constraints &currentConstraints,
                             const bool &randomize = false) const;

    Constraints initializeState() const;

};

} // namespace oILAB

#include "GbMesoStateEnsembleImplementation.h"
#endif //OILAB_GBMESOSTATES_H
