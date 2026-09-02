//
// Created by Nikhil Chandra Admal on 12/11/23.
//

#ifndef OILAB_MESOSTATE_H
#define OILAB_MESOSTATE_H

#include "Gb.h"
#include <string>

namespace oILAB {

/*! \brief Which set of (translation, shift) pairs GbShifts enumerates.
 *
 *  \p Flat is the original scheme and the default, so existing behaviour is unchanged: DSCL
 *  translations are taken from a parallelepiped whose edges follow the CSL box directions, which
 *  restricts attention to translations whose CSL shifts lie along the boundary.
 *
 *  \p Full drops that restriction. Translations are taken from a ball \f$|t|\le t_{max}\f$
 *  intersected with a slab \f$|t\cdot\hat n|\le t_{\perp max}\f$ about the boundary, and each one
 *  yields a single shift \f$s\f$. Because the CSL shift is no longer forced into the boundary
 *  plane, the resulting mesostates can be non-flat.
 */
enum class GbShiftSearch
{
    Flat,
    Full
};

template <int dim> class GbShifts {
  using VectorDimD = typename LatticeCore<dim>::VectorDimD;
  using VectorDimI = typename LatticeCore<dim>::VectorDimI;

protected:
  // static std::vector<LatticeVector<dim>> getGbCslVectors(const Gb<dim>& gb,
  // const ReciprocalLatticeVector<dim>& axis);
  static std::vector<std::pair<LatticeVector<dim>, VectorDimD>>
  getbShiftPairs(const Gb<dim> &gb,
                 const std::vector<LatticeVector<dim>> &gbCslVectors,
                 const double& tMax,
                 const double& sPerpMax);

  /*! \brief Enumerates (translation, shift) pairs without restricting the CSL shift to the
   *  boundary plane, so that non-flat mesostates become reachable.
   *
   *  For every DSCL translation \f$t\f$ with \f$|t|\le t_{max}b\f$ and
   *  \f$|t\cdot\hat n|\le t_{\perp max}b\f$, the shift is
   *  \f[ s = \left(\Lambda_{\mathcal A}t - t/2\right) \bmod \mathcal C_{\text{mesostate}} \f]
   *  and the pair is kept when \f$|s\cdot\hat n|\le s_{\perp max}b/2\f$.  All lengths scale with
   *  \f$b\f$, the shortest lattice vector of \f$\mathcal A\f$.
   *
   *  The reduction is taken modulo the *mesostate* CSL cell, not the primitive one: reducing
   *  modulo the primitive cell confines every shift to a fraction of the box, leaving most of the
   *  boundary without nodes.  The whole expression is reduced, rather than just
   *  \f$\Lambda_{\mathcal A}t\f$, so that every shift lands strictly inside the cell.
   *
   *  \f$x_{\mathcal A}=s-t/2\f$ and \f$x_{\mathcal B}=s+t/2\f$ are lattice vectors of
   *  \f$\mathcal A\f$ and \f$\mathcal B\f$ for *any* DSCL translation, because
   *  \f$\Lambda_{\mathcal B}t\in\mathcal A\f$ and \f$\Lambda_{\mathcal A}t\in\mathcal B\f$ by
   *  the definition of the shift tensors, and \f$\Lambda_{\mathcal A}+\Lambda_{\mathcal B}=I\f$.
   *  So widening the search cannot break the lattice-membership requirement downstream.
   *
   *  @param tPerpMax half-thickness, in units of \f$b\f$, of the slab about the boundary that the
   *         translations are drawn from.  Values \f$\ge t_{max}\f$ leave the ball untouched, since
   *         \f$|t|\le t_{max}b\f$ already implies \f$|t\cdot\hat n|\le t_{max}b\f$.
   *  @param oneTranslationPerSite when true, keeps only the shortest translation among those that
   *         share a shift.  Distinct translations can fold onto the same shift, and two such pairs
   *         engaged together put two nodes at one point, which GbMesoState rejects as a clash.
   *         Leaving this false keeps the larger set and lets those signatures be rejected
   *         downstream instead.
   *  @param filename if non-empty, the surviving pairs are written there, one per line.
   */
  static std::vector<std::pair<LatticeVector<dim>, VectorDimD>>
  getNonFlatShiftPairs(const Gb<dim> &gb,
                       const std::vector<LatticeVector<dim>> &gbCslVectors,
                       const double& tMax,
                       const double& sPerpMax,
                       const double& tPerpMax,
                       const bool& oneTranslationPerSite,
                       const std::string& filename);

public:
  const Gb<dim> &gb;
  const ReciprocalLatticeVector<dim> &axis;
  const std::vector<LatticeVector<dim>> gbCslVectors;
  std::vector<std::pair<LatticeVector<dim>, VectorDimD>> tShiftPairs;
  /*! Which enumeration produced \p tShiftPairs. */
  const GbShiftSearch search;

  /*! @param search \p GbShiftSearch::Flat (default) reproduces the original enumeration exactly;
   *         \p GbShiftSearch::Full uses getNonFlatShiftPairs() instead.
   *  @param tPerpMax slab half-thickness in units of \f$b\f$; ignored unless \p search is
   *         \p Full.  The default is larger than any sensible \p tMax, so the slab is inert
   *         until asked for.
   *  @param oneTranslationPerSite ignored unless \p search is \p Full; see
   *         getNonFlatShiftPairs().
   *  @param filename ignored unless \p search is \p Full; where to write the pair list.
   */
  explicit GbShifts(const Gb<dim> &gb, const ReciprocalLatticeVector<dim> &axis,
                    const std::vector<LatticeVector<dim>> &gbCslVectors,
                    const double& tMax= 1,
                    const double& sPerpMax= 1,
                    const GbShiftSearch& search= GbShiftSearch::Flat,
                    const double& tPerpMax= 1.0e300,
                    const bool& oneTranslationPerSite= false,
                    const std::string& filename= "translationsNonFlat.txt");

    };
    } // namespace oILAB
#endif //OILAB_MESOSTATE_H
