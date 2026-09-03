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
    Full,
    Sites
};

/*! \brief One coincidence node: an atom of \f$\mathcal A\f$ and an atom of
 *  \f$\mathcal B\f$ brought together at a point.
 *
 *  The two grains are displaced independently, so the coincidence point need not be the
 *  midpoint of the pair that meets there -- which is what the \p Flat and \p Full searches
 *  assume, and what makes the sites off the boundary plane unreachable to them.  The jump
 *  \f$\textbf t=\textbf x_{\mathcal B}-\textbf x_{\mathcal A}\f$ is a DSCL vector either
 *  way; only its division between the grains differs.
 */
template <int dim>
struct GbNode
{
    using VectorDimD = typename LatticeCore<dim>::VectorDimD;

    /*! The coincidence point.  A point of \f$\tfrac{1}{2}\mathcal D\f$, the lattice of
     *  midpoints of \f$\mathcal A\f$ and \f$\mathcal B\f$ atoms. */
    VectorDimD s;
    /*! The atom of \f$\mathcal A\f$ brought to \p s. */
    LatticeVector<dim> xA;
    /*! The atom of \f$\mathcal B\f$ brought to \p s. */
    LatticeVector<dim> xB;

    GbNode(const VectorDimD& s_, const LatticeVector<dim>& xA_, const LatticeVector<dim>& xB_) :
        s(s_), xA(xA_), xB(xB_) {}

    /*! Displacement of grain \f$\mathcal A\f$ at this node. */
    VectorDimD uA() const { return s - xA.cartesian(); }
    /*! Displacement of grain \f$\mathcal B\f$ at this node.  Not \f$-\textbf u_{\mathcal A}\f$
     *  in general. */
    VectorDimD uB() const { return s - xB.cartesian(); }
    /*! The jump across the boundary, \f$\textbf u_{\mathcal A}-\textbf u_{\mathcal B}\f$. */
    VectorDimD t()  const { return xB.cartesian() - xA.cartesian(); }
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

  /*! \brief Enumerates coincidence nodes site-first: every point of
   *  \f$\tfrac{1}{2}\mathcal D\f$ within \p slabHalfThickness of the flat boundary, paired
   *  with every atom of \f$\mathcal A\f$ and of \f$\mathcal B\f$ that can reach it by moving
   *  no further than \p dMax.
   *
   *  This inverts the order the other two searches use.  They enumerate translations and derive
   *  the coincidence point as a midpoint, which confines it to \f$\tfrac{1}{2}\mathcal D\f$
   *  points that happen to be midpoints of an admissible pair, and -- because the shift is then
   *  reduced into a cell a whole number of CSL planes thick -- collapses to the boundary plane
   *  itself whenever the requested layer is thinner than one such plane.  Choosing the point
   *  first and letting the two grains move by different amounts removes both restrictions.
   *
   *  @param slabHalfThickness how far off the flat boundary, in Angstrom, sites are taken from
   *  @param dMax how far one atom may move, in Angstrom, to reach its site
   *  @param dropInvertedNodes discard nodes whose atom of \f$\mathcal B\f$ lies below its atom
   *         of \f$\mathcal A\f$ along the boundary normal, i.e. \f$\textbf t\cdot\hat n<0\f$.
   *         The two grains occupy fixed sides, so such a node turns the boundary inside out
   *         locally: box() keeps an atom of \f$\mathcal A\f$ below facet A and one of
   *         \f$\mathcal B\f$ above facet B, and with the facets swapped the slab between them
   *         is claimed by both grains and each keeps material well past where the other begins.
   *         Every discarded node has a mirror partner with the two atoms exchanged, so this
   *         removes wrong-handed duplicates rather than reachable states.
   *  @param dropZeroJumpNodes discard nodes whose two atoms are already the same point.  Their
   *         jump vanishes, so both grains are displaced identically and the boundary is left
   *         undeformed; engaged alone each is simply the flat boundary again.  They remain
   *         meaningful in a state that engages several nodes, where such a node pins the surface
   *         at a point the two grains already share, so this is a flag rather than a deletion.
   *  @param filename if non-empty, the nodes are written there
   */
  static std::vector<GbNode<dim>>
  getSiteNodes(const Gb<dim> &gb,
               const std::vector<LatticeVector<dim>> &gbCslVectors,
               const double& slabHalfThickness,
               const double& dMax,
               const bool& dropInvertedNodes,
               const bool& dropZeroJumpNodes,
               const std::string& filename);

public:
  const Gb<dim> &gb;
  const ReciprocalLatticeVector<dim> &axis;
  const std::vector<LatticeVector<dim>> gbCslVectors;
  std::vector<std::pair<LatticeVector<dim>, VectorDimD>> tShiftPairs;
  /*! The coincidence nodes.  Populated only by \p GbShiftSearch::Sites; the other two searches
   *  describe their states through \p tShiftPairs instead. */
  std::vector<GbNode<dim>> nodes;
  /*! Which enumeration produced \p tShiftPairs or \p nodes. */
  const GbShiftSearch search;

  /*! @param search \p GbShiftSearch::Flat (default) reproduces the original enumeration exactly;
   *         \p GbShiftSearch::Full uses getNonFlatShiftPairs() instead.
   *  @param tPerpMax slab half-thickness in units of \f$b\f$; ignored unless \p search is
   *         \p Full.  The default is larger than any sensible \p tMax, so the slab is inert
   *         until asked for.
   *  @param oneTranslationPerSite ignored unless \p search is \p Full; see
   *         getNonFlatShiftPairs().
   *  @param filename where to write the pair or node list; ignored when \p search is \p Flat.
   *  @param slabHalfThickness used only by \p Sites: how far off the boundary, in Angstrom,
   *         coincidence points are taken from.
   *  @param dMax used only by \p Sites: how far one atom may move, in Angstrom, to reach its
   *         coincidence point.
   *  @param dropInvertedNodes used only by \p Sites; see getSiteNodes().
   *  @param dropZeroJumpNodes used only by \p Sites; see getSiteNodes().
   */
  explicit GbShifts(const Gb<dim> &gb, const ReciprocalLatticeVector<dim> &axis,
                    const std::vector<LatticeVector<dim>> &gbCslVectors,
                    const double& tMax= 1,
                    const double& sPerpMax= 1,
                    const GbShiftSearch& search= GbShiftSearch::Flat,
                    const double& tPerpMax= 1.0e300,
                    const bool& oneTranslationPerSite= false,
                    const std::string& filename= "translationsNonFlat.txt",
                    const double& slabHalfThickness= 1.0,
                    const double& dMax= 1.5,
                    const bool& dropInvertedNodes= true,
                    const bool& dropZeroJumpNodes= true);

    };
    } // namespace oILAB
#endif //OILAB_MESOSTATE_H
