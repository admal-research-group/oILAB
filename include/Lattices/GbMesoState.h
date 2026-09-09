//
// Created by Nikhil Chandra Admal on 5/26/24.
//

#ifndef OILAB_GBMESOSTATE_H
#define OILAB_GBMESOSTATE_H

#include "Gb.h"
#include "GbContinuum.h"
#include "LatticeCore.h"
#include "OrderedTuplet.h"
#include "../Math/PeriodicFunction.h"

namespace oILAB {

/*! Class template that defines a GB mesostate.
 *
 */
template <int dim> class GbMesoState : public GbContinuum<dim> {
    using VectorDimD = LatticeCore<dim>::VectorDimD;
    using XuPairs = std::deque<std::pair<VectorDimD,VectorDimD>>;

  /*!
   * \brief Returns the cartesian coordinates of the CSL vectors that define a
   * mesostate's GB.
   * @param mesoStateCslVectors - a vector (size = \p dim) of CSL vectors that
   * define the box of the mesostate.
   * @return Cartesian coordinates of the \p dim-1 grain boundary CSL vectors
   */
  // ensure that the input is of the right dimension
  static Eigen::Matrix<double, dim, dim - 1> getMesoStateGbDomain(const std::vector<LatticeVector<dim>> &mesoStateCslVectors);

  /*!
   * \brief Returns the nodes \f$\textbf x\f$ and their displacements \f$\textbf u\f$ of the faceted boundaries of the
   * grains that form the GB.
   * @param gb - grain boundary
   * @param mesoStateCslVectors - the box vectors of the mesostate
   * @param engagedTsPairs - a deque of engaged (t,s) pairs (translation vector \p t, shift vector \p s).
   * \p t is a DSCL vector while \p s is expressed in Cartesian coordinates
   * @return A pair of deques for the boundaries of the two undeformed grains that are deformed and glued to form the GB.
   * Each deque contains `(x,u)` pairs.
   */
  static std::pair<XuPairs,XuPairs> getFacetedSurfaces(const Gb<dim> &gb,
      const std::vector<LatticeVector<dim>> &mesoStateCslVectors,
      const std::deque<std::pair<LatticeVector<dim>, VectorDimD>> &engagedTsPairs);

  /*! \brief The same, for nodes whose two grains are displaced independently.
   *
   * Each node carries its own \f$\textbf x_{\mathcal A}\f$, \f$\textbf x_{\mathcal B}\f$
   * and coincidence point, so the displacements \f$\textbf u_{\mathcal A}=\textbf s-\textbf
   * x_{\mathcal A}\f$ and \f$\textbf u_{\mathcal B}=\textbf s-\textbf x_{\mathcal B}\f$
   * are read off rather than split evenly. Both still carry their grain onto \f$\textbf s\f$,
   * which is all the gluing requires.
   */
  static std::pair<XuPairs,XuPairs> getFacetedSurfaces(const Gb<dim> &gb,
      const std::vector<LatticeVector<dim>> &mesoStateCslVectors,
      const std::deque<GbNode<dim>> &engagedNodes);

  /*!
   * \brief Returns the two in-plane period vectors of the mesostate box in Cartesian
   * coordinates, laid out as \f$\{p_{1x},p_{1y},p_{1z},p_{2x},p_{2y},p_{2z}\}\f$, as
   * expected by GbContinuum/GbFacet.
   * @param cslVectors the box vectors of the mesostate
   */
  static std::vector<double> getMesoStateBoxDim(const std::vector<LatticeVector<dim>> &cslVectors);

public:
  /*! \brief Where one \f$(\textbf t,\textbf s)\f$ pair puts its two nodes.
   *
   *  \p xA and \p xB are the node positions of the two grains, wrapped into the bicrystal box,
   *  and \p u is half the translation -- grain A carries \f$+\textbf u\f$ and grain B
   *  \f$-\textbf u\f$, so both nodes deform onto the CSL shift \f$\textbf s\f$.  \p siteA
   *  and \p siteB are the same two positions as integer coordinates of lattices
   *  \f$\mathcal A\f$ and \f$\mathcal B\f$, which is what decides whether two pairs collide.
   */
  struct NodePlacement {
    VectorDimD xA;
    VectorDimD xB;
    VectorDimD u;
    OrderedTuplet<dim> siteA;
    OrderedTuplet<dim> siteB;
  };

  /*!
   * \brief Returns the two nodes a single \f$(\textbf t,\textbf s)\f$ pair contributes.
   *
   * A pair depends on nothing but itself, so its two sites can be worked out before any
   * mesostate is built.  Two engaged pairs collide -- and the mesostate is rejected -- exactly
   * when they share \p siteA or share \p siteB, so a caller that knows the sites can avoid
   * enumerating the colliding combinations rather than constructing and discarding them.
   * getFacetedSurfaces() uses this same function, so the two can never disagree.
   *
   * @param gb grain boundary
   * @param mesoStateCslVectors the box vectors of the mesostate
   * @param t translation vector (a DSCL vector)
   * @param s shift vector, in Cartesian coordinates
   * @throws std::runtime_error if either node fails to land on its lattice
   */
  static NodePlacement nodePlacement(const Gb<dim> &gb,
                                     const std::vector<LatticeVector<dim>> &mesoStateCslVectors,
                                     const LatticeVector<dim> &t,
                                     const VectorDimD &s);

  /*! Species written for the atoms the mesostate brings into coincidence, so that the boundary
   *  the construction built can be picked out of the relaxed structure.  Grain \f$\mathcal A\f$
   *  is 1 and \f$\mathcal B\f$ is 2; both members of a coincident pair carry this instead, so
   *  whichever of the two survives the overlap removal still identifies the boundary plane. */
  static constexpr int coincidenceType = 3;

  /*!
   * Grain boundary
   */
  const Gb<dim> &gb;

  /*!
   * Grain boundary tilt axis
   */
  const ReciprocalLatticeVector<dim> &axis;

  /*!
   * A vector (size = \p dim) of CSL vectors that define the box of the
   * mesostate. The second and third vectors should be parallel to the grain
   * boundary, while the third vector should be out of the grain boundary plane.
   */
  const std::vector<LatticeVector<dim>> &mesoStateCslVectors;

  /*!
   * @param engagedTsPairs a deque of pairs (translation vector \f$\textbf t\f$, shift
   * vector \f$\textbf s\f$) that defines a mesostate. Translating lattice
    * \f$\mathcal A\f$ by \f$\textbf t/2\f$ and lattice
     * \f$\mathcal B\f$ by \f$-\textbf t/2\f$ results in a CSL shift of \f$\textbf
   * s\f$.
   */
  const std::deque<std::pair<LatticeVector<dim>, VectorDimD>> engagedTsPairs;

  /*!
   * The engaged coincidence nodes, when the mesostate was built from them. Empty for a
   * mesostate built from \p engagedTsPairs, and vice versa.
   */
  const std::deque<GbNode<dim>> engagedNodes;

  explicit GbMesoState(
      const Gb<dim> &gb,
      const ReciprocalLatticeVector<dim> &axis,
      const std::deque<std::pair<LatticeVector<dim>, VectorDimD>>& engagedTsPairs,
      const std::vector<LatticeVector<dim>> &mesoStateCslVectors);

  /*! Builds the mesostate from coincidence nodes, so that the two grains may be displaced by
   * different amounts at each node. */
  explicit GbMesoState(
      const Gb<dim> &gb,
      const ReciprocalLatticeVector<dim> &axis,
      const std::deque<GbNode<dim>>& engagedNodes,
      const std::vector<LatticeVector<dim>> &mesoStateCslVectors);

  /*!
   * \brief Calculate the energy of a mesostate using lammps
   * @param minimize when true the configuration is relaxed in LAMMPS before its energy is read,
   * so the returned energy is that of the minimized configuration.  False (the default) reports
   * the energy of the as-constructed configuration.
   * @param configFile a deformed configuration already written by box(). Writing one is by far
   * the most expensive step of building a mesostate -- the displacement field is evaluated at
   * every atom -- so a caller that has written the configuration for its own output should hand
   * the path over rather than have it computed a second time. Empty (the default) writes a
   * scratch copy, which is the behaviour when the caller has none.
   * @param minimizedDumpFile if non-empty, where to write the configuration as LAMMPS leaves it.
   * @param tetherHalfWidth when positive, atoms within this distance of the boundary are held to
   * their as-constructed positions by a harmonic spring during the relaxation.
   * @param tetherStiffness spring constant of that restraint, in eV/Angstrom^2.
   * @param springEnergy if non-null, receives the energy stored in the restraint.
   * @param unminimizedEnergy if non-null, receives the boundary energy before relaxation.
   * @return (density, energy) of the mesostate
   */
  /*! \brief What both relaxations of one mesostate cost.
   *
   *  A tether answers "what does the state the enumeration built cost", by holding the boundary
   *  atoms where the construction put them; a free minimisation answers "what does the boundary
   *  this state leads to cost", by letting them go.  The two are different questions and a
   *  faceting study wants both, so both are reported, together with the energy before either
   *  relaxation. */
  /*! \brief A mesostate's atoms in memory, laid out exactly as read_oILAB_output() returns
   *  them: \p atoms one row per atom -- species, x, y, z, radius -- \p box holding the three
   *  cell vectors as its columns, and \p origin the cell origin.
   *
   *  box() used to hand its result on only as a file, which the coincidence count then re-read
   *  and LAMMPS read again.  LAMMPS takes its atoms in memory now, so the file is needed only
   *  when a state is being kept to look at, and a survey pass can skip writing it entirely --
   *  measured at 23% of box().  Matching the reader's layout means nothing downstream has to
   *  care which way the atoms arrived. */
  struct Configuration
  {
    Eigen::MatrixXd atoms;
    Eigen::Matrix3d box;
    Eigen::Vector3d origin;
  };

  struct Relaxations
  {
    double density   = 0.0;   //!< atoms in the boundary region; the same for both relaxations,
                              //!< since the overlap removal and the group definitions precede them
    double unrelaxed = 0.0;   //!< the configuration as constructed, before any relaxation
    double tethered  = 0.0;   //!< relaxed with the boundary atoms restrained
    double spring    = 0.0;   //!< energy the restraint had to store to hold them
    double full      = 0.0;   //!< relaxed with nothing held
  };

  /*! \brief Run both relaxations of this mesostate, each from the as-constructed configuration.
   *
   *  The two are separate LAMMPS runs over the same input, not one run continued: a free
   *  minimisation started from the tethered result would be exploring the neighbourhood of the
   *  tethered structure rather than of the state itself, and the two energies would no longer be
   *  comparable.  With \p tetherHalfWidth at zero there is only one relaxation to run and the
   *  tethered figures repeat the free ones.
   *
   *  \p chainRelaxations does the two in one invocation instead: relax tethered, release the
   *  restraint, relax again from there.  That is one LAMMPS start-up rather than two and skips
   *  the descent the free run would repeat, at the cost of asking a different question -- the
   *  free minimum reached from the tethered structure need not be the one reached from the
   *  as-constructed structure.  Measure the difference on the boundary at hand before trusting
   *  it; on some landscapes it is nothing and on others it is not.
   *
   *  @param configFile a deformed configuration already written by box(); empty writes a scratch
   *         copy.  Writing one evaluates the displacement field at every atom and dominates the
   *         cost of a state, so a caller that has one should hand it over.
   *  @param tetheredDumpFile where to leave the tethered structure, if anywhere.
   *  @param fullDumpFile where to leave the freely relaxed structure, if anywhere. */
  Relaxations relaxations(const std::string &lmpLocation,
                          const std::string &potentialName,
                          const std::string &configFile = "",
                          const double &tetherHalfWidth = 0.0,
                          const double &tetherStiffness = 1.0,
                          const std::string &tetheredDumpFile = "",
                          const std::string &fullDumpFile = "",
                          const bool &chainRelaxations = false) const;

  /*! As above, from a configuration already in memory -- what box() hands back when it is asked
   *  for one.  This is the form the sweep uses; the overload above reads a file and delegates
   *  here, and exists for callers that still have only a path. */
  Relaxations relaxations(const std::string &lmpLocation,
                          const std::string &potentialName,
                          const Configuration &configuration,
                          const double &tetherHalfWidth,
                          const double &tetherStiffness,
                          const std::string &tetheredDumpFile,
                          const std::string &fullDumpFile,
                          const bool &chainRelaxations) const;

  // std::pair<double,double> densityEnergy() const;
  std::tuple<double, double>
  densityEnergy(const std::string &lmpLocation,
                const std::string &potentialName,
                const bool &minimize = false,
                const std::string &configFile = "",
                const std::string &minimizedDumpFile = "",
                const double &tetherHalfWidth = 0.0,
                const double &tetherStiffness = 1.0,
                double *springEnergy = nullptr,
                double *unminimizedEnergy = nullptr) const;

  /*! This function outputs/prints a grain boundary mesostate
   * @param filename name of the file to be written to
   * @param atomsExpelled if non-null, receives the number of atoms the deformation carried out
   * of their own grain and which were therefore left out of both configurations.
   * @param dropUnengagedCoincidences when true, every coincident group the state did not engage
   * is removed entirely -- both of the atoms that met there -- so that the configuration holds
   * exactly the coincidences the signature names and an unengaged site is empty.  This changes
   * the structure: it leaves a vacancy where the pair was, so the atom count, the density and
   * the relaxed structure differ from the same state built untouched, and energies from a run
   * with this on are not comparable with energies from a run with it off.
   * @param atomsDropped if non-null, receives how many atoms that removed.
   */
  typename std::enable_if<dim == 3, void>::type
  /*! @param filename where to write the two extended-XYZ files, \p _reference0 for the
   *         undeformed configuration and \p _reference1 for the deformed one.  EMPTY writes
   *         neither, which is what a survey pass wants: the files cost about a quarter of this
   *         function and nothing reads them unless the state is being kept.
   *  @param deformedConfiguration if non-null, receives the deformed configuration in memory,
   *         ready for the coincidence count and for LAMMPS without a file in between. */
  box(const std::string &filename, int* atomsExpelled = nullptr,
      const bool& dropUnengagedCoincidences = false, int* atomsDropped = nullptr,
      Configuration* deformedConfiguration = nullptr) const;
    };
    } // namespace oILAB

#include "GbMesoStateImplementation.h"
#endif //OILAB_GBMESOSTATE_H
