//
// Created by himanshu on 9/1/26.
//

#ifndef OILAB_GBFACET_H
#define OILAB_GBFACET_H

#include <array>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_triangle_primitive_3.h>

/*! \brief A triangulated, doubly-periodic faceted surface carrying a nodal displacement field.
 *
 *  A GbFacet is built from a cloud of nodes \f$\{(\mathbf x_i,\mathbf u_i)\}\f$ together with an
 *  externally supplied periodic connectivity (\p Topology).  Taking the connectivity from
 *  outside is what allows the two facets of a grain boundary to share it, and sharing it is
 *  what glues the boundary together: the deformed nodes \f$\mathbf x_i+\mathbf u_i\f$ of the
 *  two grains coincide by construction, so if both facets are triangulated with the *same*
 *  connectivity their deformed surfaces are identical triangle-for-triangle, i.e. one
 *  continuous surface.  Triangulating the two clouds independently does *not* achieve this:
 *  the two reference clouds differ (by the translation vectors \f$\mathbf t_i\f$), so
 *  independent Delaunay triangulations produce different connectivity and the deformed
 *  surfaces then agree only at the nodes.
 */
class GbFacet
{
public:

    /*! Whether a facet announces itself on stdout when it is built.
     *
     *  Off by default.  A sweep builds three facets per mesostate -- the two grains' and the
     *  deformed surface -- so over a run of any size the announcements are tens of thousands of
     *  lines, interleaved from every thread, and they bury whatever the run was actually
     *  reporting.  Turn it on for a single state, where the mesh is the thing being looked at. */
    static bool announceConstruction;

    /*! \brief Whether the solid angles are memoised across states.
     *
     *  The panels the quadrature integrates over are built from the REFERENCE node positions --
     *  only their displacements come from the state -- so every solid angle in the image sum is
     *  pure geometry, and the field is linear in the nodal displacements.  That makes the whole
     *  image sum for one (query point, triangle) pair reusable by any other state that puts the
     *  same triangle at the same place, which sweeps do constantly: they draw their nodes from
     *  one fixed pool of coincidence sites and evaluate at one fixed set of lattice positions.
     *
     *  Measured on sigma5 (310).  Over 1119 states, 83% of lookups hit and box() falls from 83.0
     *  to 43.7 ms/state; over 26,879 states, 97.7% hit -- a reuse factor of 43.6 -- and box()
     *  falls from 109.2 to 42.5 ms/state, 2.6x, for 2.3 million entries and about 350 MB.  The
     *  reuse grows with the sweep because the pool of distinct triangles is bounded while the
     *  number of states is not, so the bigger the run the more this returns.  The image sum is
     *  96% of displacement(), which is why halving it halves the field.
     *
     *  THAT BOUNDEDNESS IS A PROPERTY OF FLAT BOUNDARIES, NOT OF THE METHOD.  It holds while
     *  the states of a sweep share one shape: the whole 26,879-state flat sweep interned 2,144
     *  triangles, 0.08 per state.  It fails as soon as the shape varies, because the faceted
     *  surface is triangulated THROUGH the engaged coincidence points -- so shape and
     *  engagement are one variable, and every state brings its own triangulation.  Measured on
     *  the faceted sigma5 (310) sweep at slabHalfThickness 0.4: 5.3 new triangles per state,
     *  fitting triangles = 18.7*states^0.90, which is 39 million of them and some 4 GB by the
     *  end of a 10.7M-state run -- sixty-six times the per-state growth of the flat case.
     *
     *  Hence the cap on triangleInterner, and hence facetedGbTreeStrategy.txt: choosing a shape
     *  first and enumerating the coincidence points that lie on it would let the states of a
     *  branch share one triangulation, which is what would restore the premise this paragraph
     *  opens with.
     *
     *  It is exact, not an approximation: over both sweeps every state came out byte-identical
     *  to the direct sum, all four energies included.  Set OILAB_FACET_MEMO=0 to select the
     *  direct sum, which is what that was checked against, and OILAB_FACET_MEMO_VERIFY=1 to
     *  recompute on every hit and report any disagreement. */
    inline static const bool memoiseSolidAngles =
        std::getenv("OILAB_FACET_MEMO") == nullptr
        || std::string(std::getenv("OILAB_FACET_MEMO")) != "0";

    /*! \brief Empties the solid-angle memo if the box it was filled for has changed.
     *
     *  Only the periods are checked here.  The quadrature settings differ between the facets of
     *  a single mesostate -- deformedSurface uses refinement 1 -- so they sit in the memo's key
     *  instead; resetting on them empties the table several times per state.  A no-op when
     *  nothing has changed, and called once per facet rather than once per query. */
    static void resetSolidAngleMemo(const Eigen::Vector3d& period1,
                                    const Eigen::Vector3d& period2);

    /*! Reports the memo's hit rate and footprint on exit.  Set OILAB_FACET_MEMO_STATS=1. */
    inline static const bool reportMemoStatistics =
        std::getenv("OILAB_FACET_MEMO_STATS") != nullptr;

    /*! \brief One line on what the memo is doing, for a sweep to print beside its progress.
     *
     *  The exit report alone cannot answer the question a slowing sweep actually asks, which is
     *  not "what was the hit rate overall" but "what is it doing NOW".  A run whose rate has
     *  halved has a cumulative hit rate dominated by the fast early states, so the interval
     *  figures here -- everything "since last" -- are the ones that carry the signal.
     *
     *  Four things, because between them they separate the candidate explanations for a sweep
     *  that decays: the interval hit rate (the pool of distinct geometry outgrowing the table),
     *  entries against capacity (saturation), the share of lookups still taking a shard's writer
     *  lock (insert contention, which persists only while the table has room), and the interner
     *  sizes (which are bounded by nothing and grow for the whole sweep).
     *
     *  Cheap: a handful of relaxed atomic loads plus one shared lock per interner.  Call it from
     *  a progress line, not from the hot path. */
    static std::string memoStatistics();

    /*! \brief How far away, in face radii, an image has to be before its solid angle is taken
     *  as a point dipole instead of the exact closed form.
     *
     *  Near the field point the exact form is the only thing that will do -- the solid angle is
     *  strongly peaked, and the dipole is an expansion in (radius/distance)^2.  Far away the
     *  expansion is excellent and costs a dot product instead of three square roots and an
     *  atan2, and it needs no subdivision.
     *
     *  This is not the same thing as subtracting the dipole and restoring the tail in closed
     *  form.  That was tried and it fails: the analytic total 2*pi*meanDispProjected is the
     *  value of the dipole lattice sum only when every image is far from the field point, and
     *  the L = 0 term goes as 1/h^2 for an atom sitting on the sheet.  It moved atoms by up to
     *  6.7 A.  Replacing the kernel where it is accurate is safe; restoring a far-field constant
     *  the near field never reaches is not.
     *
     *  Sixteen radii, chosen by measurement.  At K = 53 on one state the exact kernel gives a
     *  maximum displacement error of 2.298e-4 A in 1887 ms and this gives 2.334e-4 A in 766 ms
     *  -- indistinguishable, and 2.5 times faster.  Eight radii is too aggressive: it holds up
     *  to K = 32 but puts a floor near 6e-4 A that grows with K, since ever more of the sum is
     *  being approximated.  Four and two are worse still, 3.7e-3 and 1.3e-2.
     *
     *  The point of it is reach rather than speed: the image sum converges as 1/K -- measured
     *  e(K) ~ 0.05 to 0.07 A / K depending on the state -- so 10^-3 A needs K around 70, and
     *  cheap distant images are what makes that affordable.  Zero or negative uses the exact
     *  form everywhere. */
    inline static double dipoleRadii = [] {
        const char* v= std::getenv("OILAB_FACET_DIPOLE_RADII");
        return v ? std::atof(v) : 16.0;
    }();

    /*! \brief What one triangle contributes to one query point, summed over every periodic image.
     *
     *  \p w are the coefficients of the triangle's three nodal displacements -- the field being
     *  linear in them -- and \p omega the bare solid angle the same sum accumulates, which the
     *  far-field closure needs so it can subtract exactly what the explicit shells covered.
     *  Geometry only: nothing here depends on the state. */
    struct TriangleWeights
    {
        std::array<double,3> w{{0.0,0.0,0.0}};
        double omega= 0.0;
    };

    /*! One corner of a triangle.  \p node indexes the point cloud, and \p o1,\p o2 are
     * periodic image offsets: the corner sits at
     * \f$\mathbf x_{node}+o_1\mathbf p_1+o_2\mathbf p_2\f$.
     */
    struct Corner
    {
        std::size_t node;
        int o1;
        int o2;
    };

    using Triangle = std::array<Corner,3>;

    /*! Periodic connectivity: a list of consistently oriented triangles. */
    using Topology = std::vector<Triangle>;

    /*! \brief Periodic Delaunay triangulation of \p points, projected into the plane normal
     *  to \p normal, on the torus spanned by the two in-plane period vectors.
     *
     *  Evaluate this on the DEFORMED nodes \f$\mathbf x_i+\mathbf u_i\f$ and hand the result
     *  to both facets of the boundary.  Triangulating in the deformed configuration is what
     *  makes the glued surface well shaped, and it fixes the orientation: the returned
     *  triangles are oriented so that their normals point along \p normal in the deformed
     *  configuration.
     *
     *  @param points nodes to triangulate, in Cartesian coordinates
     *  @param box_dim the two in-plane period vectors, laid out as
     *         \f$\{p_{1x},p_{1y},p_{1z},p_{2x},p_{2y},p_{2z}\}\f$
     *  @param normal normal of the (nominally flat) GB plane
     *  @return \f$2n\f$ consistently oriented triangles, one per abstract torus face
     */
    static Topology triangulate(const std::vector<Eigen::Vector3d>& points,
                                const std::vector<double>& box_dim,
                                const Eigen::Vector3d& normal);

    /*! @param point_cloud one row per node, holding \f$\{x,y,z,u_x,u_y,u_z\}\f$
     *  @param box_dim the two in-plane period vectors, as in triangulate()
     *  @param topology connectivity, normally shared with the opposing facet
     *  @param normalSense +1 to keep the orientation \p topology was built with, -1 to reverse
     *         it.  The two facets of a boundary share a topology but face opposite ways: each
     *         grain must sit on the positive side of its own facet, otherwise the solid angle,
     *         and with it the sign of the displacement, comes out backwards for that grain.
     *  @param imageShells how many shells of periodic images to integrate explicitly; the
     *         remainder is closed off analytically.  The closure absorbs the *mean* nodal
     *         displacement exactly, so this only controls how well the variation about that mean
     *         is resolved far from the facet.  That part converges as 1/imageShells -- each distant
     *         image contributes only its dipole moment -- and the cost grows as (2n+1)^2, so 4 is a
     *         compromise: it leaves a floor of order 1e-3 in units of the displacement, well below
     *         what the near-field refinement resolves.
     *  @param refinement how finely each triangle is subdivided for the quadrature.  The
     *         triangulation joins the nodes, so its triangles are as large as the boundary's
     *         node spacing -- far too coarse to resolve a displacement that varies from node to
     *         node.  Each face is split into \p refinement^2 sub-triangles by cutting every edge
     *         into \p refinement parts, and the displacement is interpolated linearly across the
     *         face and sampled at each sub-triangle.  This refines only the quadrature: the faces
     *         are planar, so subdividing them changes no geometry, and the side test, the surface
     *         proximity test and the exported mesh all continue to use the unsubdivided faces.
     *         \p refinement=1 reproduces the original one-value-per-face rule exactly.
     *
     *         The subdivision is applied to the nearest 3x3 block of images only.  Beyond that
     *         the kernel varies little across a whole face, so the face's mean displacement is
     *         already the right weight for it and subdividing would buy nothing; restricting the
     *         refinement this way keeps the cost near that of the unrefined sum instead of
     *         multiplying the whole image sum by refinement^2.  Neither the total solid angle nor
     *         the rigid-translation case is affected by where the split falls, because the
     *         sub-triangles tile their parent exactly and a uniform displacement interpolates to
     *         the same constant on all of them.
     */
    GbFacet(const std::vector<std::vector<double>>& point_cloud,
            const std::vector<double>& box_dim,
            const Topology& topology,
            const int& normalSense = 1,
            const int& imageShells = 4,
            const int& refinement = 4);

    /*! \brief Displacement at a point \f$\mathbf x\f$ outside the facet.
     *
     *  \f[ \mathbf u(\mathbf x)=\frac{1}{2\pi}\int_{S_\infty}
     *      \mathbf u_S(\mathbf y)\,
     *      \frac{(\mathbf x-\mathbf y)\cdot\hat{\mathbf n}(\mathbf y)}
     *            {|\mathbf x-\mathbf y|^3}\,\mathrm dS(\mathbf y) \f]
     *
     *  i.e. the nodal displacement weighted by the solid angle the facet subtends at
     *  \f$\mathbf x\f$, integrated over \f$S_\infty\f$: the facet *and all of its periodic
     *  images*.  Only the images make the field periodic, and only they make it converge to the
     *  right limit -- a single copy of a finite patch subtends a solid angle that dies off with
     *  distance instead of tending to the \f$2\pi\f$ of a sheet.
     *
     *  The \f$1/2\pi\f$ is fixed by the case the construction has to reproduce: a doubly
     *  periodic sheet subtends exactly \f$\pm2\pi\f$ (independently of how it is faceted), so a
     *  uniform nodal displacement \f$\mathbf b\f$ must give the rigid translation
     *  \f$\pm\mathbf b\f$, which is what a grain translated by \f$\pm\mathbf t/2\f$ does.
     *
     *  At a node of the facet the integral is not used at all: that node's own displacement, which
     *  came from the \f$(\mathbf t,\mathbf s)\f$ pair, is returned instead.  The quadrature cannot
     *  supply it -- the kernel is singular on the surface and the solid angle jumps from
     *  \f$+2\pi\f$ to \f$-2\pi\f$ across it, so the field there is undefined, not merely coarse.
     *
     *  \warning Between the nodes the kernel is still singular on the surface and the quadrature
     *  does not resolve it, so accuracy degrades for a point lying on the facet that is not one of
     *  its nodes.
     */
    Eigen::Vector3d displacement(const Eigen::Vector3d& x) const;

    /*! How close a point must be to a node -- allowing for whole in-plane periods -- before
     *  displacement() hands back that node's own displacement instead of integrating. */
    static constexpr double nodeTolerance = 1.0e-6;

    /*! Solid angle subtended at \p x by the facet together with all of its periodic images:
     *  \f$+2\pi\f$ on the positive side of the facet, \f$-2\pi\f$ on the negative side.  The
     *  magnitude is a property of the sheet, not of its shape, so only the side has to be
     *  determined. */
    double solidAngle(const Eigen::Vector3d& x) const;

    /*! \brief Whether \p x lies on the facet itself, within \p tolerance.
     *
     *  The side test answers +1 or -1 and has no third state, so a point sitting exactly on the
     *  surface gets pushed arbitrarily to one side.  The nodes of the facet are lattice sites of
     *  the grain, so those atoms are precisely the ones the boundary brings into coincidence and
     *  they have to be recognised rather than sorted onto a side.  Periodic images are covered:
     *  the query runs against the same replicated soup the side test uses. */
    bool isOnSurface(const Eigen::Vector3d& x, const double& tolerance = 1.0e-6) const;

    /*! +1 if \p x lies on the positive side of the facet (the side its normals point towards),
     *  -1 otherwise.
     *
     *  Decided by counting how many times a ray leaving \p x along the facet normal crosses the
     *  surface: an even number of crossings puts \p x on the positive side, an odd number on the
     *  negative side.  A nearest-point test is not good enough -- for a point down inside a
     *  corrugation the nearest point can lie on the flank of a neighbouring ridge, whose normal
     *  faces away, and the test then reports the wrong side. */
    int sideOf(const Eigen::Vector3d& x) const;

    /*! Area-weighted mean of the nodal displacement, i.e. the rigid translation the field tends
     *  to far from the facet. */
    const Eigen::Vector3d& meanDisplacement() const;

    double signedDistanceAlongNormal(const Eigen::Vector3d& x,const Eigen::Vector3d& normal) const;

    void export_to_vtp(const std::string& filename) const;

    /*! Mesh vertex positions, one row per mesh vertex.  A node of the point cloud appears
     * once per periodic image referenced by the topology, so this is longer than the cloud. */
    const Eigen::MatrixXd& vertices() const;

    /*! Nodal displacements, row-aligned with vertices(). */
    const Eigen::MatrixXd& displacements() const;

    /*! Triangle connectivity indexing into vertices(), one row per face.  Two facets built
     * from the same Topology have identical face matrices. */
    const Eigen::MatrixXi& faces() const;

    /*! Face normals, one row per face. */
    const Eigen::MatrixXd& faceNormals() const;

    /*! vertices() + displacements(): where this facet's surface actually sits.  Two facets
     * that glue into one continuous surface return the same matrix here. */
    Eigen::MatrixXd deformedVertices() const;

private:

    struct MeshIntegrationData
    {
        Eigen::MatrixXd vertices;
        Eigen::MatrixXd displacements;
        Eigen::MatrixXi faces;
        Eigen::MatrixXd normals;
    };

    struct IntegrationCache
    {
        /*! Quadrature panels: the corners of every sub-triangle, one row each.  A face
         *  contributes refinement^2 of them, laid out consecutively. */
        Eigen::MatrixXd Q0;
        Eigen::MatrixXd Q1;
        Eigen::MatrixXd Q2;

        /*! The displacement carried by each sub-triangle: the linear interpolant of the face's
         *  three nodal values, evaluated at the sub-triangle's centroid.  Because the interpolant
         *  is linear, that centroid value is also its exact mean over the sub-triangle, so this
         *  is a midpoint rule and not merely a sample. */
        Eigen::MatrixXd panelDisp;

        /*! The same interpolation as \p panelDisp, but kept as the three barycentric
         *  coefficients rather than applied.  panelDisp is what the direct sum multiplies by;
         *  these are what the memo stores, since they are the part that does not depend on the
         *  state.  Row \p p sums to one. */
        Eigen::MatrixXd panelBary;

        /*! The unsubdivided faces, used for the distant images and for the mean displacement
         *  and net orientation, which are properties of the surface rather than of the
         *  quadrature. */
        Eigen::MatrixXd F0;
        Eigen::MatrixXd F1;
        Eigen::MatrixXd F2;
        Eigen::MatrixXd faceMeanDisp;

        Eigen::MatrixXd normals;
        Eigen::VectorXd faceAreas;

        /*! Per face: centroid, vector area S = 1/2 (B-A) x (C-A), and the radius of the
         *  smallest sphere about the centroid holding the face.  Expanding the closed form
         *  about the centroid gives N -> 2 (rho . S) and D -> 4 R^3, so
         *  Omega -> -(rho . S)/R^3 -- the linear term vanishing because the centroid is where
         *  it is, which makes the error second order in (radius / distance). */
        Eigen::MatrixXd faceCentroid;
        Eigen::MatrixXd faceVectorArea;
        Eigen::VectorXd faceRadius;

        /*! Each face's identity as a small integer -- see the interning in GbFacet.cpp.
         *  Assigned once when the facet is built, so the nine coordinates are hashed once per
         *  face per mesostate instead of once per face per atom. */
        Eigen::VectorXi faceIdentity;
    };

    /*! The closest-point structure behind sideOf() and signedDistanceAlongNormal().
     *
     * A plain triangle soup, not a Surface_mesh.  The faces of a cut-open periodic surface do not
     * in general glue up into a manifold patch -- each face is represented by whichever periodic
     * copy keeps it compact, and neighbouring faces need not agree on which copy of a shared
     * vertex they use -- so requiring a manifold here would reject perfectly good geometry.  A
     * soup has no such constraint, and closest-point queries do not need connectivity.
     *
     * The soup covers a 3x3 block of period images so that a query point anywhere in or beside
     * the base cell finds its true nearest point rather than one across a seam.
     */
    struct MeshBundle
    {
        using Kernel    = CGAL::Exact_predicates_inexact_constructions_kernel;
        using Point     = Kernel::Point_3;
        using Vector    = Kernel::Vector_3;
        using Triangle  = Kernel::Triangle_3;
        using Iterator  = std::vector<Triangle>::const_iterator;
        using Primitive = CGAL::AABB_triangle_primitive_3<Kernel,Iterator>;
        using Traits    = CGAL::AABB_traits_3<Kernel,Primitive>;
        using Tree      = CGAL::AABB_tree<Traits>;

        std::vector<Triangle> soup;
        /*! Row of \p meshData.faces that each triangle of \p soup came from. */
        std::vector<int> soupToRow;
        std::unique_ptr<Tree> tree;
    };

    MeshBundle bundle;

    MeshIntegrationData meshData;

    /*! This facet's quadrature settings as a small integer, for the memo key.
     *
     *  Declared BEFORE integrationCache deliberately.  Members are initialised in declaration
     *  order, and build_integration_cache() -- which runs while integrationCache is being
     *  initialised -- assigns this; declared after, its own initialiser would run second and
     *  wipe the value.  Mutable because that function is const. */
    mutable int settingsIdentity= 0;

    IntegrationCache integrationCache;

    Eigen::Vector3d period1;
    Eigen::Vector3d period2;

    /*! Area-weighted mean face normal: the facet's overall orientation, already carrying the
     *  normalSense the facet was built with. */
    Eigen::Vector3d planeNormal;

    /*! Centroid of the mesh vertices, used to fold a query point back over the base cell. */
    Eigen::Vector3d centroid;

    /*! How far the surface reaches either side of \p centroid along \p planeNormal: the lowest
     *  and highest signed height of any vertex.
     *
     *  sideOf() decides which side of the sheet a point is on by shooting rays upward along
     *  planeNormal and counting crossings, which is the only way to answer it in general -- but
     *  not for a point outside this band.  The sheet is piecewise linear between its vertices and
     *  its periodic images are translated in plane only, so no part of it lies outside the band;
     *  a point above it cannot be reached by an upward ray, and a point below it is crossed
     *  exactly once.  Both answers follow from the height alone.
     *
     *  This is an exact short circuit, not an approximation, and it takes most of the queries:
     *  the boundary occupies a couple of Angstrom of a crystal tens of Angstrom thick, and every
     *  atom outside that band is classified by one dot product.  Measured at 35.7 -> 24.6
     *  core-seconds of box() over a 1119-state sweep, with every state unchanged. */
    double bandLow= 0.0, bandHigh= 0.0;

    /*! Area-weighted mean nodal displacement over one period cell. */
    Eigen::Vector3d meanDisp;

    /*! The same mean weighted by projected rather than true area -- what the far-field tail
     *  actually sees.  See the closure in displacement(). */
    Eigen::Vector3d meanDispProjected;

    int refinement;

    int imageShells;



    IntegrationCache build_integration_cache() const;

    /*! Translates \p x by whole periods so that it sits over the base cell.  The field is
     *  periodic, so this is exact; it is what makes the truncated image sum periodic too, and it
     *  keeps queries inside the region the triangle soup covers. */
    Eigen::Vector3d foldIntoCell(const Eigen::Vector3d& x) const;

    /*! One pass of the quadrature over the base facet, evaluated at \p y.  Returns both the
     *  displacement-weighted integral and the bare solid angle, so that the far-field closure
     *  can subtract exactly what the explicit shells already contributed.
     *  @param subdivided integrate over the refined panels rather than the whole faces. */
    void accumulate(const Eigen::Vector3d& y,
                    Eigen::Vector3d& weighted,
                    double& omega,
                    const bool& subdivided) const;

    /*! The image sum for one triangle at one query point, from the memo if it is there.
     *  \p x must already be folded into the base cell. */
    TriangleWeights triangleWeights(const Eigen::Vector3d& x, const int& face,
                                    const int& pointIdentity) const;

    /*! The same sum, computed rather than looked up. */
    TriangleWeights computeTriangleWeights(const Eigen::Vector3d& x, const int& face) const;
};

#endif
