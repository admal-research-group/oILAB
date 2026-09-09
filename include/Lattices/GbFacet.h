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

    IntegrationCache integrationCache;

    Eigen::Vector3d period1;
    Eigen::Vector3d period2;

    /*! Area-weighted mean face normal: the facet's overall orientation, already carrying the
     *  normalSense the facet was built with. */
    Eigen::Vector3d planeNormal;

    /*! Centroid of the mesh vertices, used to fold a query point back over the base cell. */
    Eigen::Vector3d centroid;

    /*! Area-weighted mean nodal displacement over one period cell. */
    Eigen::Vector3d meanDisp;

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
    TriangleWeights triangleWeights(const Eigen::Vector3d& x, const int& face) const;

    /*! The same sum, computed rather than looked up. */
    TriangleWeights computeTriangleWeights(const Eigen::Vector3d& x, const int& face) const;
};

#endif
