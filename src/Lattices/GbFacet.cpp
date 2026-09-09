//
// Created by himanshu on 9/1/26.
//

//GB Facet.cpp
// Created by Himanshu Joshi on 6/5/26.
//

#include "../../include/Lattices/GbFacet.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>

#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>

// Quiet by default; see the declaration for why.
bool GbFacet::announceConstruction = false;
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Periodic_2_triangulation_vertex_base_2.h>
#include <CGAL/Periodic_2_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_data_structure_2.h>

namespace
{
    /*! Orthonormal in-plane frame \f$\{e_1,e_2\}\f$ with \f$e_1\times e_2=\hat n\f$. */
    void planeFrame(const Eigen::Vector3d& normal,
                    Eigen::Vector3d& n_hat, Eigen::Vector3d& e1, Eigen::Vector3d& e2)
    {
        if (normal.norm() < 1e-14)
            throw std::runtime_error("GbFacet: the GB normal is a zero vector.");
        n_hat = normal.normalized();
        const Eigen::Vector3d ref = (std::abs(n_hat(0)) < 0.9) ? Eigen::Vector3d::UnitX()
                                                               : Eigen::Vector3d::UnitY();
        e1 = (ref - ref.dot(n_hat) * n_hat).normalized();
        e2 = n_hat.cross(e1).normalized();
    }

    void periodVectors(const std::vector<double>& box_dim,
                       Eigen::Vector3d& pv1, Eigen::Vector3d& pv2)
    {
        if (box_dim.size() < 6)
            throw std::runtime_error("GbFacet: box_dim must hold the two in-plane period "
                                     "vectors as {p1x,p1y,p1z,p2x,p2y,p2z}.");
        pv1 << box_dim[0], box_dim[1], box_dim[2];
        pv2 << box_dim[3], box_dim[4], box_dim[5];
    }

    /*! Canonical key of an abstract torus face.  Each of the three rotations is written out with
     *  the leading corner's offset subtracted, and the smallest is kept.  Two periodic copies of
     *  the same face then collapse onto the same key, so this deduplicates the 9-sheeted covering
     *  exactly, using integer node indices only -- no floating-point coordinate matching.
     *
     *  All three rotations have to be tried.  Anchoring on the corner with the smallest node index
     *  does not pin the rotation when a face carries the same node at more than one corner, which
     *  is what happens for a cell holding only one or two nodes: there the corners are periodic
     *  images of a single node, the three rotations produce three different keys, and the face gets
     *  counted more than once.
     */
    std::array<long,9> faceKey(const GbFacet::Triangle& t)
    {
        std::array<long,9> best{};
        bool haveBest = false;
        for (int b = 0; b < 3; ++b) {
            const GbFacet::Triangle r{t[b], t[(b+1)%3], t[(b+2)%3]};
            const int dx = r[0].o1;
            const int dy = r[0].o2;
            std::array<long,9> k{};
            for (int i = 0; i < 3; ++i) {
                k[3*i    ] = static_cast<long>(r[i].node);
                k[3*i + 1] = r[i].o1 - dx;
                k[3*i + 2] = r[i].o2 - dy;
            }
            if (!haveBest || k < best) { best = k; haveBest = true; }
        }
        return best;
    }
}

GbFacet::Topology GbFacet::triangulate(const std::vector<Eigen::Vector3d>& points,
                                       const std::vector<double>& box_dim,
                                       const Eigen::Vector3d& normal)
{
    // One node is enough.  A single node tiles the torus with two triangles (V=1, E=3, F=2, so
    // V-E+F=0), and the periodic offsets place its images -- there is no need for three nodes in
    // the cell.  The real invariant is the 2n face count checked at the end.
    if (points.empty())
        throw std::runtime_error("GbFacet::triangulate: need at least one point.");

    Eigen::Vector3d n_hat, e1, e2;
    planeFrame(normal, n_hat, e1, e2);

    Eigen::Vector3d pv1, pv2;
    periodVectors(box_dim, pv1, pv2);

    // Map the period parallelogram to the unit-square torus [0,1)^2 via M = [p1 | p2],
    // expressed in the in-plane frame.
    Eigen::Matrix2d M;
    M.col(0) << pv1.dot(e1), pv1.dot(e2);
    M.col(1) << pv2.dot(e1), pv2.dot(e2);
    if (std::abs(M.determinant()) < 1e-12)
        throw std::runtime_error("GbFacet::triangulate: the two period vectors are parallel or "
                                 "do not span the GB plane.");
    const Eigen::Matrix2d M_inv = M.inverse();

    using K2  = CGAL::Exact_predicates_inexact_constructions_kernel;
    using GT  = CGAL::Periodic_2_Delaunay_triangulation_traits_2<K2>;
    using PVb = CGAL::Periodic_2_triangulation_vertex_base_2<GT>;
    using Vb  = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, GT, PVb>;
    using Fb  = CGAL::Periodic_2_triangulation_face_base_2<GT>;
    using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
    using PDT = CGAL::Periodic_2_Delaunay_triangulation_2<GT, TDS>;

    const std::size_t n = points.size();
    std::vector<Eigen::Vector2d> uv(n);          // canonical torus coordinates, in [0,1)^2
    std::vector<std::array<int,2>> wrap(n);      // periods added to reach that canonical position
    std::vector<std::pair<PDT::Point, std::size_t>> pts2d;
    pts2d.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const Eigen::Vector2d q(points[i].dot(e1), points[i].dot(e2));
        const Eigen::Vector2d c = M_inv * q;
        // An input node need not lie in the base cell.  Record how many periods were added to
        // bring it there, because the offsets CGAL reports are relative to the canonical position,
        // not to the position the node was handed in at.
        double u0 = c(0) - std::floor(c(0));
        double u1 = c(1) - std::floor(c(1));
        int w0 = static_cast<int>(-std::floor(c(0)));
        int w1 = static_cast<int>(-std::floor(c(1)));
        // The domain is half open, and CGAL requires the point to be inside it.  For a coordinate
        // a hair below zero, floor() gives -1 and the subtraction rounds to exactly 1, landing the
        // point on the excluded face; fold it back to the other side and take the period with it.
        // Left unhandled this rejects a fifth of all mesostates -- and, with CGAL's own checks
        // compiled out, corrupted the triangulation instead of rejecting anything.
        if (!(u0 < 1.0)) { u0 = 0.0; --w0; }
        if (!(u1 < 1.0)) { u1 = 0.0; --w1; }
        if (u0 < 0.0) u0 = 0.0;
        if (u1 < 0.0) u1 = 0.0;
        uv[i] << u0, u1;
        wrap[i] = {w0, w1};
        pts2d.push_back({PDT::Point(uv[i](0), uv[i](1)), i});
    }

    PDT pdt(GT::Iso_rectangle_2(0, 0, 1, 1));
    pdt.insert(pts2d.begin(), pts2d.end());

    // A node silently dropped here would be a node missing from the surface, so insist that
    // every one survived.  Nodes collide only if two of them project onto the same in-plane
    // point modulo the periods, which means the mesostate itself is degenerate.
    if (pdt.number_of_vertices() != n)
        throw std::runtime_error("GbFacet::triangulate: " + std::to_string(n) + " nodes collapsed to "
                                 + std::to_string(pdt.number_of_vertices()) + " after projection onto "
                                 "the GB plane -- two nodes share an in-plane position modulo the periods.");

    // CGAL keeps the triangulation in a 9-sheeted covering, so faces_begin() walks 9 copies of
    // every abstract torus face.  Take the vertex indices and offsets straight off each face and
    // deduplicate by canonical key.
    Topology topology;
    std::set<std::array<long,9>> seen;
    for (auto f = pdt.faces_begin(); f != pdt.faces_end(); ++f) {
        Triangle t;
        for (int i = 0; i < 3; ++i) {
            const auto off = pdt.get_offset(f, i);
            const std::size_t id = f->vertex(i)->info();
            // get_offset() locates the corner relative to the node's canonical position, so add
            // back the wrap to express it relative to the position in `points`.  Skipping this
            // puts any node that was handed in from outside the base cell one period out of place,
            // which lifts its faces to cell-spanning, orientation-reversed triangles.
            t[i] = Corner{id, off.x() + wrap[id][0], off.y() + wrap[id][1]};
        }

        if (!seen.insert(faceKey(t)).second) continue;
        topology.push_back(t);
    }

    // Orient the faces consistently.  This cannot be done face by face from the geometry: a
    // periodic Delaunay triangulation of few points emits triangles whose three nodes are exactly
    // collinear in projection, and a collapsed triangle has no orientation of its own.  Consistent
    // orientation is a combinatorial property of the surface, so propagate it across shared edges
    // instead, and only then use the geometry -- once, on a face that does have area -- to decide
    // which way round the whole surface should face.
    {
        // Translation-invariant descriptor of a directed edge: the two node indices plus the
        // offset difference between its endpoints.  Two faces meet along a torus edge when their
        // descriptors agree; if they agree in the *same* direction, one of them is flipped.
        using EdgeKey = std::array<long,4>;
        auto edgeKey = [](const Corner& a, const Corner& b, bool& forward) -> EdgeKey {
            const long d1 = b.o1 - a.o1;
            const long d2 = b.o2 - a.o2;
            const bool asIs = (a.node < b.node)
                           || (a.node == b.node && (d1 > 0 || (d1 == 0 && d2 > 0)));
            forward = asIs;
            if (asIs) return EdgeKey{(long)a.node,(long)b.node, d1, d2};
            return EdgeKey{(long)b.node,(long)a.node,-d1,-d2};
        };

        // edge -> the faces on it, with the direction each one traverses it in
        std::map<EdgeKey,std::vector<std::pair<std::size_t,bool>>> edges;
        for (std::size_t f = 0; f < topology.size(); ++f)
            for (int i = 0; i < 3; ++i) {
                bool forward = false;
                const EdgeKey k = edgeKey(topology[f][i], topology[f][(i+1)%3], forward);
                edges[k].push_back({f, forward});
            }
        for (const auto& [k, incident] : edges)
            if (incident.size() != 2)
                throw std::runtime_error("GbFacet::triangulate: an edge is shared by "
                                         + std::to_string(incident.size())
                                         + " faces; the triangulation is not a surface.");

        std::vector<char> visited(topology.size(), 0);
        std::vector<std::size_t> stack{0};
        visited[0] = 1;
        std::size_t reached = 1;
        while (!stack.empty()) {
            const std::size_t f = stack.back();
            stack.pop_back();
            for (int i = 0; i < 3; ++i) {
                bool forwardF = false;
                const EdgeKey k = edgeKey(topology[f][i], topology[f][(i+1)%3], forwardF);
                for (const auto& [g, forwardG] : edges.at(k)) {
                    if (g == f || visited[g]) continue;
                    // Neighbours must traverse their shared edge in opposite directions.
                    if (forwardG == forwardF) std::swap(topology[g][1], topology[g][2]);
                    visited[g] = 1;
                    ++reached;
                    stack.push_back(g);
                }
            }
        }
        if (reached != topology.size())
            throw std::runtime_error("GbFacet::triangulate: the triangulation is disconnected ("
                                     + std::to_string(reached) + " of "
                                     + std::to_string(topology.size()) + " faces reached).");

        // Global sense: a Delaunay triangle projects to a positive-area triangle when it faces
        // along n_hat, and e1 x e2 = n_hat.  Pick the largest face so the test is well conditioned.
        double best = 0.0;
        double bestSigned = 0.0;
        for (const auto& t : topology) {
            Eigen::Vector2d q[3];
            for (int i = 0; i < 3; ++i) {
                const Eigen::Vector3d p = points[t[i].node] + t[i].o1 * pv1 + t[i].o2 * pv2;
                q[i] << p.dot(e1), p.dot(e2);
            }
            const double a = 0.5 * ((q[1]-q[0])(0)*(q[2]-q[0])(1) - (q[1]-q[0])(1)*(q[2]-q[0])(0));
            if (std::abs(a) > best) { best = std::abs(a); bestSigned = a; }
        }
        if (best < 1e-14)
            throw std::runtime_error("GbFacet::triangulate: every triangle is collapsed in "
                                     "projection, so the surface has no orientation.");
        if (bestSigned < 0.0)
            for (auto& t : topology) std::swap(t[1], t[2]);
    }

    // Re-anchor every face so that its centroid falls in the base period cell.  The dedupe key is
    // translation invariant, so it is free to pick any image as the representative, and anchoring
    // on the lowest node index picks images that do not tile the cell -- the representatives then
    // overlap in places and leave gaps in others, which makes the ray crossing count in sideOf()
    // come out wrong.  Anchoring on the centroid makes them a genuine fundamental domain.
    {
        Eigen::Matrix<double,3,2> Pp;
        Pp.col(0) = pv1;
        Pp.col(1) = pv2;
        for (auto& t : topology) {
            Eigen::Vector3d C = Eigen::Vector3d::Zero();
            for (int i = 0; i < 3; ++i)
                C += points[t[i].node] + t[i].o1 * pv1 + t[i].o2 * pv2;
            C /= 3.0;
            const Eigen::Vector2d cc = Pp.colPivHouseholderQr().solve(C);
            const int s1 = static_cast<int>(std::floor(cc(0)));
            const int s2 = static_cast<int>(std::floor(cc(1)));
            for (int i = 0; i < 3; ++i) { t[i].o1 -= s1; t[i].o2 -= s2; }
        }
    }

    // A triangulation of n points on a torus has exactly 2n faces (V - E + F = 0, E = 3F/2).
    if (topology.size() != 2 * n)
        throw std::runtime_error("GbFacet::triangulate: got " + std::to_string(topology.size())
                                 + " faces for " + std::to_string(n) + " nodes, expected "
                                 + std::to_string(2 * n) + ".");

    return topology;
}

// Constructor
GbFacet::GbFacet(const std::vector<std::vector<double>>& point_cloud,
                 const std::vector<double>& box_dim,
                 const Topology& topology,
                 const int& normalSense,
                 const int& imageShells,
                 const int& refinement)
    : refinement(refinement), imageShells(imageShells)
{
    if (point_cloud.empty())
        throw std::runtime_error("GbFacet: need at least one point.");
    if (topology.empty())
        throw std::runtime_error("GbFacet: empty topology.");
    if (normalSense != 1 && normalSense != -1)
        throw std::runtime_error("GbFacet: normalSense must be +1 or -1.");
    if (refinement < 1)
        throw std::runtime_error("GbFacet: refinement must be >= 1.");
    if (imageShells < 0)
        throw std::runtime_error("GbFacet: imageShells must be >= 0.");

    Eigen::Vector3d pv1, pv2;
    periodVectors(box_dim, pv1, pv2);
    period1 = pv1;
    period2 = pv2;

    // Enumerate mesh vertices: one per (node, o1, o2) corner the topology actually references.
    // The iteration order depends only on the topology, so two facets sharing a topology get
    // identical vertex numbering and an identical face matrix.
    std::map<std::array<long,3>, int> cornerToVertex;
    std::vector<Eigen::Vector3d> pos, disp;
    pos.reserve(topology.size() * 3);
    disp.reserve(topology.size() * 3);

    auto vertexOf = [&](const Corner& c) -> int {
        const std::array<long,3> key{static_cast<long>(c.node), c.o1, c.o2};
        const auto it = cornerToVertex.find(key);
        if (it != cornerToVertex.end()) return it->second;

        if (c.node >= point_cloud.size())
            throw std::runtime_error("GbFacet: the topology references node "
                                     + std::to_string(c.node) + " but the point cloud holds only "
                                     + std::to_string(point_cloud.size()) + " nodes.");
        if (point_cloud[c.node].size() < 6)
            throw std::runtime_error("GbFacet: point cloud row " + std::to_string(c.node)
                                     + " must hold {x,y,z,ux,uy,uz}.");

        const auto& r = point_cloud[c.node];
        const int v = static_cast<int>(pos.size());
        // The displacement field is periodic, so every periodic image of a node carries the
        // node's displacement: image x_i + o1*p1 + o2*p2 deforms to (x_i+u_i) + o1*p1 + o2*p2.
        pos.push_back(Eigen::Vector3d(r[0], r[1], r[2]) + c.o1 * pv1 + c.o2 * pv2);
        disp.push_back(Eigen::Vector3d(r[3], r[4], r[5]));
        cornerToVertex[key] = v;
        return v;
    };

    meshData.faces.resize(static_cast<int>(topology.size()), 3);
    for (std::size_t f = 0; f < topology.size(); ++f)
        for (int j = 0; j < 3; ++j)
            meshData.faces(static_cast<int>(f), j) = vertexOf(topology[f][j]);

    const int nVertices = static_cast<int>(pos.size());
    meshData.vertices.resize(nVertices, 3);
    meshData.displacements.resize(nVertices, 3);
    for (int i = 0; i < nVertices; ++i) {
        meshData.vertices.row(i)      = pos[i].transpose();
        meshData.displacements.row(i) = disp[i].transpose();
    }

    // Face normals straight off the (already consistently oriented) triangles.
    //
    // A face can be collapsed here even though it is well shaped in the configuration the
    // topology was triangulated in: the reference facet is a stepped surface, and three nodes
    // that span a proper triangle once deformed may be exactly collinear beforehand.  Such a
    // face has no reference area, so it must contribute nothing rather than abort the build --
    // its normal is left at zero, which zeroes its weight in every surface integral.  The face
    // is kept in the face matrix so that both facets keep identical connectivity.
    const int nFaces = static_cast<int>(meshData.faces.rows());
    meshData.normals.resize(nFaces, 3);
    std::vector<bool> degenerate(nFaces, false);
    int nDegenerate = 0;
    for (int f = 0; f < nFaces; ++f) {
        const Eigen::Vector3d p0 = meshData.vertices.row(meshData.faces(f,0));
        const Eigen::Vector3d p1 = meshData.vertices.row(meshData.faces(f,1));
        const Eigen::Vector3d p2 = meshData.vertices.row(meshData.faces(f,2));
        const Eigen::Vector3d nrm = (p1 - p0).cross(p2 - p0);
        const double len = nrm.norm();
        if (len < 1e-14) {
            meshData.normals.row(f).setZero();
            degenerate[f] = true;
            ++nDegenerate;
        }
        else
            meshData.normals.row(f) = (normalSense * nrm / len).transpose();
    }

    // Closest-point structure: a triangle soup over a 3x3 block of period images (see MeshBundle).
    bundle.soup.clear();
    bundle.soupToRow.clear();
    bundle.soup.reserve(static_cast<std::size_t>(nFaces) * 9);
    bundle.soupToRow.reserve(static_cast<std::size_t>(nFaces) * 9);
    for (int m = -1; m <= 1; ++m)
        for (int k = -1; k <= 1; ++k) {
            const Eigen::Vector3d T = m * pv1 + k * pv2;
            for (int f = 0; f < nFaces; ++f) {
                if (degenerate[f]) continue;   // a zero-area triangle has no well defined closest point
                const Eigen::Vector3d p0 = meshData.vertices.row(meshData.faces(f,0)) + T.transpose();
                const Eigen::Vector3d p1 = meshData.vertices.row(meshData.faces(f,1)) + T.transpose();
                const Eigen::Vector3d p2 = meshData.vertices.row(meshData.faces(f,2)) + T.transpose();
                bundle.soup.emplace_back(MeshBundle::Point(p0(0),p0(1),p0(2)),
                                         MeshBundle::Point(p1(0),p1(1),p1(2)),
                                         MeshBundle::Point(p2(0),p2(1),p2(2)));
                bundle.soupToRow.push_back(f);
            }
        }
    if (bundle.soup.empty())
        throw std::runtime_error("GbFacet: every face is degenerate.");

    bundle.tree = std::make_unique<MeshBundle::Tree>(bundle.soup.cbegin(), bundle.soup.cend());
    bundle.tree->accelerate_distance_queries();

    integrationCache = build_integration_cache();

    // Area-weighted mean nodal displacement over one period cell.  Far from the facet the solid
    // angle spreads uniformly over the images, so this is the rigid translation the field tends
    // to -- and it is what the far-field closure in displacement() uses.
    const double totalArea = integrationCache.faceAreas.sum();
    if (totalArea < 1e-14)
        throw std::runtime_error("GbFacet: the facet has zero total area.");
    meanDisp = (integrationCache.faceMeanDisp.array().colwise()
                    * integrationCache.faceAreas.array())
                   .colwise().sum().transpose() / totalArea;

    // Overall orientation of the facet, for the ray direction in sideOf().  Taken from the faces
    // themselves so that it automatically carries whichever sense the facet was built with.
    const Eigen::Vector3d meanNormal =
        (meshData.normals.array().colwise() * integrationCache.faceAreas.array())
            .colwise().sum().transpose() / totalArea;
    if (meanNormal.norm() < 1e-6)
        throw std::runtime_error("GbFacet: the facet has no net orientation -- its face normals "
                                 "cancel, so it cannot be treated as a sheet.");
    planeNormal = meanNormal.normalized();
    centroid    = meshData.vertices.colwise().mean().transpose();

    if (!announceConstruction) return;
    std::cout << "Mesh built with " << nVertices << " vertices (" << point_cloud.size()
              << " nodes + periodic copies), " << nFaces << " faces";
    if (nDegenerate > 0)
        std::cout << " (" << nDegenerate << " with zero area in this configuration)";
    if (refinement > 1)
        std::cout << ", integrated over " << integrationCache.Q0.rows()
                  << " panels (" << refinement << "x refinement)";
    std::cout << "." << std::endl;
}

const Eigen::MatrixXd& GbFacet::vertices() const      { return meshData.vertices; }
const Eigen::MatrixXd& GbFacet::displacements() const { return meshData.displacements; }
const Eigen::MatrixXi& GbFacet::faces() const         { return meshData.faces; }
const Eigen::MatrixXd& GbFacet::faceNormals() const   { return meshData.normals; }

const Eigen::Vector3d& GbFacet::meanDisplacement() const { return meanDisp; }

Eigen::MatrixXd GbFacet::deformedVertices() const
{
    return meshData.vertices + meshData.displacements;
}

GbFacet::IntegrationCache GbFacet::build_integration_cache() const
{
    IntegrationCache cache;
    const auto& data = meshData;

    cache.F0 = data.vertices(data.faces.col(0), Eigen::placeholders::all);
    cache.F1 = data.vertices(data.faces.col(1), Eigen::placeholders::all);
    cache.F2 = data.vertices(data.faces.col(2), Eigen::placeholders::all);
    const Eigen::MatrixXd& P0 = cache.F0;
    const Eigen::MatrixXd& P1 = cache.F1;
    const Eigen::MatrixXd& P2 = cache.F2;

    const Eigen::MatrixXd D0 = data.displacements(data.faces.col(0), Eigen::placeholders::all);
    const Eigen::MatrixXd D1 = data.displacements(data.faces.col(1), Eigen::placeholders::all);
    const Eigen::MatrixXd D2 = data.displacements(data.faces.col(2), Eigen::placeholders::all);

    cache.normals = data.normals;

    const Eigen::MatrixXd E1 = P1 - P0;
    const Eigen::MatrixXd E2 = P2 - P0;

    Eigen::MatrixXd Cross(E1.rows(),3);
    Cross.col(0) = E1.col(1).cwiseProduct(E2.col(2)) - E1.col(2).cwiseProduct(E2.col(1));
    Cross.col(1) = E1.col(2).cwiseProduct(E2.col(0)) - E1.col(0).cwiseProduct(E2.col(2));
    Cross.col(2) = E1.col(0).cwiseProduct(E2.col(1)) - E1.col(1).cwiseProduct(E2.col(0));

    cache.faceAreas   = 0.5 * Cross.rowwise().norm();
    cache.faceMeanDisp = (D0 + D1 + D2) / 3.0;

    // Subdivide every face for the quadrature.  The triangulation only joins the nodes, so its
    // triangles span the node spacing of the boundary; taking one displacement value per triangle
    // resolves nothing of how the displacement varies between neighbouring nodes.  Cutting each
    // edge into `refinement` parts splits a face into refinement^2 similar sub-triangles, and the
    // displacement is interpolated linearly across the face and sampled on each of them.
    //
    // Two properties of the parent-level scheme survive this, and they are the reason it is done
    // by subdivision rather than by a higher-order panel formula:
    //
    //   * the sub-triangles tile the parent exactly, so their solid angles sum to the parent's --
    //     the total solid angle, and with it the far-field closure in displacement(), is
    //     unchanged at every refinement;
    //   * a uniform nodal displacement interpolates to the same constant on every sub-triangle,
    //     so the rigid-translation case stays exact for any refinement.
    //
    // The interpolant is linear, so its value at a sub-triangle's centroid is also its exact mean
    // over that sub-triangle: this is a midpoint rule, second order in the panel size, not a
    // point sample.
    const int nFaces   = static_cast<int>(P0.rows());
    const int perFace  = refinement*refinement;
    const int nPanels  = nFaces*perFace;

    cache.Q0.resize(nPanels,3);
    cache.Q1.resize(nPanels,3);
    cache.Q2.resize(nPanels,3);
    cache.panelDisp.resize(nPanels,3);

    // Barycentric lattice point (i,j) of a face, and the displacement interpolated there.
    const double k = static_cast<double>(refinement);
    auto latticePoint = [&](const Eigen::RowVector3d& a, const Eigen::RowVector3d& b,
                            const Eigen::RowVector3d& c, int i, int j)
    {
        return Eigen::RowVector3d(a + (i/k)*(b-a) + (j/k)*(c-a));
    };

    for (int f = 0; f < nFaces; ++f)
    {
        const Eigen::RowVector3d p0 = P0.row(f), p1 = P1.row(f), p2 = P2.row(f);
        const Eigen::RowVector3d d0 = D0.row(f), d1 = D1.row(f), d2 = D2.row(f);
        int panel = f*perFace;

        auto emit = [&](int i0,int j0, int i1,int j1, int i2,int j2)
        {
            cache.Q0.row(panel) = latticePoint(p0,p1,p2,i0,j0);
            cache.Q1.row(panel) = latticePoint(p0,p1,p2,i1,j1);
            cache.Q2.row(panel) = latticePoint(p0,p1,p2,i2,j2);
            // Mean of the interpolant at the three corners == its value at the centroid,
            // the interpolant being linear.
            cache.panelDisp.row(panel) = ( latticePoint(d0,d1,d2,i0,j0)
                                         + latticePoint(d0,d1,d2,i1,j1)
                                         + latticePoint(d0,d1,d2,i2,j2) ) / 3.0;
            ++panel;
        };

        // Upward sub-triangles, one per lattice cell, and the downward ones that fill the gaps
        // between them: refinement*(refinement+1)/2 and refinement*(refinement-1)/2, which is
        // refinement^2 together.
        for (int i = 0; i + 0 <= refinement-1; ++i)
            for (int j = 0; i + j <= refinement-1; ++j)
                emit(i,j, i+1,j, i,j+1);
        for (int i = 0; i + 0 <= refinement-2; ++i)
            for (int j = 0; i + j <= refinement-2; ++j)
                emit(i+1,j, i+1,j+1, i,j+1);

        if (panel != (f+1)*perFace)
            throw std::runtime_error("GbFacet: subdivision produced "
                                     + std::to_string(panel - f*perFace) + " panels for a face, "
                                       "expected " + std::to_string(perFace) + ".");
    }

    return cache;
}

// Public functions
void GbFacet::accumulate(const Eigen::Vector3d& y,
                         Eigen::Vector3d& weighted,
                         double& omega,
                         const bool& subdivided) const
{
    // Each panel's contribution is its solid angle at y, evaluated in closed form by the
    // van Oosterom-Strackee construction:
    //
    //     tan(Omega/2) = a.(b x c) / ( |a||b||c| + (a.b)|c| + (a.c)|b| + (b.c)|a| )
    //
    // with a, b, c the corners measured from y.  The two-argument atan2 is required, not atan:
    // the denominator turns negative once the triangle subtends more than a hemisphere, and the
    // signed numerator carries the orientation, so the result is a signed solid angle.
    //
    // This replaces a 7-point quadrature of the same integral.  The kernel goes as 1/|y-Q|^3 and
    // no fixed rule resolves it close to the surface -- the quadrature lost 57% of the solid angle
    // at a tenth of a triangle's size and 98% at a 250th, which is what threw atoms hundreds of
    // Angstrom out of the box.  The closed form is exact at every distance, and cheaper.
    //
    // The panels are the subdivided faces built in build_integration_cache(), each carrying the
    // linear interpolant of its face's nodal displacements sampled at its own centroid.  The
    // displacement is still piecewise constant, but now over panels far smaller than the node
    // spacing rather than over whole node-to-node triangles, and it is bounded whatever the
    // refinement because |Omega| <= 2*pi.
    const auto& cache = integrationCache;
    const Eigen::MatrixXd& A0 = subdivided ? cache.Q0 : cache.F0;
    const Eigen::MatrixXd& A1 = subdivided ? cache.Q1 : cache.F1;
    const Eigen::MatrixXd& A2 = subdivided ? cache.Q2 : cache.F2;
    const Eigen::MatrixXd& AD = subdivided ? cache.panelDisp : cache.faceMeanDisp;
    const int numPanels = static_cast<int>(A0.rows());

    weighted.setZero();
    omega = 0.0;

    for (int f = 0; f < numPanels; ++f)
    {
        const Eigen::Vector3d a = A0.row(f).transpose() - y;
        const Eigen::Vector3d b = A1.row(f).transpose() - y;
        const Eigen::Vector3d c = A2.row(f).transpose() - y;
        const double aNorm = a.norm();
        const double bNorm = b.norm();
        const double cNorm = c.norm();

        const double numerator   = a.dot(b.cross(c));
        const double denominator = aNorm*bNorm*cNorm
                                 + a.dot(b)*cNorm + a.dot(c)*bNorm + b.dot(c)*aNorm;

        // Negated to match the sign this code integrates with, h = (y - P0).n rather than the
        // (P - y).n of the standard solid angle.  A collapsed panel gives a zero numerator and so
        // contributes nothing, which is what a panel with no area should do.
        const double solidAngleOfPanel = -2.0*std::atan2(numerator, denominator);

        omega    += solidAngleOfPanel;
        weighted += solidAngleOfPanel * AD.row(f).transpose();
    }
}

bool GbFacet::isOnSurface(const Eigen::Vector3d& x, const double& tolerance) const
{
    const Eigen::Vector3d folded = foldIntoCell(x);
    const double squaredDistance = CGAL::to_double(
        bundle.tree->squared_distance(MeshBundle::Point(folded(0), folded(1), folded(2))));
    return squaredDistance <= tolerance*tolerance;
}

Eigen::Vector3d GbFacet::foldIntoCell(const Eigen::Vector3d& x) const
{
    Eigen::Matrix<double,3,2> P;
    P.col(0) = period1;
    P.col(1) = period2;
    const Eigen::Vector3d d = x - centroid;
    const Eigen::Vector2d c = P.colPivHouseholderQr().solve(d);
    return x - P * Eigen::Vector2d(std::round(c(0)), std::round(c(1)));
}

int GbFacet::sideOf(const Eigen::Vector3d& x) const
{
    // The soup only covers a 3x3 block of images, so bring the query point over the base cell.
    const Eigen::Vector3d folded = foldIntoCell(x);

    // A ray that grazes an edge shared by two faces is counted by both of them, so its parity
    // comes out even whichever side the point is on.  That is not a rare accident: a cell holding
    // few nodes is tiled by a couple of large triangles whose shared edge runs right across it.
    // Shoot three rays and take the majority -- a degenerate hit cannot catch all three at once.
    const Eigen::Vector3d reference = (std::abs(planeNormal(0)) < 0.9) ? Eigen::Vector3d::UnitX()
                                                                      : Eigen::Vector3d::UnitY();
    const Eigen::Vector3d inPlane1 = (reference - reference.dot(planeNormal)*planeNormal).normalized();
    const Eigen::Vector3d inPlane2 = planeNormal.cross(inPlane1).normalized();

    // The three rays are separated by moving their ORIGINS in plane, not only by tilting their
    // directions.  A worse degeneracy than a grazed edge is a ray that leaves through a vertex:
    // the nodes are lattice sites, so an atom of the same lattice column sits directly above a
    // node and projects onto it exactly, and the vertex it then passes through is shared by every
    // triangle meeting there.  Tilting the direction by an angle displaces the ray by that angle
    // times the distance travelled -- microscopic next to the node spacing, so all three rays
    // still hit the vertex and the majority confirms the wrong answer rather than escaping it.
    // Offsetting the origin by a fixed fraction of the period clears the vertex outright, while
    // staying orders of magnitude below the node spacing, so it cannot move the point across the
    // surface.  The three offsets are 120 degrees apart so that no single feature can catch them
    // all; the small tilts are kept, since they separate the rays further along their path.
    const double offset = 1.0e-3*std::min(period1.norm(), period2.norm());
    constexpr double tilt = 1.0e-4;
    const Eigen::Vector3d origins[3] = {
        folded + offset*inPlane1,
        folded - 0.5*offset*inPlane1 + 0.866*offset*inPlane2,
        folded - 0.5*offset*inPlane1 - 0.866*offset*inPlane2 };
    const Eigen::Vector3d directions[3] = {
        planeNormal + tilt*inPlane1,
        planeNormal - tilt*inPlane1 + 0.5*tilt*inPlane2,
        planeNormal + 0.7*tilt*inPlane2 };

    int votes = 0;
    for (int r = 0; r < 3; ++r) {
        const MeshBundle::Kernel::Ray_3 ray(
            MeshBundle::Point(origins[r](0), origins[r](1), origins[r](2)),
            MeshBundle::Vector(directions[r](0), directions[r](1), directions[r](2)));
        votes += (bundle.tree->number_of_intersected_primitives(ray) % 2 == 0) ? 1 : -1;
    }
    return (votes >= 0) ? 1 : -1;
}

double GbFacet::solidAngle(const Eigen::Vector3d& x) const
{
    // The solid angle of a doubly periodic sheet is +-2*pi whatever its shape, so only the side
    // needs deciding.  Summing the truncated image series to decide it is not safe: for a point
    // far from the facet, or down inside a corrugation, the explicit shells can carry the wrong
    // sign entirely.
    return sideOf(x) * 2.0 * M_PI;
}

Eigen::Vector3d GbFacet::displacement(const Eigen::Vector3d& x0) const
{
    // A point sitting on a node of the facet is one of the sites the boundary was built from, and
    // its displacement is the one the (t,s) pair prescribes.  The quadrature cannot supply it: the
    // kernel is singular on the surface, and the solid angle jumps from +2*pi to -2*pi across it,
    // so the field is genuinely undefined there rather than merely inaccurate.  Away from the
    // nodes the quadrature below takes over.
    {
        Eigen::Matrix<double,3,2> periods;
        periods.col(0) = period1;
        periods.col(1) = period2;
        for (int v = 0; v < meshData.vertices.rows(); ++v) {
            const Eigen::Vector3d offset =
                x0 - Eigen::Vector3d(meshData.vertices.row(v));
            // the node repeats on the in-plane lattice, so match modulo whole periods
            const Eigen::Vector2d coefficients = periods.colPivHouseholderQr().solve(offset);
            const Eigen::Vector2d whole(std::round(coefficients(0)), std::round(coefficients(1)));
            if ((offset - periods*whole).norm() < nodeTolerance)
                return meshData.displacements.row(v);
        }
    }

    // Centre the window of images on the query point, not on the base cell: a window fixed in
    // space includes a different set of images depending on where x0 sits, which would leave the
    // truncated field only approximately periodic.  Folding first makes it periodic exactly.
    const Eigen::Vector3d x = foldIntoCell(x0);

    // Translating an image of the facet by T is the same as translating the field point by -T,
    // so the images are summed by re-evaluating the base facet at shifted points.
    Eigen::Vector3d weighted = Eigen::Vector3d::Zero();
    double omega = 0.0;
    // Only the nearest block of images is integrated over the refined panels.  Further out the
    // kernel barely varies across a whole face, so the face mean is already the right weight and
    // the subdivision would only cost time; the sub-triangles tile their parent exactly, so the
    // solid angle -- and with it the closure below -- is the same either way.
    for (int m = -imageShells; m <= imageShells; ++m)
        for (int n = -imageShells; n <= imageShells; ++n) {
            Eigen::Vector3d w;
            double o;
            const bool near = (std::abs(m) <= 1 && std::abs(n) <= 1);
            accumulate(x - (m*period1 + n*period2), w, o, near);
            weighted += w;
            omega    += o;
        }

    // Far-field closure.  The explicit shells miss the rest of the infinite sheet, and the tail
    // converges only as 1/imageShells, so truncating it outright is not good enough.  Beyond the
    // explicit shells the kernel varies slowly across a period cell, so those images see only the
    // mean nodal displacement; and the total solid angle of the whole sheet is known exactly to be
    // +-2*pi whatever the faceting.  So replace the tail by meanDisp times the solid angle still
    // unaccounted for.  For a uniform nodal displacement this is exact for any imageShells, and in
    // general only the *variation* of the displacement is left to the explicit shells.
    const double omegaSheet = solidAngle(x);

    return (weighted + meanDisp * (omegaSheet - omega)) / (2.0 * M_PI);
}

double GbFacet::signedDistanceAlongNormal(const Eigen::Vector3d& x,
                                           const Eigen::Vector3d& n) const {
    MeshBundle::Point query(x(0), x(1), x(2));
    auto closest = bundle.tree->closest_point(query);
    Eigen::Vector3d diff(x(0) - CGAL::to_double(closest.x()),
                         x(1) - CGAL::to_double(closest.y()),
                         x(2) - CGAL::to_double(closest.z()));
    return diff.dot(n);
}

void GbFacet::export_to_vtp(const std::string& output_filename) const {
    const auto& data = meshData;
    const int nVertices = static_cast<int>(data.vertices.rows());
    const int nFaces    = static_cast<int>(data.faces.rows());

    std::ofstream os(output_filename);
    if (!os) {
        std::cerr << "Error: Could not open file for VTP export: " << output_filename << std::endl;
        return;
    }
    // Round-trip the doubles exactly: at the default 6 significant digits the exported nodes of
    // the two facets no longer deform onto each other, so the glued surface looks discontinuous
    // to whatever reads the file.
    os << std::setprecision(17);

    os << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
       << "  <PolyData>\n"
       << "    <Piece NumberOfPoints=\"" << nVertices << "\" NumberOfPolys=\"" << nFaces << "\">\n";

    os << "      <Points>\n"
       << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < nVertices; ++i)
        os << "          " << data.vertices(i,0) << " " << data.vertices(i,1) << " " << data.vertices(i,2) << "\n";
    os << "        </DataArray>\n"
       << "      </Points>\n";

    os << "      <PointData Vectors=\"displacement\">\n"
       << "        <DataArray type=\"Float64\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < nVertices; ++i)
        os << "          " << data.displacements(i,0) << " " << data.displacements(i,1) << " " << data.displacements(i,2) << "\n";
    os << "        </DataArray>\n"
       << "      </PointData>\n";

    os << "      <Polys>\n"
       << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 0; i < nFaces; ++i)
        os << "          " << data.faces(i,0) << " " << data.faces(i,1) << " " << data.faces(i,2) << "\n";
    os << "        </DataArray>\n"
       << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 1; i <= nFaces; ++i)
        os << "          " << 3*i << "\n";
    os << "        </DataArray>\n"
       << "      </Polys>\n"
       << "    </Piece>\n"
       << "  </PolyData>\n"
       << "</VTKFile>\n";

    std::cout << "Mesh exported to " << output_filename
              << " (" << nVertices << " vertices, " << nFaces << " faces)" << std::endl;
}

