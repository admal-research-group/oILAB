//
// Analytic checks on GbFacet: the periodic triangulation and the solid-angle displacement field.
//
// The load-bearing check is the uniform-jump one.  A doubly periodic sheet subtends exactly
// +-2*pi at any point off it, whatever its shape, so a uniform nodal displacement b must give
// the rigid translation +-b -- exactly, with no dependence on how many image shells are summed.
// That is the case the mesostate construction has to reproduce (a grain translated by +-t/2), and
// it is sensitive to almost everything that can go wrong: the 1/2*pi normalisation, the sign of
// the normals, the periodic offsets, the far-field closure and the side test.
//
// It does have one blind spot: an error that enters the displacement-weighted integral and the
// bare solid angle the same way cancels in the ratio.  A face lifted to the wrong periodic image
// is exactly such an error, so the tiling check below tests the geometry directly.
//

#include "../../include/Lattices/GbFacet.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(const bool ok, const std::string& what, const double& value = 0.0)
{
    std::cout << (ok ? "  [ ok ] " : "  [FAIL] ") << what;
    if (value != 0.0 || !ok) std::cout << "   (" << std::scientific << value << std::fixed << ")";
    std::cout << std::endl;
    if (!ok) ++failures;
}

/*! A periodic point cloud on a cell of side \p a: a K x K grid, jittered in plane so that the
 *  triangulation is not degenerate, corrugated out of plane by \p zAmp, and carrying either a
 *  uniform or a spatially varying displacement.
 *
 *  The in-plane jitter deliberately pushes some nodes outside the base cell (negative
 *  coordinates), because that is what exposes a mishandled periodic offset.
 */
std::vector<std::vector<double>> cloud(const int K, const double a, const double zAmp,
                                       const Eigen::Vector3d& b, const bool varyU)
{
    std::vector<std::vector<double>> pc;
    for (int i = 0; i < K; ++i)
        for (int j = 0; j < K; ++j) {
            const double x = (i + 0.13*std::sin(3.0*j)) * a/K;
            const double y = (j + 0.11*std::cos(2.0*i)) * a/K;
            const double z = zAmp * std::sin(2*M_PI*i/K) * std::cos(2*M_PI*j/K);
            Eigen::Vector3d u = b;
            if (varyU) u = b * (1.0 + 0.7*std::sin(2*M_PI*i/K)*std::cos(2*M_PI*j/K));
            pc.push_back({x, y, z, u(0), u(1), u(2)});
        }
    return pc;
}

std::vector<Eigen::Vector3d> deformed(const std::vector<std::vector<double>>& pc)
{
    std::vector<Eigen::Vector3d> d;
    d.reserve(pc.size());
    for (const auto& r : pc) d.emplace_back(r[0]+r[3], r[1]+r[4], r[2]+r[5]);
    return d;
}

} // namespace

int main()
{
    const double a = 10.0;
    const std::vector<double> box{a,0,0, 0,a,0};
    const Eigen::Vector3d n(0,0,1);
    const Eigen::Vector3d b(0.3,-0.2,0.15);
    const double tol = 1.0e-12;

    std::cout << std::fixed << std::setprecision(4);

    try {

    // ---- 1. A uniform jump gives the rigid translation +-b, exactly, at any shell count.
    for (const double zAmp : {0.0, 1.5})
        for (const int shells : {0,1,2,4}) {
            const auto pc = cloud(5, a, zAmp, b, false);
            const GbFacet f(pc, box, GbFacet::triangulate(deformed(pc), box, n), +1, shells);

            // Stay clear of the corrugation band: a point inside it can legitimately be on either
            // side of the surface, so +z alone would not say which side it is on.
            double worst = 0.0;
            for (const double d : {2.0, 7.0, 25.0})
                for (const double sx : {0.0, 3.3}) {
                    worst = std::max(worst, (f.displacement({sx,1.1, d}) - b).norm());
                    worst = std::max(worst, (f.displacement({sx,1.1,-d}) + b).norm());
                }
            check(worst <= tol, "uniform jump -> rigid +-b, zAmp=" + std::to_string(zAmp)
                                + " shells=" + std::to_string(shells), worst);
        }

    // ---- 2. The representatives tile one period cell exactly once.
    //
    // Catches a face lifted to the wrong periodic image, which the uniform-jump check cannot see:
    // such a face spans the cell and is orientation-reversed, so the signed areas still sum to the
    // cell area while the absolute areas overshoot it, and the reversed faces show up as normals
    // pointing against the rest.
    for (const double zAmp : {0.0, 1.5}) {
        const auto pc = cloud(5, a, zAmp, b, false);
        const GbFacet f(pc, box, GbFacet::triangulate(deformed(pc), box, n), +1, 1);

        const auto& V = f.vertices();
        const auto& F = f.faces();
        double signedArea = 0.0, absArea = 0.0;
        for (int i = 0; i < F.rows(); ++i) {
            const Eigen::Vector3d p0 = V.row(F(i,0)), p1 = V.row(F(i,1)), p2 = V.row(F(i,2));
            const double A = 0.5 * ((p1-p0)(0)*(p2-p0)(1) - (p1-p0)(1)*(p2-p0)(0));
            signedArea += A;
            absArea    += std::abs(A);
        }
        check(std::abs(signedArea - a*a) <= 1e-9,
              "signed projected area == cell area, zAmp=" + std::to_string(zAmp),
              std::abs(signedArea - a*a));
        check(std::abs(absArea - a*a) <= 1e-9,
              "no face lifted to the wrong image, zAmp=" + std::to_string(zAmp),
              std::abs(absArea - a*a));
        check(f.faceNormals().col(2).minCoeff() > 0.0,
              "every face normal faces along the given normal, zAmp=" + std::to_string(zAmp),
              f.faceNormals().col(2).minCoeff());
    }

    // ---- 3. Reversing the normal sense reverses the field, and the solid angle is +-2*pi.
    {
        const auto pc = cloud(5, a, 1.5, b, false);
        const auto topology = GbFacet::triangulate(deformed(pc), box, n);
        const GbFacet fp(pc, box, topology, +1, 4);
        const GbFacet fm(pc, box, topology, -1, 4);
        const Eigen::Vector3d x(1.0, 2.0, 6.0);

        check((fp.displacement(x) + fm.displacement(x)).norm() <= tol,
              "reversing normalSense reverses the field",
              (fp.displacement(x) + fm.displacement(x)).norm());
        check(std::abs(fp.solidAngle(x) - 2*M_PI) <= tol,
              "solid angle == +2*pi on the positive side",
              std::abs(fp.solidAngle(x) - 2*M_PI));
        check(fp.sideOf(x) == 1 && fm.sideOf(x) == -1, "sideOf follows normalSense");
    }

    // ---- 4. A varying jump: the field is exactly periodic, converges in the shell count, and
    //         tends to the mean nodal displacement far from the facet.
    {
        const auto pc = cloud(6, a, 1.5, b, true);
        const auto topology = GbFacet::triangulate(deformed(pc), box, n);
        const Eigen::Vector3d x(2.0, 3.0, 4.0);

        double previousStep = std::numeric_limits<double>::max();
        Eigen::Vector3d previous = Eigen::Vector3d::Zero();
        for (const int shells : {1,2,4,8}) {
            const GbFacet f(pc, box, topology, +1, shells);
            const Eigen::Vector3d u = f.displacement(x);

            const double periodicity = std::max(
                (f.displacement(x + Eigen::Vector3d(a,0,0)) - u).norm(),
                (f.displacement(x + Eigen::Vector3d(0,a,0)) - u).norm());
            check(periodicity <= tol,
                  "field is periodic, shells=" + std::to_string(shells), periodicity);

            if (shells > 1) {
                const double step = (u - previous).norm();
                check(step < previousStep,
                      "shell " + std::to_string(shells) + " refines the previous value", step);
                previousStep = step;
            }
            previous = u;
        }

        const GbFacet f(pc, box, topology, +1, 8);
        const double farError = (f.displacement({2,3,400}) - f.meanDisplacement()).norm();
        check(farError <= 1e-3, "far field tends to the mean nodal displacement", farError);
    }

    // ---- 5. A three-node cell, as small a mesostate as the ensemble produces.  The periodic
    //         Delaunay triangulation of so few points emits triangles whose nodes are exactly
    //         collinear in projection, which have no orientation of their own -- the orientation
    //         has to come from the neighbouring faces instead.
    {
        const std::vector<double> smallBox{8.0,0,0, 0,3.6,0};
        std::vector<std::vector<double>> pc{
            {0.0, 0.0, 0.0,  0.10, 0.00, 0.05},
            {2.7, 1.2, 0.4, -0.05, 0.10, 0.00},
            {5.3, 2.4, -0.4, 0.00,-0.10, 0.05}};
        const GbFacet f(pc, smallBox, GbFacet::triangulate(deformed(pc), smallBox, n), +1, 4);

        check(f.faces().rows() == 6, "three nodes on a torus give 6 faces",
              static_cast<double>(f.faces().rows()));
        const double periodicity =
            (f.displacement({1.0,0.5,9.0}) - f.displacement({9.0,0.5,9.0})).norm();
        check(periodicity <= tol, "three-node field is periodic", periodicity);
        check(f.sideOf({1.0,0.5, 30.0}) ==  1 &&
              f.sideOf({1.0,0.5,-30.0}) == -1, "three-node side test");
    }

    }
    catch (const std::exception& e) {
        std::cout << "  [FAIL] threw: " << e.what() << std::endl;
        ++failures;
    }

    std::cout << (failures ? "testGbFacet FAILED with " : "testGbFacet passed, ")
              << failures << " failure(s)." << std::endl;
    return failures ? -1 : 0;
}
