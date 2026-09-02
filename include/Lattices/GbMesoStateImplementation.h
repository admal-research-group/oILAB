//
// Created by Nikhil Chandra Admal on 5/27/24.
//

#ifndef OILAB_GBMESOSTATEIMPLEMENTATION_H
#define OILAB_GBMESOSTATEIMPLEMENTATION_H

#include <cassert>
#include "../IO/Lammps.h"
#include <cfloat>
#include <cmath>
#include <iostream>
#include <set>
#include "GbMesoState.h"
#include "OrderedTuplet.h"

namespace oILAB {

template<int dim>
std::vector<double> GbMesoState<dim>::getMesoStateBoxDim(const std::vector<LatticeVector<dim>>& cslVectors)
{
    static_assert(dim == 3, "getMesoStateBoxDim requires dim=3");
    // Store the actual in-plane period vectors (Cartesian) so GbFacet can build correct
    // periodic images.  Format: {p1_x, p1_y, p1_z, p2_x, p2_y, p2_z}.
    const auto p1 = cslVectors[1].cartesian();
    const auto p2 = cslVectors[2].cartesian();
    return {p1(0), p1(1), p1(2), p2(0), p2(1), p2(2)};
}

template <int dim>
GbMesoState<dim>::GbMesoState(
    const Gb<dim> &gb,
    const ReciprocalLatticeVector<dim> &axis,
    const std::deque<std::pair<LatticeVector<dim>, VectorDimD>>& engagedTsPairs,
    const std::vector<LatticeVector<dim>> &mesoStateCslVectors) try:
    /*init*/ GbContinuum<dim>(getMesoStateGbDomain(mesoStateCslVectors),
                              getFacetedSurfaces(gb,mesoStateCslVectors,engagedTsPairs),
                              gb.nA.cartesian(),
                              getMesoStateBoxDim(mesoStateCslVectors)),
    /*init*/ gb(gb),
    /*init*/ axis(axis),
    /*init*/ mesoStateCslVectors(mesoStateCslVectors),
    /*init*/ engagedTsPairs(engagedTsPairs) {}
    catch(std::runtime_error& e)
    {
        throw;
    }


    template<int dim>
    Eigen::Matrix<double, dim,dim-1> GbMesoState<dim>::getMesoStateGbDomain(const std::vector<LatticeVector<dim>>& mesoStateCslVectors)
    {
        Eigen::Matrix<double, dim,dim-1> mesoStateGbDomain;
        for(int i=1; i<dim; ++i)
            mesoStateGbDomain.col(i-1)= mesoStateCslVectors[i].cartesian();
        return mesoStateGbDomain;
    }

    template<int dim>
    typename GbMesoState<dim>::NodePlacement
    GbMesoState<dim>::nodePlacement(const Gb<dim>& gb,
                                    const std::vector<LatticeVector<dim>>& mesoStateCslVectors,
                                    const LatticeVector<dim>& t,
                                    const VectorDimD& s)
    {
        std::vector<LatticeVector<dim>> bicrystalBoxVectors(mesoStateCslVectors);
        bicrystalBoxVectors[0]= 2*mesoStateCslVectors[0];
        VectorDimD shift;
        shift << -0.5-FLT_EPSILON,-FLT_EPSILON,-FLT_EPSILON;

        NodePlacement node;
        node.u= t.cartesian()/2;
        node.xA= s-node.u;
        node.xB= s+node.u;

        // modulo w.r.t the bicrystal box
        LatticeVector<dim>::modulo(node.xA,bicrystalBoxVectors,shift);
        LatticeVector<dim>::modulo(node.xB,bicrystalBoxVectors,shift);

        // A pair whose node misses its lattice cannot take part in any mesostate.  Throwing lets
        // the caller drop that one pair; the run used to exit(0) here, taking the whole sweep
        // down with a success status.
        try {
            node.siteA << gb.bc.A.latticeVector(node.xA);
        }
        catch(const std::runtime_error& e) {
            throw std::runtime_error("xA is not a lattice vector of A: " + std::string(e.what()));
        }
        try {
            node.siteB << gb.bc.B.latticeVector(node.xB);
        }
        catch(const std::runtime_error& e) {
            throw std::runtime_error("xB is not a lattice vector of B: " + std::string(e.what()));
        }
        return node;
    }

    template<int dim>
    std::pair<typename GbMesoState<dim>::XuPairs, typename GbMesoState<dim>::XuPairs>
    GbMesoState<dim>::getFacetedSurfaces(const Gb<dim>& gb,
                                 const std::vector<LatticeVector<dim>>& mesoStateCslVectors,
                                 const std::deque<std::pair<LatticeVector<dim>,VectorDimD>>& engagedTsPairs)
    {
        XuPairs xuPairsA, xuPairsB;

        std::set<OrderedTuplet<dim>> xAIntegerCoordsSet, xBIntegerCoordsSet;
        for(const auto& [t,s] : engagedTsPairs)
        {
            const NodePlacement node= nodePlacement(gb,mesoStateCslVectors,t,s);

            // Two engaged pairs may not put nodes on the same site of either lattice.  A caller
            // that enumerates over nodePlacement()'s sites never reaches this throw; it is kept
            // as the assertion that the enumeration and the construction still agree.
            const bool insertedA= xAIntegerCoordsSet.insert(node.siteA).second;
            const bool insertedB= xBIntegerCoordsSet.insert(node.siteB).second;
            if(!insertedA || !insertedB)
                throw std::runtime_error("Clash in constraints.");

            xuPairsA.emplace_back(node.xA,node.u);
            xuPairsB.emplace_back(node.xB,-node.u);
        }
        return std::make_pair(xuPairsA,xuPairsB);
    }

    /*-------------------------------------*/
    template<int dim>
    std::tuple<double,double> GbMesoState<dim>::densityEnergy(const std::string& lmpLocation,
                                                             const std::string& potentialName,
                                                             const bool& minimize) const
    {
        box("temp" + std::to_string(omp_get_thread_num()));
        std::pair<double,double> densityEnergyPair= energy(lmpLocation,
                                                           "temp" + std::to_string(omp_get_thread_num()) + "_reference1.txt",
                                                           potentialName,
                                                           minimize);


        return {densityEnergyPair.first,densityEnergyPair.second};

    }


    template<int dim>
    typename std::enable_if<dim==3,void>::type
    GbMesoState<dim>::box(const std::string& name) const
    {
        const auto& config= gb.bc.box(mesoStateCslVectors,0);
        std::vector<LatticeVector<3>> boxVectors;
        boxVectors.push_back(this->mesoStateCslVectors[0]);
        boxVectors.push_back(this->mesoStateCslVectors[1]);
        boxVectors.push_back(this->mesoStateCslVectors[2]);

        std::vector<VectorDimD> referenceConfigA, deformedConfigA;
        std::vector<VectorDimD> referenceConfigB, deformedConfigB;

        // Fold positions back into the box along the two periodic directions.  The header written
        // below declares PBC="F T T": box vectors 1 and 2 are periodic, while the first spans the
        // two grains and is a free surface.  Only the in-plane part may be wrapped -- folding the
        // out-of-plane direction would carry an atom of one grain into the other.
        //
        // This matters for the deformed configuration: a node's displacement can carry it past a
        // period boundary, and LAMMPS reads that file, so an atom left outside the box it is
        // given is an error rather than a wrapped image.  The reference configuration is folded
        // by the same rule, so that both files hold every atom inside the box they declare.
        Eigen::Matrix<double,dim,dim-1> inPlanePeriods;
        for (int i=1; i<dim; ++i)
            inPlanePeriods.col(i-1)= boxVectors[i].cartesian();

        const auto wrapIntoBox= [&inPlanePeriods](const VectorDimD& x)
        {
            // The least-squares solve returns the coordinates of x in the (generally
            // non-orthogonal) in-plane basis; what it cannot represent is the out-of-plane part,
            // which stays in the residual and is therefore left untouched.
            const Eigen::Matrix<double,dim-1,1> coordinates=
                inPlanePeriods.colPivHouseholderQr().solve(x);
            // Snap before flooring: an atom that belongs at coordinate 1 -- the periodic image of
            // 0 -- routinely lands a rounding error below it, and a bare floor() would leave it
            // sitting on the far face of the box instead of at its origin.
            Eigen::Matrix<double,dim-1,1> whole;
            for (int i=0; i<dim-1; ++i)
                whole(i)= std::floor(coordinates(i) + FLT_EPSILON);
            return VectorDimD(x - inPlanePeriods*whole);
        };

        for (const auto &latticeVector: config) {
            VectorDimD x;
            // Both lattices fill the whole box, so each has to be cut back to the grain it
            // actually occupies: an atom of A above facet A, or of B below facet B, lies on the
            // far side of the boundary and is discarded.  Without this the two grains
            // interpenetrate and the box holds roughly twice the atoms it should.
            if (&(latticeVector.lattice) == &(gb.bc.A) && this->inGrainA(latticeVector.cartesian())) {
                x= latticeVector.cartesian() + this->displacement(latticeVector.cartesian(),1);
                referenceConfigA.push_back(wrapIntoBox(latticeVector.cartesian()));
                deformedConfigA.push_back(wrapIntoBox(x));
            }
            else if (&(latticeVector.lattice) == &(gb.bc.B) && this->inGrainB(latticeVector.cartesian())) {
                x= latticeVector.cartesian() + this->displacement(latticeVector.cartesian(),2);
                referenceConfigB.push_back(wrapIntoBox(latticeVector.cartesian()));
                deformedConfigB.push_back(wrapIntoBox(x));
            }
        }

        int nAtoms= referenceConfigA.size()+referenceConfigB.size();

        // The out-of-plane direction is the one the two grains stack along, so an atom carried
        // past a face by its displacement cannot be folded back the way the in-plane ones are --
        // that would drop it into the other grain.  The box is stretched to hold it instead.
        //
        // The stretch is symmetric about the boundary: energy() rebuilds the LAMMPS cell as
        // [-w/2, +w/2] from the box width alone and ignores the origin written here, so a box
        // grown to one side would put those same atoms back outside the cell LAMMPS reads.
        //
        // Returns the factor the first box vector (and the origin with it) is scaled by; exactly
        // 1 when nothing sits outside, so a configuration that never needed the room is written
        // byte-for-byte as before.
        const auto outOfPlaneScale= [&boxVectors](const std::vector<VectorDimD>& grainA,
                                                  const std::vector<VectorDimD>& grainB)
        {
            Eigen::Matrix<double,dim,dim> basis;
            basis.col(0)= 2*boxVectors[0].cartesian();
            for (int i=1; i<dim; ++i)
                basis.col(i)= boxVectors[i].cartesian();
            // Row 0 of the inverse reads off the coefficient along the first box vector alone,
            // which is the only one that can be out of range once the in-plane fold has run.
            const Eigen::Matrix<double,1,dim> outOfPlaneCoordinate= basis.inverse().row(0);
            const VectorDimD origin= -boxVectors[0].cartesian();

            double overshoot= 0.0;
            for (const auto* positions : {&grainA, &grainB})
                for (const auto& position : *positions)
                {
                    const double a= outOfPlaneCoordinate.dot(position-origin);
                    overshoot= std::max(overshoot, std::max(-a, a-1.0));
                }
            return (overshoot > 0.0 ? 1.0 + 2.0*overshoot : 1.0);
        };
        const double referenceScale= outOfPlaneScale(referenceConfigA, referenceConfigB);
        const double deformedScale=  outOfPlaneScale(deformedConfigA,  deformedConfigB);

        std::string referenceFile= name + "_reference0.txt";
        std::string deformedFile= name + "_reference1.txt";

        std::ofstream reference, deformed;
        reference.open(referenceFile);
        deformed.open(deformedFile);
        if (!reference || !deformed) std::cerr << "Unable to open files";
        reference << nAtoms << std::endl; deformed << nAtoms << std::endl;
        reference << "Lattice=\" "; deformed << "Lattice=\" ";

        reference << std::setprecision(15) << (referenceScale*2*boxVectors[0].cartesian()).transpose() << " ";
        deformed << std::setprecision(15) << (deformedScale*2*boxVectors[0].cartesian()).transpose() << " ";
        reference << std::setprecision(15) << (boxVectors[1].cartesian()).transpose() << " ";
        deformed << std::setprecision(15) << (boxVectors[1].cartesian()).transpose() << " ";
        reference << std::setprecision(15) << (boxVectors[2].cartesian()).transpose();
        deformed << std::setprecision(15) << (boxVectors[2].cartesian()).transpose();
        reference << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"F T T\" origin=\" "; deformed << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\" F T T\" origin=\" ";
        reference << std::setprecision(15) << (-referenceScale * boxVectors[0].cartesian()).transpose() << "\"" << std::endl;
        deformed << std::setprecision(15) << (-deformedScale * boxVectors[0].cartesian()).transpose() << "\"" << std::endl;

        for(const auto& position : referenceConfigA)
            reference << 1 << " " << std::setprecision(15) << position.transpose() << "  " << 0.05 << std::endl;
        for(const auto& position : referenceConfigB)
            reference << 2 << " " << std::setprecision(15) << position.transpose() << "  " << 0.05 << std::endl;
        for(const auto& position : deformedConfigA)
            deformed << 1 << " " << std::setprecision(15) << position.transpose() << "  " << 0.05 << std::endl;
        for(const auto& position : deformedConfigB)
            deformed << 2 << " " << std::setprecision(15) << position.transpose() << "  " << 0.05 << std::endl;

        reference.close();
        deformed.close();
    }

 } // namespace oILAB

#endif
