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

template <int dim>
GbMesoState<dim>::GbMesoState(
    const Gb<dim> &gb,
    const ReciprocalLatticeVector<dim> &axis,
    const std::deque<GbNode<dim>>& engagedNodes,
    const std::vector<LatticeVector<dim>> &mesoStateCslVectors) try:
    /*init*/ GbContinuum<dim>(getMesoStateGbDomain(mesoStateCslVectors),
                              getFacetedSurfaces(gb,mesoStateCslVectors,engagedNodes),
                              gb.nA.cartesian(),
                              getMesoStateBoxDim(mesoStateCslVectors)),
    /*init*/ gb(gb),
    /*init*/ axis(axis),
    /*init*/ mesoStateCslVectors(mesoStateCslVectors),
    /*init*/ engagedNodes(engagedNodes) {}
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

    template<int dim>
    std::pair<typename GbMesoState<dim>::XuPairs, typename GbMesoState<dim>::XuPairs>
    GbMesoState<dim>::getFacetedSurfaces(const Gb<dim>& gb,
                                 const std::vector<LatticeVector<dim>>& mesoStateCslVectors,
                                 const std::deque<GbNode<dim>>& engagedNodes)
    {
        XuPairs xuPairsA, xuPairsB;

        std::vector<LatticeVector<dim>> bicrystalBoxVectors(mesoStateCslVectors);
        bicrystalBoxVectors[0]= 2*mesoStateCslVectors[0];
        VectorDimD shift;
        shift << -0.5-FLT_EPSILON,-FLT_EPSILON,-FLT_EPSILON;

        // Remove the common-mode translation.  Far from the boundary each grain is translated
        // by the mean of its nodal displacements, and only the DIFFERENCE of the two means is
        // physical -- it is the rigid-body translation between the grains, a microscopic degree
        // of freedom of the boundary.  Their average is a rigid translation of the whole
        // bicrystal, which is unobservable under periodic boundaries but is not harmless here:
        // the LAMMPS energy windows are defined in absolute box coordinates, so a drifting
        // configuration shifts which atoms are counted as boundary and which as bulk.  The
        // symmetric split removed it by construction; independent displacements do not.
        //
        // Subtracting the same vector from both displacements leaves every jump
        // u_A - u_B exactly as it was, and moves both deformed surfaces together, so they stay
        // glued.  Only the coincidence points move, by -drift.
        VectorDimD drift= VectorDimD::Zero();
        for(const auto& node : engagedNodes)
            drift+= 0.5*(node.uA()+node.uB());
        if (!engagedNodes.empty())
            drift/= static_cast<double>(engagedNodes.size());

        std::set<OrderedTuplet<dim>> xAIntegerCoordsSet, xBIntegerCoordsSet;
        for(const auto& node : engagedNodes)
        {
            // Read the displacements off the node before wrapping: they are properties of the
            // node, not of which periodic image its atoms are drawn in.  Wrapping moves an atom
            // by a box vector and carries its deformed position with it, which is a shift the
            // two facets then differ by -- an in-plane period, which GbContinuum removes when it
            // glues them.
            const VectorDimD uA= node.uA()-drift;
            const VectorDimD uB= node.uB()-drift;
            VectorDimD xA= node.xA.cartesian();
            VectorDimD xB= node.xB.cartesian();
            LatticeVector<dim>::modulo(xA,bicrystalBoxVectors,shift);
            LatticeVector<dim>::modulo(xB,bicrystalBoxVectors,shift);

            OrderedTuplet<dim> xAIntegerCoords, xBIntegerCoords;
            try {
                xAIntegerCoords << gb.bc.A.latticeVector(xA);
            }
            catch(const std::runtime_error& e) {
                throw std::runtime_error("xA is not a lattice vector of A: " + std::string(e.what()));
            }
            try {
                xBIntegerCoords << gb.bc.B.latticeVector(xB);
            }
            catch(const std::runtime_error& e) {
                throw std::runtime_error("xB is not a lattice vector of B: " + std::string(e.what()));
            }

            // Two engaged nodes may not use the same atom of either grain.
            const bool insertedA= xAIntegerCoordsSet.insert(xAIntegerCoords).second;
            const bool insertedB= xBIntegerCoordsSet.insert(xBIntegerCoords).second;
            if(!insertedA || !insertedB)
                throw std::runtime_error("Clash in constraints.");

            xuPairsA.emplace_back(xA,uA);
            xuPairsB.emplace_back(xB,uB);
        }
        return std::make_pair(xuPairsA,xuPairsB);
    }

    /*-------------------------------------*/
    template<int dim>
    typename GbMesoState<dim>::Relaxations
    GbMesoState<dim>::relaxations(const std::string& lmpLocation,
                                  const std::string& potentialName,
                                  const std::string& configFile,
                                  const double& tetherHalfWidth,
                                  const double& tetherStiffness,
                                  const std::string& tetheredDumpFile,
                                  const std::string& fullDumpFile,
                                  const bool& chainRelaxations) const
    {
        // Only a path to hand on, so the atoms have to come off disk.  Building the state has to
        // happen first if there is not even a path.  The sweep does not come through here: it
        // passes the configuration box() built, and never writes a file at all.
        Configuration configuration;
        if (configFile.empty())
            box("", nullptr, false, nullptr, &configuration);
        else {
            const auto [atoms, cellBox, origin]= read_oILAB_output(configFile);
            configuration.atoms= atoms;
            configuration.box= cellBox;
            configuration.origin= origin;
        }
        return relaxations(lmpLocation, potentialName, configuration, tetherHalfWidth,
                           tetherStiffness, tetheredDumpFile, fullDumpFile, chainRelaxations);
    }

    /*-------------------------------------*/
    template<int dim>
    typename GbMesoState<dim>::Relaxations
    GbMesoState<dim>::relaxations(const std::string& lmpLocation,
                                  const std::string& potentialName,
                                  const Configuration& configuration,
                                  const double& tetherHalfWidth,
                                  const double& tetherStiffness,
                                  const std::string& tetheredDumpFile,
                                  const std::string& fullDumpFile,
                                  const bool& chainRelaxations) const
    {
        const auto& atoms= configuration.atoms;
        const auto& cellBox= configuration.box;
        const auto& origin= configuration.origin;

        Relaxations result;

        if (chainRelaxations && tetherHalfWidth > 0.0)
        {
            // One invocation for both: LAMMPS relaxes against the restraint, reports, releases
            // it, and relaxes again from there.  The tethered structure goes to the first dump
            // and the freely relaxed one to the second, as when the two are run separately.
            double spring= 0.0, unrelaxedEnergy= 0.0, freeEnergy= 0.0;
            const auto tethered= energy(lmpLocation, atoms, cellBox, origin, potentialName, true,
                                        tetheredDumpFile, tetherHalfWidth, tetherStiffness,
                                        &spring, &unrelaxedEnergy,
                                        true, fullDumpFile, &freeEnergy);
            result.density  = tethered.first;
            result.tethered = tethered.second;
            result.spring   = spring;
            result.unrelaxed= unrelaxedEnergy;
            result.full     = freeEnergy;
            return result;
        }

        if (tetherHalfWidth > 0.0 && bothRelaxationsInOneInvocation)
        {
            // Both relaxations in one LAMMPS invocation.  They are still two independent runs --
            // the second starts with its own clear and read_data, from the configuration the
            // construction produced, not from where the tether left the atoms -- so the answers
            // are those of two separate invocations.  What is saved is a process launch, a parse
            // of the potential file and a neighbour-list build, which is most of what an
            // invocation costs; the physics is untouched.
            double spring= 0.0, unrelaxedEnergy= 0.0, freeEnergy= 0.0;
            const auto tethered= energy(lmpLocation, atoms, cellBox, origin, potentialName, true,
                                        tetheredDumpFile, tetherHalfWidth, tetherStiffness,
                                        &spring, &unrelaxedEnergy,
                                        false, fullDumpFile, &freeEnergy, true);
            result.density  = tethered.first;
            result.tethered = tethered.second;
            result.spring   = spring;
            result.unrelaxed= unrelaxedEnergy;
            result.full     = freeEnergy;
            return result;
        }

        // One invocation per relaxation: the free one first, then the tethered one if there is
        // a tether.  Slower by a process launch and a potential parse per state, and the answer
        // is the same -- which is what makes it the check on the shared-invocation path above.
        double ignoredSpring= 0.0, unrelaxed= 0.0;
        const auto full= energy(lmpLocation, atoms, cellBox, origin, potentialName, true,
                                fullDumpFile, 0.0, 1.0, &ignoredSpring, &unrelaxed);
        result.density  = full.first;
        result.full     = full.second;
        result.unrelaxed= unrelaxed;
        result.tethered = result.full;
        result.spring   = 0.0;
        if (tetherHalfWidth > 0.0) {
            double spring= 0.0, unrelaxedAgain= 0.0;
            const auto tethered= energy(lmpLocation, atoms, cellBox, origin, potentialName, true,
                                        tetheredDumpFile, tetherHalfWidth, tetherStiffness,
                                        &spring, &unrelaxedAgain);
            result.tethered= tethered.second;
            result.spring  = spring;
            if (std::abs(tethered.first - result.density) > 1.0e-9)
                throw std::runtime_error("the tethered and free relaxations of one mesostate "
                                         "disagree on its density");
        }
        return result;
    }

    /*-------------------------------------*/
    template<int dim>
    std::tuple<double,double> GbMesoState<dim>::densityEnergy(const std::string& lmpLocation,
                                                             const std::string& potentialName,
                                                             const bool& minimize,
                                                             const std::string& configFile,
                                                             const std::string& minimizedDumpFile,
                                                             const double& tetherHalfWidth,
                                                             const double& tetherStiffness,
                                                             double* springEnergy,
                                                             double* unminimizedEnergy) const
    {
        // Writing the configuration means evaluating the displacement field at every atom, which
        // dominates the cost of a mesostate -- so it is done here only when the caller has not
        // already written one.
        std::string deformedFile= configFile;
        if (deformedFile.empty()) {
            box("temp" + std::to_string(omp_get_thread_num()));
            deformedFile= "temp" + std::to_string(omp_get_thread_num()) + "_reference1.txt";
        }
        std::pair<double,double> densityEnergyPair= energy(lmpLocation,
                                                           deformedFile,
                                                           potentialName,
                                                           minimize,
                                                           minimizedDumpFile,
                                                           tetherHalfWidth,
                                                           tetherStiffness,
                                                           springEnergy,
                                                           unminimizedEnergy);


        return {densityEnergyPair.first,densityEnergyPair.second};

    }


    template<int dim>
    typename std::enable_if<dim==3,void>::type
    GbMesoState<dim>::box(const std::string& name, int* atomsExpelled,
                          const bool& dropUnengagedCoincidences, int* atomsDropped,
                          Configuration* deformedConfiguration) const
    {
        const auto& config= gb.bc.box(mesoStateCslVectors,0);
        std::vector<LatticeVector<3>> boxVectors;
        boxVectors.push_back(this->mesoStateCslVectors[0]);
        boxVectors.push_back(this->mesoStateCslVectors[1]);
        boxVectors.push_back(this->mesoStateCslVectors[2]);

        std::vector<VectorDimD> referenceConfigA, deformedConfigA;
        std::vector<VectorDimD> referenceConfigB, deformedConfigB;
        // Atoms the deformation carries out of their own grain, counted so that a caller can see
        // that it happened: it changes the atom count, and with it the density the energy is
        // reported against.
        int expelled= 0;
        // Species written for each atom: coincidenceType for the atoms the mesostate brings
        // together, the grain's own type for the rest.
        std::vector<int> speciesA, speciesB;

        // Which atoms are the ones brought into coincidence.  They are identified by the same
        // wrapped integer coordinates the clash rule uses, because an atom of the configuration
        // and the node describing it need not be given in the same periodic image.  Marking them
        // is what lets the boundary the construction actually built be picked out of the relaxed
        // structure: after LAMMPS fuses each coincident pair, the survivor still carries this
        // species, so the atoms forming the boundary plane remain identifiable.
        std::vector<LatticeVector<dim>> bicrystalBoxVectors(mesoStateCslVectors);
        bicrystalBoxVectors[0]= 2*mesoStateCslVectors[0];
        VectorDimD nodeShift;
        nodeShift << -0.5-FLT_EPSILON,-FLT_EPSILON,-FLT_EPSILON;
        const auto siteKey= [&bicrystalBoxVectors,&nodeShift]
                            (const Lattice<dim>& lattice, const VectorDimD& x)
        {
            VectorDimD wrapped= x;
            LatticeVector<dim>::modulo(wrapped,bicrystalBoxVectors,nodeShift);
            OrderedTuplet<dim> key;
            key << lattice.latticeVector(wrapped);
            return key;
        };
        std::set<OrderedTuplet<dim>> coincidentA, coincidentB;
        for (const auto& [x,u] : this->xuPairsOfFacetedSurfaces.first)
            coincidentA.insert(siteKey(gb.bc.A,x));
        for (const auto& [x,u] : this->xuPairsOfFacetedSurfaces.second)
            coincidentB.insert(siteKey(gb.bc.B,x));

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
                // The cut above is made at the REFERENCE facet, which is the right cut to make
                // there -- but the displacement can carry an atom that was on its own side of
                // that facet through the deformed surface and into the other grain, where it has
                // no business being.  It is dropped from both configurations rather than from the
                // deformed one alone: the two are the same atoms seen before and after, and a
                // state does not contain an atom that its own deformation expels.
                if (!this->inGrainAAfterDeformation(x)) { ++expelled; continue; }
                referenceConfigA.push_back(wrapIntoBox(latticeVector.cartesian()));
                deformedConfigA.push_back(wrapIntoBox(x));
                speciesA.push_back(coincidentA.count(siteKey(gb.bc.A,latticeVector.cartesian()))
                                   ? coincidenceType : 1);
            }
            else if (&(latticeVector.lattice) == &(gb.bc.B) && this->inGrainB(latticeVector.cartesian())) {
                x= latticeVector.cartesian() + this->displacement(latticeVector.cartesian(),2);
                if (!this->inGrainBAfterDeformation(x)) { ++expelled; continue; }
                referenceConfigB.push_back(wrapIntoBox(latticeVector.cartesian()));
                deformedConfigB.push_back(wrapIntoBox(x));
                speciesB.push_back(coincidentB.count(siteKey(gb.bc.B,latticeVector.cartesian()))
                                   ? coincidenceType : 2);
            }
        }

        if (atomsExpelled) *atomsExpelled= expelled;

        // ---- drop the coincidences this state did not engage ---------------------------
        // The displacement field cannot be aimed at the engaged nodes alone: wherever it happens
        // to close the gap between two other atoms of the two grains, those meet as well.  The
        // configuration then holds coincidences the signature never named, and shows coincidence
        // points the state never chose.
        //
        // A coincidence the state did not engage is removed entirely: both of the atoms that
        // met there are left out, so the site is empty rather than holding the one atom a
        // fusion would have left.
        //
        // This is not free.  Leaving one atom would have been pure bookkeeping -- the survivor
        // is what the overlap removal produces anyway -- but removing both puts a vacancy in the
        // boundary, so the atom count, the density and the relaxed structure all differ from
        // what the same state gives untouched.  That is the point: an unengaged site should look
        // unengaged, and an atom sitting on it looks exactly like one that was engaged.  The
        // cost is that these states are no longer the states the enumeration nominally built,
        // and their energies are not comparable with a run that kept the atoms.
        //
        // A group holding an engaged atom is the state's own coincidence and stays; it loses
        // only the ordinary atoms that drifted onto it.
        int dropped= 0;
        if (dropUnengagedCoincidences)
        {
            const std::size_t countA= deformedConfigA.size();
            std::vector<VectorDimD> deformed(deformedConfigA);
            deformed.insert(deformed.end(), deformedConfigB.begin(), deformedConfigB.end());
            std::vector<int> species(speciesA);
            species.insert(species.end(), speciesB.begin(), speciesB.end());
            const std::size_t total= deformed.size();

            // Grouped with the cutoff the overlap removal uses, so what is grouped here is what
            // LAMMPS would have fused.  The in-plane coordinates are taken once per atom rather
            // than once per pair: the pair loop is then arithmetic, and the periodic images are
            // handled by rounding the difference of those coordinates.
            const auto inPlaneSolver= inPlanePeriods.colPivHouseholderQr();
            std::vector<Eigen::Matrix<double,dim-1,1>> inPlane(total);
            std::vector<VectorDimD> outOfPlane(total);
            for (std::size_t i=0; i<total; ++i) {
                inPlane[i]= inPlaneSolver.solve(deformed[i]);
                outOfPlane[i]= deformed[i] - inPlanePeriods*inPlane[i];
            }

            std::vector<int> parent(total);
            for (std::size_t i=0; i<total; ++i) parent[i]= (int)i;
            const auto root= [&parent](int i)
            { while (parent[i]!=i) { parent[i]= parent[parent[i]]; i= parent[i]; } return i; };

            for (std::size_t i=0; i<total; ++i)
                for (std::size_t j=i+1; j<total; ++j)
                {
                    Eigen::Matrix<double,dim-1,1> difference= inPlane[i]-inPlane[j];
                    for (int k=0; k<dim-1; ++k) difference(k)-= std::round(difference(k));
                    const VectorDimD separation=
                        inPlanePeriods*difference + (outOfPlane[i]-outOfPlane[j]);
                    if (separation.norm() >= lammpsOverlapCutoff) continue;
                    const int ri= root((int)i), rj= root((int)j);
                    if (ri!=rj) parent[ri]= rj;
                }

            std::vector<int> groupSize(total,0), groupEngaged(total,0);
            for (std::size_t i=0; i<total; ++i) {
                const int r= root((int)i);
                ++groupSize[r];
                if (species[i]==coincidenceType) ++groupEngaged[r];
            }

            std::vector<char> keep(total,1);
            for (std::size_t i=0; i<total; ++i) {
                const int r= root((int)i);
                if (groupSize[r] < 2) continue;                    // not a coincidence at all
                if (groupEngaged[r] > 0) {
                    // The state's own coincidence: keep what it engaged, drop what drifted in.
                    if (species[i]!=coincidenceType) keep[i]= 0;
                }
                else keep[i]= 0;      // nobody engaged this: the whole coincidence goes
            }

            // An atom the field carries onto the boundary surface without meeting another one
            // there is the same problem in a milder form.  It forms no coincidence, so the
            // grouping above leaves it -- but it stands on the surface exactly where an engaged
            // atom would, so the boundary shows a point of contact the state never chose, and
            // the site it fills should be empty.  The side tests keep such an atom deliberately,
            // on the reasoning that atoms on the facet are the ones the grains are brought
            // together at; that reasoning is about coincidences and does not cover this.
            //
            // Only atoms of the state's own coincidences are allowed to stand on the surface.
            // Judged with the cutoff the overlap removal uses, not with a floating-point
            // epsilon: an atom a few thousandths of an angstrom off the surface is standing on
            // the boundary as surely as one exactly on it, and the two have to be treated alike
            // or the rule catches only the ones that happen to land exactly.
            for (std::size_t i=0; i<total; ++i) {
                if (!keep[i] || species[i]==coincidenceType) continue;
                if (this->onDeformedSurface(deformed[i], lammpsOverlapCutoff)) keep[i]= 0;
            }

            std::vector<VectorDimD> keptReferenceA, keptDeformedA, keptReferenceB, keptDeformedB;
            std::vector<int> keptSpeciesA, keptSpeciesB;
            for (std::size_t i=0; i<total; ++i) {
                if (!keep[i]) { ++dropped; continue; }
                if (i < countA) {
                    keptReferenceA.push_back(referenceConfigA[i]);
                    keptDeformedA.push_back(deformedConfigA[i]);
                    keptSpeciesA.push_back(speciesA[i]);
                }
                else {
                    const std::size_t k= i-countA;
                    keptReferenceB.push_back(referenceConfigB[k]);
                    keptDeformedB.push_back(deformedConfigB[k]);
                    keptSpeciesB.push_back(speciesB[k]);
                }
            }
            referenceConfigA= std::move(keptReferenceA);
            deformedConfigA=  std::move(keptDeformedA);
            speciesA=         std::move(keptSpeciesA);
            referenceConfigB= std::move(keptReferenceB);
            deformedConfigB=  std::move(keptDeformedB);
            speciesB=         std::move(keptSpeciesB);
        }
        if (atomsDropped) *atomsDropped= dropped;

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

        // The deformed configuration in memory, in the layout read_oILAB_output() produces, so
        // that the coincidence count and LAMMPS can take it directly.  The three cell vectors go
        // in as columns and the species leads each row, exactly as the file spells them.
        if (deformedConfiguration)
        {
            Configuration& out= *deformedConfiguration;
            out.atoms.resize(nAtoms,5);
            out.atoms.setZero();
            for (std::size_t i=0; i<deformedConfigA.size(); ++i) {
                out.atoms(i,0)= speciesA[i];
                out.atoms.row(i).segment(1,3)= deformedConfigA[i].transpose();
                out.atoms(i,4)= 0.05;
            }
            const std::size_t offset= deformedConfigA.size();
            for (std::size_t i=0; i<deformedConfigB.size(); ++i) {
                out.atoms(offset+i,0)= speciesB[i];
                out.atoms.row(offset+i).segment(1,3)= deformedConfigB[i].transpose();
                out.atoms(offset+i,4)= 0.05;
            }
            out.box.col(0)= deformedScale*2*boxVectors[0].cartesian();
            out.box.col(1)= boxVectors[1].cartesian();
            out.box.col(2)= boxVectors[2].cartesian();
            out.origin= -deformedScale*boxVectors[0].cartesian();
        }

        // No name, no files.  A survey pass wants the numbers, not the configuration, and
        // writing two extended-XYZ files at fifteen digits is about a quarter of this function.
        if (name.empty()) return;

        std::string referenceFile= name + "_reference0.txt";
        std::string deformedFile= name + "_reference1.txt";

        // Seventeen digits throughout, atoms and cell vectors alike, because that is what a
        // double survives a round trip through decimal at.  At fifteen the file was not the
        // configuration the construction produced but a close decimal neighbour of it, and
        // reading it back moved the freely relaxed energy of 42 of 1119 states by up to 5.6e-6
        // eV -- the free minimisation being soft enough to amplify the last bits, while the
        // as-constructed and tethered energies were unmoved.  The sweep does not read these
        // files any longer, so this matters for anything that re-examines a state afterwards:
        // at seventeen digits it reproduces the sweep exactly, checked over all 1119 states.
        std::ofstream reference, deformed;
        reference.open(referenceFile);
        deformed.open(deformedFile);
        if (!reference || !deformed) std::cerr << "Unable to open files";
        reference << nAtoms << std::endl; deformed << nAtoms << std::endl;
        reference << "Lattice=\" "; deformed << "Lattice=\" ";

        reference << std::setprecision(17) << (referenceScale*2*boxVectors[0].cartesian()).transpose() << " ";
        deformed << std::setprecision(17) << (deformedScale*2*boxVectors[0].cartesian()).transpose() << " ";
        reference << std::setprecision(17) << (boxVectors[1].cartesian()).transpose() << " ";
        deformed << std::setprecision(17) << (boxVectors[1].cartesian()).transpose() << " ";
        reference << std::setprecision(17) << (boxVectors[2].cartesian()).transpose();
        deformed << std::setprecision(17) << (boxVectors[2].cartesian()).transpose();
        reference << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"F T T\" origin=\" "; deformed << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\" F T T\" origin=\" ";
        reference << std::setprecision(17) << (-referenceScale * boxVectors[0].cartesian()).transpose() << "\"" << std::endl;
        deformed << std::setprecision(17) << (-deformedScale * boxVectors[0].cartesian()).transpose() << "\"" << std::endl;

        for(std::size_t i=0; i<referenceConfigA.size(); ++i)
            reference << speciesA[i] << " " << std::setprecision(17)
                      << referenceConfigA[i].transpose() << "  " << 0.05 << std::endl;
        for(std::size_t i=0; i<referenceConfigB.size(); ++i)
            reference << speciesB[i] << " " << std::setprecision(17)
                      << referenceConfigB[i].transpose() << "  " << 0.05 << std::endl;
        for(std::size_t i=0; i<deformedConfigA.size(); ++i)
            deformed << speciesA[i] << " " << std::setprecision(17)
                     << deformedConfigA[i].transpose() << "  " << 0.05 << std::endl;
        for(std::size_t i=0; i<deformedConfigB.size(); ++i)
            deformed << speciesB[i] << " " << std::setprecision(17)
                     << deformedConfigB[i].transpose() << "  " << 0.05 << std::endl;

        reference.close();
        deformed.close();
    }

 } // namespace oILAB

#endif
