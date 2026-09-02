//
// Created by Nikhil Chandra Admal on 5/26/24.
//

#ifndef OILAB_GBMESOSTATEENSEMBLEIMPLEMENTATION_H
#define OILAB_GBMESOSTATEENSEMBLEIMPLEMENTATION_H

#include <cassert>
#include <deque>
#include "LatticeVector.h"
#include "../Utilities/randomInteger.h"

namespace oILAB {
template <int dim>
GbMesoStateEnsemble<dim>::GbMesoStateEnsemble(const Gb<dim>& gb,
                                              const ReciprocalLatticeVector<dim>& axis,
                                              std::vector<LatticeVector<dim>>& ensembleCslVectors,
                                              const double& tMax,
                                              const double& sPerpMax,
                                              const GbShiftSearch& search,
                                              const double& tPerpMax,
                                              const bool& oneTranslationPerSite,
                                              const std::string& filename):
/*init*/ GbShifts<dim>(gb, axis,
                       std::vector<LatticeVector<dim>>(ensembleCslVectors.begin() + 1,
                                                       ensembleCslVectors.end()),
                       tMax,
                       sPerpMax,
                       search,
                       tPerpMax,
                       oneTranslationPerSite,
                       filename),
/*init*/ ensembleCslVectors(ensembleCslVectors){
        std::cout << "--------------------GBMesoStateEnsemble class construction "
                     "---------------------------"
                  << std::endl;
        std::cout << std::endl;

        std::cout << "Ensemble CSL vectors:" << std::endl;
        for (const auto &latticeVector : ensembleCslVectors)
          std::cout << latticeVector.cartesian().transpose() << std::endl;
        std::cout << std::endl;
    }

    /*-------------------------------------*/
    template<int dim>
    std::deque<std::pair<LatticeVector<dim>,typename GbMesoStateEnsemble<dim>::VectorDimD>>
        GbMesoStateEnsemble<dim>::getEngagedTsPairs(const std::vector<std::pair<LatticeVector<dim>,typename GbMesoStateEnsemble<dim>::VectorDimD>>& tShiftPairs,
                                                    const typename GbMesoStateEnsemble<dim>::Constraints& constraints)
    {
        using VectorDimD = typename GbMesoStateEnsemble<dim>::VectorDimD;
        assert(tShiftPairs.size() == constraints.size());

        std::deque<std::pair<LatticeVector<dim>,VectorDimD>> engagedTsPairs;
        for(int i=0; i<tShiftPairs.size(); ++i) {
            if (constraints(i) == 1) {
                const auto& t= tShiftPairs[i].first;
                auto s= tShiftPairs[i].second;
                engagedTsPairs.push_back(std::make_pair(t,s));
            }
        }
        return engagedTsPairs;
    }

    /*-------------------------------------*/
    template<int dim>
    std::map<typename GbMesoStateEnsemble<dim>::Constraints,GbMesoState<dim>> GbMesoStateEnsemble<dim>::collectMesoStates(const std::string& filename) const
    {
        std::deque<Constraints> constraintsEnsemble(enumerateConstraints(this->tShiftPairs.size()));
        std::cout << "Number of mesostates in the ensemble = " << constraintsEnsemble.size() << std::endl;
        std::cout << std::endl;
        std::cout << "------------------------------" << std::endl;
        std::map<Constraints,GbMesoState<dim>> mesoStates;

        int count= -1;
        for(const Constraints& constraints : constraintsEnsemble)
        {
            try {
                if (constraints(0) != 1) continue;
                //mesoStates.emplace_back(constructMesoState(constraints));
                mesoStates.emplace(constraints,constructMesoState(constraints));
                count++;
                std::cout << "Constructing mesostate " << count << " of " << constraintsEnsemble.size() << std::endl;
                std::cout << "Mesostate signature:  " << constraints.transpose() << std::endl;
                if (!filename.empty())
                    //mesoStates.back().box(filename + std::to_string(count));
                    mesoStates.at(constraints).box(filename + std::to_string(count));
            }
            catch(std::runtime_error& e)
            {
                std::cout << e.what() << std::endl;
            }

        }
        return mesoStates;
    }

    /*-------------------------------------*/
    template<int dim>
    GbMesoState<dim> GbMesoStateEnsemble<dim>::constructMesoState(const Constraints& constraints) const {
        auto engagedTsPairs(getEngagedTsPairs(this->tShiftPairs,constraints));
        try {
            GbMesoState<dim> mesostate(this->gb, this->axis, engagedTsPairs, ensembleCslVectors);
            return mesostate;
        }
        catch(std::runtime_error& e)
        {
            throw(e);
        }
    }

    /*-------------------------------------*/
    template<int dim>
    std::deque<typename GbMesoStateEnsemble<dim>::Constraints> GbMesoStateEnsemble<dim>::enumerateConstraints(const int& size)
    {
        std::deque<Constraints> constraintsEnsemble;

        auto tuples= XTuplet::generate_tuples(2,size);

        for (auto& tuple : tuples) {
            // always engage the (t,s)=(0,0) constraint
            if (tuple(0) == 1)
                constraintsEnsemble.push_back(tuple);
        }
        return constraintsEnsemble;
    }

    /*-------------------------------------*/
    template<int dim>
    typename GbMesoStateEnsemble<dim>::Constraints GbMesoStateEnsemble<dim>::sampleNewState(const Constraints& currentConstraints,
                                                                                            const bool& randomize) const
    {
        // new mesostate construction
        Constraints newConstraints(this->tShiftPairs.size());
        bool msConstructionSuccess = false;
        while (!msConstructionSuccess) {
            try {
                newConstraints = currentConstraints;
                // alter the constraints
                if (randomize){
                    for (auto &elem: newConstraints)
                        elem = random<int>(0, 1);
                } else {
                    int randomSpot = random<int>(0, this->tShiftPairs.size() - 1);
                    newConstraints(randomSpot) = random<int>(0, 1);
                }
                GbMesoState<dim> temp(this->gb,
                                      this->axis,
                                      getEngagedTsPairs(this->tShiftPairs, newConstraints),
                                      ensembleCslVectors);
                msConstructionSuccess = true;
            }
            catch (std::runtime_error &e) {
                throw(e);
            }
        }
        return newConstraints;
    }

        /*-------------------------------------*/
    template<int dim>
    typename GbMesoStateEnsemble<dim>::Constraints GbMesoStateEnsemble<dim>::initializeState() const
    {
        Constraints initialConstraints(this->tShiftPairs.size());
        initialConstraints.setZero();
        return initialConstraints;
    }

} // namespace oILAB
#endif