/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2011 Klaus Spanderen

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file fdsimpleextoustorageengine.cpp
    \brief Finite Differences extended OU engine for simple storage options
*/

#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/experimental/processes/extendedornsteinuhlenbeckprocess.hpp>
#include <ql/experimental/finitedifferences/fdmlinearoplayout.hpp>
#include <ql/experimental/finitedifferences/uniform1dmesher.hpp>
#include <ql/experimental/finitedifferences/fdminnervaluecalculator.hpp>
#include <ql/experimental/finitedifferences/fdmmeshercomposite.hpp>
#include <ql/experimental/finitedifferences/fdmsimpleprocess1dmesher.hpp>
#include <ql/experimental/finitedifferences/fdsimplebsswingengine.hpp>
#include <ql/experimental/finitedifferences/fdmsimplestoragecondition.hpp>
#include <ql/experimental/finitedifferences/fdmstepconditioncomposite.hpp>
#include <ql/experimental/finitedifferences/fdsimpleextoustorageengine.hpp>
#include <ql/experimental/finitedifferences/fdmsolverdesc.hpp>
#include <ql/experimental/finitedifferences/fdmbackwardsolver.hpp>
#include <ql/experimental/finitedifferences/fdmsimple2dextousolver.hpp>


namespace QuantLib {

    namespace {
        class FdmUnderlyingValue : public FdmInnerValueCalculator {
          public:
            FdmUnderlyingValue(const boost::shared_ptr<FdmMesher>& mesher)
            : mesher_(mesher) { }

            Real innerValue(const FdmLinearOpIterator& iter) {
                return std::exp(mesher_->location(iter, 0));
            }
            Real avgInnerValue(const FdmLinearOpIterator& iter) {
                return innerValue(iter);
            }
          private:
            const boost::shared_ptr<FdmMesher> mesher_;
        };

        class FdmStorageValue : public FdmInnerValueCalculator {
          public:
            FdmStorageValue(const boost::shared_ptr<FdmMesher>& mesher)
            : mesher_(mesher) { }

            Real innerValue(const FdmLinearOpIterator& iter) {
                const Real s = std::exp(mesher_->location(iter, 0));
                const Real v = mesher_->location(iter, 1);
                return s*v;
            }
            Real avgInnerValue(const FdmLinearOpIterator& iter) {
                return innerValue(iter);
            }
          private:
            const boost::shared_ptr<FdmMesher> mesher_;

        };
    }

    FdSimpleExtOUStorageEngine::FdSimpleExtOUStorageEngine(
            const boost::shared_ptr<ExtendedOrnsteinUhlenbeckProcess>& process,
            const boost::shared_ptr<YieldTermStructure>& rTS,
            Size tGrid, Size xGrid,
            const FdmSchemeDesc& schemeDesc)
    : process_(process),
      rTS_  (rTS),
      tGrid_(tGrid),
      xGrid_(xGrid),
      schemeDesc_(schemeDesc) {
    }

    void FdSimpleExtOUStorageEngine::calculate() const {

        QL_REQUIRE(arguments_.exercise->type() == Exercise::Bermudan,
                   "Bermudan exercise supported only");

        // 1. Layout
        const Size yGrid = arguments_.capacity/arguments_.changeRate+1;
        std::vector<Size> dim;
        dim.push_back(xGrid_);
        dim.push_back(yGrid);
        printf("yGrid is %i\n",yGrid);
        const boost::shared_ptr<FdmLinearOpLayout> layout(
                                              new FdmLinearOpLayout(dim));

        // 2. Mesher
        const Time maturity
            = rTS_->dayCounter().yearFraction(rTS_->referenceDate(),
                                              arguments_.exercise->lastDate());

        const boost::shared_ptr<Fdm1dMesher> xMesher(
                     new FdmSimpleProcess1dMesher(xGrid_, process_, maturity));

        const boost::shared_ptr<Fdm1dMesher> storageMesher(
            new Uniform1dMesher(0, arguments_.capacity,
                                arguments_.capacity/arguments_.changeRate+1));

        std::vector<boost::shared_ptr<Fdm1dMesher> > meshers;
        meshers.push_back(xMesher);
        meshers.push_back(storageMesher);
        boost::shared_ptr<FdmMesher> mesher (
                                     new FdmMesherComposite(layout, meshers));

        // 3. Calculator
        boost::shared_ptr<FdmInnerValueCalculator> storageCalculator(
                                                  new FdmStorageValue(mesher));

        // 4. Step conditions
        std::list<boost::shared_ptr<StepCondition<Array> > > stepConditions;
        std::list<std::vector<Time> > stoppingTimes;

        // 4.1 Bermudan step conditions
        std::vector<Time> exerciseTimes;
        for (Size i=0; i<arguments_.exercise->dates().size(); ++i) {
            const Time t = rTS_->dayCounter()
                           .yearFraction(rTS_->referenceDate(),
                                         arguments_.exercise->dates()[i]);

            QL_REQUIRE(t >= 0, "exercise dates must not contain past date");
            exerciseTimes.push_back(t);
        }
        stoppingTimes.push_back(exerciseTimes);
        printf("number of exervise times are %i\n",exerciseTimes.size());

        boost::shared_ptr<FdmInnerValueCalculator> underlyingCalculator(
                                               new FdmUnderlyingValue(mesher));

        stepConditions.push_back(boost::shared_ptr<StepCondition<Array> >(
            new FdmSimpleStorageCondition(exerciseTimes,
                                          mesher, underlyingCalculator,
                                          arguments_.changeRate)));

        boost::shared_ptr<FdmStepConditionComposite> conditions(
                new FdmStepConditionComposite(stoppingTimes, stepConditions));

        // 5. Boundary conditions
        std::vector<boost::shared_ptr<FdmDirichletBoundary> > boundaries;

        // 6. Solver
        FdmSolverDesc solverDesc = { mesher, boundaries, conditions,
                                     storageCalculator, maturity, tGrid_, 0 };

        boost::shared_ptr<FdmSimple2dExtOUSolver> solver(
                new FdmSimple2dExtOUSolver(
                           Handle<ExtendedOrnsteinUhlenbeckProcess>(process_),
                           rTS_, solverDesc, schemeDesc_));

        const Real x = process_->x0();
        const Real y = arguments_.load;

        results_.value = solver->valueAt(x, y);
    }
}
