/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2013 Peter Caspers

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

/*! \file gaussian1dswaptionengine.hpp
    \brief
*/

#ifndef quantlib_pricers_gaussian1d_swaption_hpp
#define quantlib_pricers_gaussian1d_swaption_hpp

#include <ql/instruments/swaption.hpp>
#include <ql/pricingengines/genericmodelengine.hpp>
#include <ql/models/shortrate/onefactormodels/gaussian1dmodel.hpp>
#include <ql/math/interpolations/cubicinterpolation.hpp>
#include <ql/payoff.hpp>


namespace QuantLib {

    //! One factor model swaption engine
    /*! \ingroup swaptionengines

        All fixed coupons with start date greater or equal to the respective
        option expiry are considered to be
        part of the exercise into right.

        \warning Cash settled swaptions are not supported
    */

    template<class Impl>
    class Gaussian1dSwaptionEngine
        : public GenericModelEngine<Gaussian1dModel<Impl>, Swaption::arguments,
                                    Swaption::results> {
      public:
        enum Probabilities {
            None,
            Naive,
            Digital
        };

        Gaussian1dSwaptionEngine(
            const boost::shared_ptr<Gaussian1dModel<Impl> > &model,
            const int integrationPoints = 64, const Real stddevs = 7.0,
            const bool extrapolatePayoff = true,
            const bool flatPayoffExtrapolation = false,
            const Handle<YieldTermStructure> &discountCurve =
                Handle<YieldTermStructure>(),
            const Probabilities probabilities = None)
            : GenericModelEngine<Gaussian1dModel<Impl>, Swaption::arguments,
                                 Swaption::results>(model),
              integrationPoints_(integrationPoints), stddevs_(stddevs),
              extrapolatePayoff_(extrapolatePayoff),
              flatPayoffExtrapolation_(flatPayoffExtrapolation),
              discountCurve_(discountCurve), probabilities_(probabilities) {

            if (!discountCurve_.empty())
                this->registerWith(discountCurve_);
        }

        void calculate() const;

      private:
        const int integrationPoints_;
        const Real stddevs_;
        const bool extrapolatePayoff_, flatPayoffExtrapolation_;
        const Handle<YieldTermStructure> discountCurve_;
        const Probabilities probabilities_;
    };

// implementation

    template <class Impl>
    void Gaussian1dSwaptionEngine<Impl>::calculate() const {

        QL_REQUIRE(this->arguments_.settlementType == Settlement::Physical,
                   "cash-settled swaptions not yet implemented ...");

        Date settlement = this->model_->termStructure()->referenceDate();

        if (this->arguments_.exercise->dates().back() <=
            settlement) { // swaption is expired, possibly generated swap is not
                          // valued
            this->results_.value = 0.0;
            return;
        }

        int idx = static_cast<int>(this->arguments_.exercise->dates().size()) - 1;
        int minIdxAlive = static_cast<int>(
            std::upper_bound(this->arguments_.exercise->dates().begin(),
                             this->arguments_.exercise->dates().end(), settlement) -
            this->arguments_.exercise->dates().begin());

        VanillaSwap swap = *this->arguments_.swap;
        Option::Type type =
            this->arguments_.type == VanillaSwap::Payer ? Option::Call : Option::Put;
        Schedule fixedSchedule = swap.fixedSchedule();
        Schedule floatSchedule = swap.floatingSchedule();

        Array npv0(2 * integrationPoints_ + 1, 0.0),
            npv1(2 * integrationPoints_ + 1, 0.0);
        Array z = this->model_->yGrid(stddevs_, integrationPoints_);
        Array p(z.size(), 0.0);

        // for probability computation
        std::vector<Array> npvp0, npvp1;
        if (probabilities_ != None) {
            for (Size i = 0; i < static_cast<Size>(idx - minIdxAlive + 2); ++i) {
                Array npvTmp0(2 * integrationPoints_ + 1, 0.0);
                Array npvTmp1(2 * integrationPoints_ + 1, 0.0);
                npvp0.push_back(npvTmp0);
                npvp1.push_back(npvTmp1);
            }
        }
        // end probabkility computation

        Date expiry1 = Null<Date>(), expiry0;
        Time expiry1Time = Null<Real>(), expiry0Time;

        do {

            if (idx == minIdxAlive - 1)
                expiry0 = settlement;
            else
                expiry0 = this->arguments_.exercise->dates()[idx];

            expiry0Time = std::max(
                this->model_->termStructure()->timeFromReference(expiry0), 0.0);

            Size j1 =
                std::upper_bound(fixedSchedule.dates().begin(),
                                 fixedSchedule.dates().end(), expiry0 - 1) -
                fixedSchedule.dates().begin();
            Size k1 =
                std::upper_bound(floatSchedule.dates().begin(),
                                 floatSchedule.dates().end(), expiry0 - 1) -
                floatSchedule.dates().begin();

            // a lazy object is not thread safe, neither is the caching
            // in gsrprocess. therefore we trigger computations here such
            // that neither lazy object recalculation nor write access
            // during caching occurs in the parallized loop below.
            // this is known to work for the gsr and markov functional
            // model implementations of Gaussian1dModel
#ifdef _OPENMP
            if (expiry1Time != Null<Real>())
                this->model_->yGrid(stddevs_, integrationPoints_, expiry1Time,
                              expiry0Time, 0.0);
            if (expiry0 > settlement) {
                for (Size l = k1; l < this->arguments_.floatingCoupons.size(); l++) {
                    this->model_->forwardRate(this->arguments_.floatingFixingDates[l],
                                        expiry0, 0.0,
                                        this->arguments_.swap->iborIndex());
                    this->model_->zerobond(this->arguments_.floatingPayDates[l], expiry0,
                                     0.0, discountCurve_);
                }
                for (Size l = j1; l < this->arguments_.fixedCoupons.size(); l++) {
                    this->model_->zerobond(this->arguments_.fixedPayDates[l], expiry0, 0.0,
                                     discountCurve_);
                }
                this->model_->numeraire(expiry0Time, 0.0, discountCurve_);
            }
#endif

#pragma omp parallel for default(shared) firstprivate(p) if(expiry0>settlement)
            for (Size k = 0; k < (expiry0 > settlement ? npv0.size() : 1);
                 k++) {

                Real price = 0.0;
                if (expiry1Time != Null<Real>()) {
                    Array yg = this->model_->yGrid(stddevs_, integrationPoints_,
                                             expiry1Time, expiry0Time,
                                             expiry0 > settlement ? z[k] : 0.0);
                    CubicInterpolation payoff0(
                        z.begin(), z.end(), npv1.begin(),
                        CubicInterpolation::Spline, true,
                        CubicInterpolation::Lagrange, 0.0,
                        CubicInterpolation::Lagrange, 0.0);
                    for (Size i = 0; i < yg.size(); i++) {
                        p[i] = payoff0(yg[i], true);
                    }
                    CubicInterpolation payoff1(
                        z.begin(), z.end(), p.begin(),
                        CubicInterpolation::Spline, true,
                        CubicInterpolation::Lagrange, 0.0,
                        CubicInterpolation::Lagrange, 0.0);
                    for (Size i = 0; i < z.size() - 1; i++) {
                        price += this->model_->gaussianShiftedPolynomialIntegral(
                            0.0, payoff1.cCoefficients()[i],
                            payoff1.bCoefficients()[i],
                            payoff1.aCoefficients()[i], p[i], z[i], z[i],
                            z[i + 1]);
                    }
                    if (extrapolatePayoff_) {
                        if (flatPayoffExtrapolation_) {
                            price += this->model_->gaussianShiftedPolynomialIntegral(
                                0.0, 0.0, 0.0, 0.0, p[z.size() - 2],
                                z[z.size() - 2], z[z.size() - 1], 100.0);
                            price += this->model_->gaussianShiftedPolynomialIntegral(
                                0.0, 0.0, 0.0, 0.0, p[0], z[0], -100.0, z[0]);
                        } else {
                            if (type == Option::Call)
                                price +=
                                    this->model_->gaussianShiftedPolynomialIntegral(
                                        0.0,
                                        payoff1.cCoefficients()[z.size() - 2],
                                        payoff1.bCoefficients()[z.size() - 2],
                                        payoff1.aCoefficients()[z.size() - 2],
                                        p[z.size() - 2], z[z.size() - 2],
                                        z[z.size() - 1], 100.0);
                            if (type == Option::Put)
                                price +=
                                    this->model_->gaussianShiftedPolynomialIntegral(
                                        0.0, payoff1.cCoefficients()[0],
                                        payoff1.bCoefficients()[0],
                                        payoff1.aCoefficients()[0], p[0], z[0],
                                        -100.0, z[0]);
                        }
                    }
                }

                npv0[k] = price;

                // for probability computation
                if (probabilities_ != None) {
                    for (Size m = 0; m < npvp0.size(); m++) {
                        Real price = 0.0;
                        if (expiry1Time != Null<Real>()) {
                            Array yg = this->model_->yGrid(
                                stddevs_, integrationPoints_, expiry1Time,
                                expiry0Time, expiry0 > settlement ? z[k] : 0.0);
                            CubicInterpolation payoff0(
                                z.begin(), z.end(), npvp1[m].begin(),
                                CubicInterpolation::Spline, true,
                                CubicInterpolation::Lagrange, 0.0,
                                CubicInterpolation::Lagrange, 0.0);
                            for (Size i = 0; i < yg.size(); i++) {
                                p[i] = payoff0(yg[i], true);
                            }
                            CubicInterpolation payoff1(
                                z.begin(), z.end(), p.begin(),
                                CubicInterpolation::Spline, true,
                                CubicInterpolation::Lagrange, 0.0,
                                CubicInterpolation::Lagrange, 0.0);
                            for (Size i = 0; i < z.size() - 1; i++) {
                                price +=
                                    this->model_->gaussianShiftedPolynomialIntegral(
                                        0.0, payoff1.cCoefficients()[i],
                                        payoff1.bCoefficients()[i],
                                        payoff1.aCoefficients()[i], p[i], z[i],
                                        z[i], z[i + 1]);
                            }
                            if (extrapolatePayoff_) {
                                if (flatPayoffExtrapolation_) {
                                    price +=
                                        this->model_
                                            ->gaussianShiftedPolynomialIntegral(
                                                  0.0, 0.0, 0.0, 0.0,
                                                  p[z.size() - 2],
                                                  z[z.size() - 2],
                                                  z[z.size() - 1], 100.0);
                                    price +=
                                        this->model_
                                            ->gaussianShiftedPolynomialIntegral(
                                                  0.0, 0.0, 0.0, 0.0, p[0],
                                                  z[0], -100.0, z[0]);
                                } else {
                                    if (type == Option::Call)
                                        price +=
                                            this->model_
                                                ->gaussianShiftedPolynomialIntegral(
                                                      0.0,
                                                      payoff1.cCoefficients()
                                                          [z.size() - 2],
                                                      payoff1.bCoefficients()
                                                          [z.size() - 2],
                                                      payoff1.aCoefficients()
                                                          [z.size() - 2],
                                                      p[z.size() - 2],
                                                      z[z.size() - 2],
                                                      z[z.size() - 1], 100.0);
                                    if (type == Option::Put)
                                        price +=
                                            this->model_
                                                ->gaussianShiftedPolynomialIntegral(
                                                      0.0,
                                                      payoff1
                                                          .cCoefficients()[0],
                                                      payoff1
                                                          .bCoefficients()[0],
                                                      payoff1
                                                          .aCoefficients()[0],
                                                      p[0], z[0], -100.0, z[0]);
                                }
                            }
                        }

                        npvp0[m][k] = price;
                    }
                }
                // end probability computation

                if (expiry0 > settlement) {
                    Real floatingLegNpv = 0.0;
                    for (Size l = k1; l < this->arguments_.floatingCoupons.size();
                         l++) {
                        floatingLegNpv +=
                            this->arguments_.nominal *
                            this->arguments_.floatingAccrualTimes[l] *
                            (this->arguments_.floatingSpreads[l] +
                             this->model_->forwardRate(
                                 this->arguments_.floatingFixingDates[l], expiry0,
                                 z[k], this->arguments_.swap->iborIndex())) *
                            this->model_->zerobond(this->arguments_.floatingPayDates[l],
                                             expiry0, z[k], discountCurve_);
                    }
                    Real fixedLegNpv = 0.0;
                    for (Size l = j1; l < this->arguments_.fixedCoupons.size(); l++) {
                        fixedLegNpv +=
                            this->arguments_.fixedCoupons[l] *
                            this->model_->zerobond(this->arguments_.fixedPayDates[l],
                                             expiry0, z[k], discountCurve_);
                    }
                    Real exerciseValue =
                        (type == Option::Call ? 1.0 : -1.0) *
                        (floatingLegNpv - fixedLegNpv) /
                        this->model_->numeraire(expiry0Time, z[k], discountCurve_);

                    // for probability computation
                    if (probabilities_ != None) {
                        if (idx == static_cast<int>(
                                       this->arguments_.exercise->dates().size()) -
                                       1) // if true we are at the latest date,
                                          // so we init
                                          // the no call probability
                            npvp0.back()[k] =
                                probabilities_ == Naive
                                    ? 1.0
                                    : 1.0 / (this->model_->zerobond(expiry0Time, 0.0,
                                                              0.0,
                                                              discountCurve_) *
                                             this->model_->numeraire(expiry0, z[k],
                                                               discountCurve_));
                        if (exerciseValue >= npv0[k]) {
                            npvp0[idx - minIdxAlive][k] =
                                probabilities_ == Naive
                                    ? 1.0
                                    : 1.0 /
                                          (this->model_->zerobond(expiry0Time, 0.0,
                                                            0.0,
                                                            discountCurve_) *
                                           this->model_->numeraire(expiry0Time, z[k],
                                                             discountCurve_));
                            for (Size ii = idx - minIdxAlive + 1;
                                 ii < npvp0.size(); ii++)
                                npvp0[ii][k] = 0.0;
                        }
                    }
                    // end probability computation

                    npv0[k] = std::max(npv0[k], exerciseValue);
                }
            }

            npv1.swap(npv0);

            // for probability computation
            if (probabilities_ != None) {
                for (Size i = 0; i < npvp0.size(); i++)
                    npvp1[i].swap(npvp0[i]);
            }
            // end probability computation

            expiry1 = expiry0;
            expiry1Time = expiry0Time;

        } while (--idx >= minIdxAlive - 1);

        this->results_.value = npv1[0] * this->model_->numeraire(0.0, 0.0, discountCurve_);

        // for probability computation
        if (probabilities_ != None) {
            std::vector<Real> prob(npvp0.size());
            for (Size i = 0; i < npvp0.size(); i++) {
                prob[i] = npvp1[i][0] *
                          (probabilities_ == Naive
                               ? 1.0
                               : this->model_->numeraire(0.0, 0.0, discountCurve_));
            }
            this->results_.additionalResults["probabilities"] = prob;
        }
        // end probability computation
    } // performCalculations
  
} // namespace QuantLib

#endif
