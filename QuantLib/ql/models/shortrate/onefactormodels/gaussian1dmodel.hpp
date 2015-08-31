/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2013, 2015 Peter Caspers

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

/*! \file gaussian1dmodel.hpp
    \brief basic interface for one factor interest rate models

    TODO as it turns out it is not optimal for all implementations to work
    with the normalized state variable y instead of the original x (e.g.
    in the monte carlo swaption engine or in the lgm model), revisit this
*/

// uncomment to enable NTL support (see below for more details and references)
// #define GAUSS1D_ENABLE_NTL

#ifndef quantlib_gaussian1dmodel_hpp
#define quantlib_gaussian1dmodel_hpp

#include <ql/models/model.hpp>
#include <ql/models/parameter.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/indexes/swapindex.hpp>
#include <ql/instruments/vanillaswap.hpp>
#include <ql/time/date.hpp>
#include <ql/time/period.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/utilities/null.hpp>
#include <ql/patterns/lazyobject.hpp>
#include <ql/patterns/curiouslyrecurring.hpp>
#include <ql/math/interpolations/cubicinterpolation.hpp>

#ifdef GAUSS1D_ENABLE_NTL
#include <boost/math/bindings/rr.hpp>
#endif

#if defined(__GNUC__) &&                                                       \
    (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 8)) || (__GNUC__ > 4))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif
#include <boost/math/special_functions/erf.hpp>
#include <boost/unordered_map.hpp>
#if defined(__GNUC__) &&                                                       \
    (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 8)) || (__GNUC__ > 4))
#pragma GCC diagnostic pop
#endif

namespace QuantLib {

/*! One factor interest rate model interface class
    The only methods that must be implemented by subclasses
    are the numeraire and zerobond methods for an input array
    of state variable values. The variable $y$ is understood
    to be the standardized (zero mean, unit variance) version
    of the model's original state variable $x$.

    NTL support may be enabled by defining GAUSS1D_ENABLE_NTL in this
    file. For details on NTL see
             http://www.shoup.net/ntl/

    \warning the variance of the state process conditional on
    $x(t)=x$ must be independent of the value of $x$

*/

template <class Impl>
class Gaussian1dModel : public CuriouslyRecurringTemplate<Impl>,
                        public TermStructureConsistentModel,
                        public LazyObject {
  public:
    const boost::shared_ptr<StochasticProcess1D> stateProcess() const;

    const Real numeraire(const Time t, const Real y = 0.0,
                         const Handle<YieldTermStructure> &yts =
                             Handle<YieldTermStructure>()) const;

    const Real zerobond(
        const Time T, const Time t = 0.0, const Real y = 0.0,
        const Handle<YieldTermStructure> &yts = Handle<YieldTermStructure>(),
        const bool adjusted = false) const;

    const Real numeraire(const Date &referenceDate, const Real y = 0.0,
                         const Handle<YieldTermStructure> &yts =
                             Handle<YieldTermStructure>()) const;

    const Real zerobond(
        const Date &maturity, const Date &referenceDate = Null<Date>(),
        const Real y = 0.0,
        const Handle<YieldTermStructure> &yts = Handle<YieldTermStructure>(),
        const bool adjusted = false) const;

    const Real zerobondOption(
        const Option::Type &type, const Date &expiry, const Date &valueDate,
        const Date &maturity, const Rate strike,
        const Date &referenceDate = Null<Date>(), const Real y = 0.0,
        const Handle<YieldTermStructure> &yts = Handle<YieldTermStructure>(),
        const Real yStdDevs = 7.0, const Size yGridPoints = 64,
        const bool extrapolatePayoff = true,
        const bool flatPayoffExtrapolation = false,
        const bool adjusted = false) const;

    const Real forwardRate(
        const Date &fixing, const Date &referenceDate = Null<Date>(),
        const Real y = 0.0,
        boost::shared_ptr<IborIndex> iborIdx = boost::shared_ptr<IborIndex>(),
        const bool adjusted = false) const;

    const Real swapRate(
        const Date &fixing, const Period &tenor,
        const Date &referenceDate = Null<Date>(), const Real y = 0.0,
        boost::shared_ptr<SwapIndex> swapIdx = boost::shared_ptr<SwapIndex>(),
        const bool adjusted = false) const;

    const Real swapAnnuity(
        const Date &fixing, const Period &tenor,
        const Date &referenceDate = Null<Date>(), const Real y = 0.0,
        boost::shared_ptr<SwapIndex> swapIdx = boost::shared_ptr<SwapIndex>(),
        const bool adjusted = false) const;

    /*! Computes the integral
    \f[ {2\pi}^{-0.5} \int_{a}^{b} p(x) \exp{-0.5*x*x} \mathrm{d}x \f]
    with
    \f[ p(x) = ax^4+bx^3+cx^2+dx+e \f].
    */
    const static Real gaussianPolynomialIntegral(const Real a, const Real b,
                                                 const Real c, const Real d,
                                                 const Real e, const Real x0,
                                                 const Real x1);

    /*! Computes the integral
    \f[ {2\pi}^{-0.5} \int_{a}^{b} p(x) \exp{-0.5*x*x} \mathrm{d}x \f]
    with
    \f[ p(x) = a(x-h)^4+b(x-h)^3+c(x-h)^2+d(x-h)+e \f].
    */
    const static Real
    gaussianShiftedPolynomialIntegral(const Real a, const Real b, const Real c,
                                      const Real d, const Real e, const Real h,
                                      const Real x0, const Real x1);

    /*! Generates a grid of values for the standardized state variable $y$
       at time $T$
        conditional on $y(t)=y$, covering yStdDevs standard deviations
       consisting of
        2*gridPoints+1 points */

    const Disposable<Array> yGrid(const Real yStdDevs, const int gridPoints,
                                  const Real T = 1.0, const Real t = 0,
                                  const Real y = 0) const;

    /*! Computes the standardized model state from the original one
        We use that the standard deviation is independent of $x$ here ! */
    const Real y(const Real x, const Time t) {
        return (x - stateProcess_->expectation(0.0, 0.0, t)) /
               stateProcess_->stdDeviation(0.0, 0.0, t);
    }

    const Real numeraireImpl(const Time t, const Real y,
                             const Handle<YieldTermStructure> &yts) const {
        QL_FAIL("no numeraire implementation given");
    }

    const Real zerobondImpl(const Time T, const Time t, const Real y,
                            const Handle<YieldTermStructure> &yts,
                            const bool adjusted) const {
        QL_FAIL("no zerobond implementation given");
    }


  private:
    // It is of great importance for performance reasons to cache underlying
    // swaps generated from indexes. In addition the indexes may only be given
    // as templates for the conventions with the tenor replaced by the actual
    // one later on.

    struct CachedSwapKey {
        const boost::shared_ptr<SwapIndex> index;
        const Date fixing;
        const Period tenor;
        const bool operator==(const CachedSwapKey &o) const {
            return index->name() == o.index->name() && fixing == o.fixing &&
                   tenor == o.tenor;
        }
    };

    struct CachedSwapKeyHasher
        : std::unary_function<CachedSwapKey, std::size_t> {
        std::size_t operator()(CachedSwapKey const &x) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, x.index->name());
            boost::hash_combine(seed, x.fixing.serialNumber());
            boost::hash_combine(seed, x.tenor.length());
            boost::hash_combine(seed, x.tenor.units());
            return seed;
        }
    };

    typedef boost::unordered_map<CachedSwapKey, boost::shared_ptr<VanillaSwap>,
                                 CachedSwapKeyHasher> CacheType;

    mutable CacheType swapCache_;

  protected:
    // we let derived classes register with the termstructure
    Gaussian1dModel(const Handle<YieldTermStructure> &yieldTermStructure)
        : TermStructureConsistentModel(yieldTermStructure) {
        registerWith(Settings::instance().evaluationDate());
    }

    virtual ~Gaussian1dModel() {}

    void performCalculations() const {
        evaluationDate_ = Settings::instance().evaluationDate();
        enforcesTodaysHistoricFixings_ =
            Settings::instance().enforcesTodaysHistoricFixings();
    }

    void generateArguments() {
        calculate();
        notifyObservers();
    }

    // retrieve underlying swap from cache if possible, otherwise
    // create it and store it in the cache
    boost::shared_ptr<VanillaSwap>
    underlyingSwap(const boost::shared_ptr<SwapIndex> &index,
                   const Date &expiry, const Period &tenor) const {

        CachedSwapKey k = {index, expiry, tenor};
        typename CacheType::iterator i = swapCache_.find(k);
        if (i == swapCache_.end()) {
            boost::shared_ptr<VanillaSwap> underlying =
                index->clone(tenor)->underlyingSwap(expiry);
            swapCache_.insert(std::make_pair(k, underlying));
            return underlying;
        }
        return i->second;
    }

    boost::shared_ptr<StochasticProcess1D> stateProcess_;
    mutable Date evaluationDate_;
    mutable bool enforcesTodaysHistoricFixings_;
};

template <class Impl>
inline const boost::shared_ptr<StochasticProcess1D>
Gaussian1dModel<Impl>::stateProcess() const {

    QL_REQUIRE(stateProcess_ != NULL, "state process not set");
    return stateProcess_;
}

template <class Impl>
inline const Real
Gaussian1dModel<Impl>::numeraire(const Time t, const Real y,
                           const Handle<YieldTermStructure> &yts) const {

    return this->impl().numeraireImpl(t, y, yts);
}

template <class Impl>
inline const Real
Gaussian1dModel<Impl>::zerobond(const Time T, const Time t, const Real y,
                          const Handle<YieldTermStructure> &yts,
                          const bool adjusted) const {
    return this->impl().zerobondImpl(T, t, y, yts, adjusted);
}

template <class Impl>
inline const Real
Gaussian1dModel<Impl>::numeraire(const Date &referenceDate, const Real y,
                           const Handle<YieldTermStructure> &yts) const {

    return numeraire(termStructure()->timeFromReference(referenceDate), y, yts);
}

template <class Impl>
inline const Real
Gaussian1dModel<Impl>::zerobond(const Date &maturity, const Date &referenceDate,
                          const Real y, const Handle<YieldTermStructure> &yts,
                          const bool adjusted) const {

    return zerobond(termStructure()->timeFromReference(maturity),
                    referenceDate != Null<Date>()
                        ? termStructure()->timeFromReference(referenceDate)
                        : 0.0,
                    y, yts, adjusted);
}

template <class Impl>
const Real Gaussian1dModel<Impl>::forwardRate(const Date &fixing,
                                        const Date &referenceDate, const Real y,
                                        boost::shared_ptr<IborIndex> iborIdx,
                                        const bool adjusted) const {

    QL_REQUIRE(iborIdx != NULL, "no ibor index given");

    calculate();

    if (fixing <= (evaluationDate_ + (enforcesTodaysHistoricFixings_ ? 0 : -1)))
        return iborIdx->fixing(fixing);

    Handle<YieldTermStructure> yts =
        iborIdx->forwardingTermStructure(); // might be empty, then use
                                            // model curve

    Date valueDate = iborIdx->valueDate(fixing);
    Date endDate = iborIdx->fixingCalendar().advance(
        valueDate, iborIdx->tenor(), iborIdx->businessDayConvention(),
        iborIdx->endOfMonth());
    // FIXME Here we should use the calculation date calendar ?
    Real dcf = iborIdx->dayCounter().yearFraction(valueDate, endDate);

    return (zerobond(valueDate, referenceDate, y, yts, adjusted) -
            zerobond(endDate, referenceDate, y, yts, adjusted)) /
           (dcf * zerobond(endDate, referenceDate, y, yts, adjusted));
}

template <class Impl>
const Real Gaussian1dModel<Impl>::swapRate(const Date &fixing, const Period &tenor,
                                     const Date &referenceDate, const Real y,
                                     boost::shared_ptr<SwapIndex> swapIdx,
                                     const bool adjusted) const {

    QL_REQUIRE(swapIdx != NULL, "no swap index given");

    calculate();

    if (fixing <= (evaluationDate_ + (enforcesTodaysHistoricFixings_ ? 0 : -1)))
        return swapIdx->fixing(fixing);

    Handle<YieldTermStructure> ytsf =
        swapIdx->iborIndex()->forwardingTermStructure();
    Handle<YieldTermStructure> ytsd =
        swapIdx->discountingTermStructure(); // either might be empty, then
                                             // use model curve

    Schedule sched, floatSched;

    boost::shared_ptr<VanillaSwap> underlying =
        underlyingSwap(swapIdx, fixing, tenor);

    sched = underlying->fixedSchedule();

    boost::shared_ptr<OvernightIndexedSwapIndex> oisIdx =
        boost::dynamic_pointer_cast<OvernightIndexedSwapIndex>(swapIdx);
    if (oisIdx != NULL) {
        floatSched = sched;
    } else {
        floatSched = underlying->floatingSchedule();
    }

    Real annuity = swapAnnuity(fixing, tenor, referenceDate, y, swapIdx,
                               adjusted); // should be fine for
                                          // overnightindexed swap indices as
                                          // well
    Rate floatleg = 0.0;
    if (ytsf.empty() && ytsd.empty()) { // simple 100-formula can be used
                                        // only in one curve setup
        floatleg =
            (zerobond(sched.dates().front(), referenceDate, y,
                      Handle<YieldTermStructure>(), adjusted) -
             zerobond(sched.calendar().adjust(sched.dates().back(),
                                              underlying->paymentConvention()),
                      referenceDate, y, Handle<YieldTermStructure>(),
                      adjusted));
    } else {
        for (Size i = 1; i < floatSched.size(); i++) {
            floatleg +=
                (zerobond(floatSched[i - 1], referenceDate, y, ytsf, adjusted) /
                     zerobond(floatSched[i], referenceDate, y, ytsf, adjusted) -
                 1.0) *
                zerobond(floatSched.calendar().adjust(
                             floatSched[i], underlying->paymentConvention()),
                         referenceDate, y, ytsd, adjusted);
        }
    }
    return floatleg / annuity;
}

template <class Impl>
const Real Gaussian1dModel<Impl>::swapAnnuity(const Date &fixing, const Period &tenor,
                                        const Date &referenceDate, const Real y,
                                        boost::shared_ptr<SwapIndex> swapIdx,
                                        const bool adjusted) const {

    QL_REQUIRE(swapIdx != NULL, "no swap index given");

    calculate();

    Handle<YieldTermStructure> ytsd =
        swapIdx->discountingTermStructure(); // might be empty, then use
                                             // model curve

    boost::shared_ptr<VanillaSwap> underlying =
        underlyingSwap(swapIdx, fixing, tenor);

    Schedule sched = underlying->fixedSchedule();

    Real annuity = 0.0;
    for (unsigned int j = 1; j < sched.size(); j++) {
        annuity += zerobond(sched.calendar().adjust(
                                sched.date(j), underlying->paymentConvention()),
                            referenceDate, y, ytsd, adjusted) *
                   swapIdx->dayCounter().yearFraction(sched.date(j - 1),
                                                      sched.date(j));
    }
    return annuity;
}

template <class Impl>
const Real Gaussian1dModel<Impl>::zerobondOption(
    const Option::Type &type, const Date &expiry, const Date &valueDate,
    const Date &maturity, const Rate strike, const Date &referenceDate,
    const Real y, const Handle<YieldTermStructure> &yts, const Real yStdDevs,
    const Size yGridPoints, const bool extrapolatePayoff,
    const bool flatPayoffExtrapolation, const bool adjusted) const {

    calculate();

    Time fixingTime = termStructure()->timeFromReference(expiry);
    Time referenceTime =
        referenceDate == Null<Date>()
            ? 0.0
            : termStructure()->timeFromReference(referenceDate);

    Array yg = yGrid(yStdDevs, yGridPoints, fixingTime, referenceTime, y);
    Array z = yGrid(yStdDevs, yGridPoints);

    Array p(yg.size());

    for (Size i = 0; i < yg.size(); i++) {
        Real expValDsc = zerobond(valueDate, expiry, yg[i], yts, adjusted);
        Real discount =
            zerobond(maturity, expiry, yg[i], yts, adjusted) / expValDsc;
        p[i] =
            std::max((type == Option::Call ? 1.0 : -1.0) * (discount - strike),
                     0.0) /
            numeraire(fixingTime, yg[i], yts) * expValDsc;
    }

    CubicInterpolation payoff(
        z.begin(), z.end(), p.begin(), CubicInterpolation::Spline, true,
        CubicInterpolation::Lagrange, 0.0, CubicInterpolation::Lagrange, 0.0);

    Real price = 0.0;
    for (Size i = 0; i < z.size() - 1; i++) {
        price += gaussianShiftedPolynomialIntegral(
            0.0, payoff.cCoefficients()[i], payoff.bCoefficients()[i],
            payoff.aCoefficients()[i], p[i], z[i], z[i], z[i + 1]);
    }
    if (extrapolatePayoff) {
        if (flatPayoffExtrapolation) {
            price += gaussianShiftedPolynomialIntegral(
                0.0, 0.0, 0.0, 0.0, p[z.size() - 2], z[z.size() - 2],
                z[z.size() - 1], 100.0);
            price += gaussianShiftedPolynomialIntegral(0.0, 0.0, 0.0, 0.0, p[0],
                                                       z[0], -100.0, z[0]);
        } else {
            if (type == Option::Call)
                price += gaussianShiftedPolynomialIntegral(
                    0.0, payoff.cCoefficients()[z.size() - 2],
                    payoff.bCoefficients()[z.size() - 2],
                    payoff.aCoefficients()[z.size() - 2], p[z.size() - 2],
                    z[z.size() - 2], z[z.size() - 1], 100.0);
            if (type == Option::Put)
                price += gaussianShiftedPolynomialIntegral(
                    0.0, payoff.cCoefficients()[0], payoff.bCoefficients()[0],
                    payoff.aCoefficients()[0], p[0], z[0], -100.0, z[0]);
        }
    }

    return numeraire(referenceTime, y, yts) * price;
}

template <class Impl>
const Real Gaussian1dModel<Impl>::gaussianPolynomialIntegral(
    const Real a, const Real b, const Real c, const Real d, const Real e,
    const Real y0, const Real y1) {

#ifdef GAUSS1D_ENABLE_NTL
    const boost::math::ntl::RR aa = 4.0 * a, ba = 2.0 * M_SQRT2 * b,
                               ca = 2.0 * c, da = M_SQRT2 * d;
    const boost::math::ntl::RR x0 = y0 * M_SQRT1_2, x1 = y1 * M_SQRT1_2;
    const boost::math::ntl::RR res =
        (0.125 * (3.0 * aa + 2.0 * ca + 4.0 * e) * boost::math::erf(x1) -
         1.0 / (4.0 * M_SQRTPI) * exp(-x1 * x1) *
             (2.0 * aa * x1 * x1 * x1 + 3.0 * aa * x1 +
              2.0 * ba * (x1 * x1 + 1.0) + 2.0 * ca * x1 + 2.0 * da)) -
        (0.125 * (3.0 * aa + 2.0 * ca + 4.0 * e) * boost::math::erf(x0) -
         1.0 / (4.0 * M_SQRTPI) * exp(-x0 * x0) *
             (2.0 * aa * x0 * x0 * x0 + 3.0 * aa * x0 +
              2.0 * ba * (x0 * x0 + 1.0) + 2.0 * ca * x0 + 2.0 * da));
    return NTL::to_double(res.value());
#else
    const Real aa = 4.0 * a, ba = 2.0 * M_SQRT2 * b, ca = 2.0 * c,
               da = M_SQRT2 * d;
    const Real x0 = y0 * M_SQRT1_2, x1 = y1 * M_SQRT1_2;
    return (0.125 * (3.0 * aa + 2.0 * ca + 4.0 * e) * boost::math::erf(x1) -
            1.0 / (4.0 * M_SQRTPI) * exp(-x1 * x1) *
                (2.0 * aa * x1 * x1 * x1 + 3.0 * aa * x1 +
                 2.0 * ba * (x1 * x1 + 1.0) + 2.0 * ca * x1 + 2.0 * da)) -
           (0.125 * (3.0 * aa + 2.0 * ca + 4.0 * e) * boost::math::erf(x0) -
            1.0 / (4.0 * M_SQRTPI) * exp(-x0 * x0) *
                (2.0 * aa * x0 * x0 * x0 + 3.0 * aa * x0 +
                 2.0 * ba * (x0 * x0 + 1.0) + 2.0 * ca * x0 + 2.0 * da));
#endif
}

template <class Impl>
const Real Gaussian1dModel<Impl>::gaussianShiftedPolynomialIntegral(
    const Real a, const Real b, const Real c, const Real d, const Real e,
    const Real h, const Real x0, const Real x1) {
    return gaussianPolynomialIntegral(
        a, -4.0 * a * h + b, 6.0 * a * h * h - 3.0 * b * h + c,
        -4 * a * h * h * h + 3.0 * b * h * h - 2.0 * c * h + d,
        a * h * h * h * h - b * h * h * h + c * h * h - d * h + e, x0, x1);
}

template <class Impl>
const Disposable<Array> Gaussian1dModel<Impl>::yGrid(const Real stdDevs,
                                               const int gridPoints,
                                               const Real T, const Real t,
                                               const Real y) const {

    // we use that the standard deviation is independent of $x$ here !

    QL_REQUIRE(stateProcess_ != NULL, "state process not set");

    Array result(2 * gridPoints + 1, 0.0);

    Real x_t, e_0_t, e_t_T, stdDev_0_t, stdDev_t_T;
    Real stdDev_0_T = stateProcess_->stdDeviation(0.0, 0.0, T);
    Real e_0_T = stateProcess_->expectation(0.0, 0.0, T);

    if (t < QL_EPSILON) {
        // stdDev_0_t = 0.0;
        stdDev_t_T = stdDev_0_T;
        // e_0_t = 0.0;
        // x_t = 0.0;
        e_t_T = e_0_T;
    } else {
        stdDev_0_t = stateProcess_->stdDeviation(0.0, 0.0, t);
        stdDev_t_T = stateProcess_->stdDeviation(t, 0.0, T - t);
        e_0_t = stateProcess_->expectation(0.0, 0.0, t);
        x_t = y * stdDev_0_t + e_0_t;
        e_t_T = stateProcess_->expectation(t, x_t, T - t);
    }

    Real h = stdDevs / ((Real)gridPoints);

    for (int j = -gridPoints; j <= gridPoints; j++) {
        result[j + gridPoints] =
            (e_t_T + stdDev_t_T * ((Real)j) * h - e_0_T) / stdDev_0_T;
    }

    return result;
}

} // namespace QuantLib

#endif
