/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2005, 2006, 2007, 2008 StatPro Italia srl
 Copyright (C) 2007, 2008 Ferdinando Ametrano
 Copyright (C) 2007 Chris Kenyon

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

/*! \file piecewiseyieldcurve.hpp
    \brief piecewise-interpolated term structure
*/

#ifndef quantlib_piecewise_yield_curve_hpp
#define quantlib_piecewise_yield_curve_hpp

#include <ql/termstructures/iterativebootstrap.hpp>
#include <ql/termstructures/yield/bootstraptraits.hpp>
#include <ql/patterns/lazyobject.hpp>
#include <ql/quote.hpp>

namespace QuantLib {

    //! Piecewise yield term structure
    /*! This term structure is bootstrapped on a number of interest
        rate instruments which are passed as a vector of handles to
        RateHelper instances. Their maturities mark the boundaries of
        the interpolated segments.

        Each segment is determined sequentially starting from the
        earliest period to the latest and is chosen so that the
        instrument whose maturity marks the end of such segment is
        correctly repriced on the curve.

        \warning The bootstrapping algorithm will raise an exception if
                 any two instruments have the same maturity date.

        \ingroup yieldtermstructures

        \test
        - the correctness of the returned values is tested by
          checking them against the original inputs.
        - the observability of the term structure is tested.
    */
    template <class Traits, class Interpolator,
              template <class> class Bootstrap = IterativeBootstrap>
    class PiecewiseYieldCurve
        : public Traits::template curve<Interpolator>::type,
          public LazyObject {
      private:
        typedef typename Traits::template curve<Interpolator>::type base_curve;
      public:
        typedef Traits traits_type;
        typedef Interpolator interpolator_type;
        //! \name Constructors
        //@{
        PiecewiseYieldCurve(
               const Date& referenceDate,
               const std::vector<boost::shared_ptr<typename Traits::helper> >&
                                                                  instruments,
               const DayCounter& dayCounter,
               const Handle<Quote>& turnOfYearEffect = Handle<Quote>(),
               Real accuracy = 1.0e-12,
               const Interpolator& i = Interpolator())
        : base_curve(referenceDate, dayCounter, i),
          instruments_(instruments),
          turnOfYearEffect_(turnOfYearEffect), accuracy_(accuracy) {
            setTurnOfYear();
            registerWith(turnOfYearEffect_);
            bootstrap_.setup(this);
        }
        PiecewiseYieldCurve(
               Natural settlementDays,
               const Calendar& calendar,
               const std::vector<boost::shared_ptr<typename Traits::helper> >&
                                                                  instruments,
               const DayCounter& dayCounter,
               const Handle<Quote>& turnOfYearEffect = Handle<Quote>(),
               Real accuracy = 1.0e-12,
               const Interpolator& i = Interpolator())
        : base_curve(settlementDays, calendar, dayCounter, i),
          instruments_(instruments),
          turnOfYearEffect_(turnOfYearEffect), accuracy_(accuracy) {
            setTurnOfYear();
            registerWith(turnOfYearEffect_);
            bootstrap_.setup(this);
        }
        //@}
        //! \name YieldTermStructure interface
        //@{
        Date maxDate() const;
        //@}
        //! \name Inspectors
        //@{
        const std::vector<Time>& times() const;
        const std::vector<Date>& dates() const;
        const std::vector<Real>& data() const;
        std::vector<std::pair<Date, Real> > nodes() const;
        //@}
        //! \name Observer interface
        //@{
        void update();
        //@}
      private:
        // methods
        void performCalculations() const;
        DiscountFactor discountImpl(Time) const;
        void setTurnOfYear();
        // data members
        std::vector<boost::shared_ptr<typename Traits::helper> > instruments_;
        Handle<Quote> turnOfYearEffect_;
        Real accuracy_;
        Date latestReference_;
        Time turnOfYear_;
        // bootstrapper classes are declared as friend to manipulate
        // the curve data. They might be passed the data instead, but
        // it would increase the complexity---which is quite high
        // enough already.
        friend class Bootstrap<PiecewiseYieldCurve<Traits,
                                                   Interpolator,
                                                   Bootstrap> >;
        friend class BootstrapError<PiecewiseYieldCurve<Traits,
                                                        Interpolator,
                                                        Bootstrap> >;
        Bootstrap<PiecewiseYieldCurve<Traits,
                                      Interpolator,
                                      Bootstrap> > bootstrap_;
    };


    // inline definitions

    template <class C, class I, template <class> class B>
    inline Date PiecewiseYieldCurve<C,I,B>::maxDate() const {
        calculate();
        return this->dates_.back();
    }

    template <class C, class I, template <class> class B>
    inline const std::vector<Time>& PiecewiseYieldCurve<C,I,B>::times() const {
        calculate();
        return this->times_;
    }

    template <class C, class I, template <class> class B>
    inline const std::vector<Date>& PiecewiseYieldCurve<C,I,B>::dates() const {
        calculate();
        return this->dates_;
    }

    template <class C, class I, template <class> class B>
    inline const std::vector<Real>& PiecewiseYieldCurve<C,I,B>::data() const {
        calculate();
        return this->data_;
    }

    template <class C, class I, template <class> class B>
    inline std::vector<std::pair<Date, Real> >
    PiecewiseYieldCurve<C,I,B>::nodes() const {
        calculate();
        return base_curve::nodes();
    }

    template <class C, class I, template <class> class B>
    inline void PiecewiseYieldCurve<C,I,B>::update() {
        base_curve::update();
        LazyObject::update();
        if (base_curve::referenceDate() != latestReference_)
            setTurnOfYear();
    }

    template <class C, class I, template <class> class B>
    inline
    DiscountFactor PiecewiseYieldCurve<C,I,B>::discountImpl(Time t) const {
        calculate();

        if ((!turnOfYearEffect_.empty()) && t>turnOfYear_) {
            QL_REQUIRE(turnOfYearEffect_->isValid(),
                       "invalid turnOfYearEffect quote");
            DiscountFactor turnOfYearEffect = turnOfYearEffect_->value();
            QL_REQUIRE(turnOfYearEffect > 0.0 && turnOfYearEffect <= 1.0,
                       "invalid turnOfYearEffect value: " << turnOfYearEffect);
            return turnOfYearEffect * base_curve::discountImpl(t);
        }

        return base_curve::discountImpl(t);
    }

    template <class C, class I, template <class> class B>
    inline void PiecewiseYieldCurve<C,I,B>::setTurnOfYear() {
        Date referenceDate = base_curve::referenceDate();
        Date turnOfYear = Date(31, December, referenceDate.year());
        turnOfYear_ = base_curve::timeFromReference(turnOfYear);
        latestReference_ = referenceDate;
    }


    // template definitions

    template <class C, class I, template <class> class B>
    void PiecewiseYieldCurve<C,I,B>::performCalculations() const {
        // just delegate to the bootstrapper
        bootstrap_.calculate();
    }

}

#endif
