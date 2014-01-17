/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Roland Lichters

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

#include <qlo/qladdindefines.hpp>
#include <qlo/credit.hpp>
#include <qlo/enumerations/factories/termstructuresfactory.hpp>

#include <ql/instruments/stock.hpp>
#include <ql/quote.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/termstructures/credit/interpolatedhazardratecurve.hpp>
#include <ql/termstructures/yield/piecewiseyieldcurve.hpp>
#include <ql/termstructures/credit/piecewisedefaultcurve.hpp>
#include <ql/math/interpolations/backwardflatinterpolation.hpp>
#include <ql/math/interpolations/loginterpolation.hpp>
#include <ql/pricingengines/credit/midpointcdsengine.hpp>
#include <ql/pricingengines/credit/isdacdsengine.hpp>
#include <ql/currencies/europe.hpp>

#include <ql/experimental/credit/riskybond.hpp>

#include <boost/algorithm/string/case_conv.hpp>

#include <ql/settings.hpp>

using boost::algorithm::to_upper_copy;

namespace QuantLibAddin {

    Issuer::Issuer(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            const boost::shared_ptr<QuantLib::DefaultProbabilityTermStructure>& dfts,
            const boost::shared_ptr<QuantLib::DefaultEventSet>& evtSet,
            bool permanent
        )
        : ObjectHandler::LibraryObject<QuantLib::Issuer>(properties, permanent) {

        std::vector<QuantLib::Issuer::key_curve_pair> curves(1, std::make_pair(
            QuantLib::NorthAmericaCorpDefaultKey(QuantLib::EURCurrency(),
                                                     QuantLib::SeniorSec, 
                                                     QuantLib::Period(),
                                                     1. // amount threshold
                                                     ),
            QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure>(dfts)
        ));
        libraryObject_ = boost::shared_ptr<QuantLib::Issuer>(new QuantLib::Issuer(curves, *evtSet));
    }

    DefaultEventSet::DefaultEventSet(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            const std::string& eventType,
            const QuantLib::Date& eventDate,
            const QuantLib::Currency& cur,
            QuantLib::Seniority sen,
            const QuantLib::Date& settlementDate,
            QuantLib::Real settledRecovery,
            bool permanent
        ) 
    : ObjectHandler::LibraryObject<QuantLib::DefaultEventSet>(properties, permanent) {
        // if no match return empty set
        libraryObject_ = boost::shared_ptr<QuantLib::DefaultEventSet>(new QuantLib::DefaultEventSet());

        // only one recovery parsed by now; bankruptcy events need the whole
        //  set or they fail to construct.
        std::map<QuantLib::Seniority, QuantLib::Real> rrs;
        rrs.insert(std::pair<QuantLib::Seniority, QuantLib::Real>(sen, settledRecovery));
        if(eventType==std::string("FailureToPayEvent")) {
            libraryObject_->insert(boost::shared_ptr<QuantLib::FailureToPayEvent> (
                new QuantLib::FailureToPayEvent(eventDate, cur, sen, 1.e7, 
                //implSettlemt, 
                settlementDate,
                rrs)));
        }else if(eventType==std::string("BankruptcyEvent")){
            libraryObject_->insert(boost::shared_ptr<QuantLib::BankruptcyEvent> (
                new QuantLib::BankruptcyEvent(eventDate, cur, sen, 
                settlementDate,
                rrs)));
        }
    }

    /* Code essentially copied from QuantLibAddIn::SimpleQuote. 
         Same considerations mentioned there apply. Here the
         Seniority of the quote is considered a fixed 
         property, not a market value.
    */
    RecoveryRateQuote::RecoveryRateQuote(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            QuantLib::Seniority sen,
            QuantLib::Real value,
            bool permanent) : Quote(properties, permanent) {
        libraryObject_ = recoveryQuote_ = boost::shared_ptr<QuantLib::RecoveryRateQuote>(
            new QuantLib::RecoveryRateQuote(value, sen));
    }


    CreditDefaultSwap::CreditDefaultSwap(
              const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
              QuantLib::Protection::Side side,
              QuantLib::Real notional,
              QuantLib::Rate upfront,
              QuantLib::Rate spread,
              const boost::shared_ptr<QuantLib::Schedule>& schedule,
              QuantLib::BusinessDayConvention paymentConvention,
              const QuantLib::DayCounter& dayCounter,
              bool settlesAccrual,
              bool paysAtDefaultTime,
              const QuantLib::Date& protectionStart,
              const QuantLib::Date& upfrontDate,
              bool rebatesAccrual,
              bool permanent)
        : Instrument(properties, permanent) {
			// dirty way to decide if this is constructed through a run only version
			if(upfrontDate == QuantLib::Null<QuantLib::Date>() && upfront == 0.) {
				libraryObject_ = boost::shared_ptr<QuantLib::CreditDefaultSwap>(
							new QuantLib::CreditDefaultSwap(side,
															notional,
															spread,
															*schedule,
															paymentConvention,
															dayCounter,
															settlesAccrual,
															paysAtDefaultTime,
															protectionStart,
															boost::shared_ptr<QuantLib::Claim>(),
															QuantLib::Actual360(true),
															rebatesAccrual));
			}else{
				libraryObject_ = boost::shared_ptr<QuantLib::CreditDefaultSwap>(
							new QuantLib::CreditDefaultSwap(side,
															notional,
															upfront,
															spread,
															*schedule,
															paymentConvention,
															dayCounter,
															settlesAccrual,
															paysAtDefaultTime,
															protectionStart,
															upfrontDate,
															boost::shared_ptr<QuantLib::Claim>(),
															QuantLib::Actual360(true),
															rebatesAccrual));
			}
    }
    
    MidPointCdsEngine::MidPointCdsEngine(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            const QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure>& defaultTS,
            QuantLib::Real recoveryRate,
            const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldTS,
            bool permanent) 
        : PricingEngine(properties, permanent) {
        libraryObject_ = boost::shared_ptr<QuantLib::PricingEngine>(new
              QuantLib::MidPointCdsEngine(defaultTS, recoveryRate, yieldTS));
    }

	IsdaCdsEngine::IsdaCdsEngine(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            const QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure>& defaultTS,
            QuantLib::Real recoveryRate,
            const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldTS,
			bool includeSttlDateFlows,
            bool permanent) 
        : PricingEngine(properties, permanent) {
        libraryObject_ = boost::shared_ptr<QuantLib::PricingEngine>(new
              QuantLib::IsdaCdsEngine(defaultTS, recoveryRate, yieldTS, includeSttlDateFlows));
    }



    SpreadCdsHelper::SpreadCdsHelper(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            const QuantLib::Handle<QuantLib::Quote>& quote,
            const QuantLib::Period& period,
            QuantLib::Natural settlementDays,
            const QuantLib::Calendar& calendar,
            QuantLib::Frequency frequency,
            QuantLib::BusinessDayConvention paymentConvention,
            QuantLib::DateGeneration::Rule rule,
            const QuantLib::DayCounter& dayCounter,
            QuantLib::Real recoveryRate,
            const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldTS,
            bool settlesAccrual,
            bool paysAtDefaultTime,
            bool includeLastDay,
            bool rebatesAccrual,
            bool useIsdaEngine,
            bool permanent) : DefaultProbabilityHelper(properties, permanent) {

        QuantLib::DayCounter lastPeriodDC;
        if(includeLastDay) 
            lastPeriodDC = QuantLib::Actual360(true);

        libraryObject_ = boost::shared_ptr<QuantLib::DefaultProbabilityHelper>(new
		       QuantLib::SpreadCdsHelper(quote,
						 period,
						 settlementDays,
						 calendar,
						 frequency,
						 paymentConvention,
						 rule,
						 dayCounter,
						 recoveryRate,
						 yieldTS,
						 settlesAccrual,
						 paysAtDefaultTime,
                         lastPeriodDC,
                         rebatesAccrual,
                         useIsdaEngine));
        // design improvement here? Use a proper factory with traits to initialize a template engine?, use a separate helper?
        if(useIsdaEngine) {
            boost::shared_ptr<QuantLib::CdsHelper> ptr2H = 
                boost::dynamic_pointer_cast<QuantLib::CdsHelper>(libraryObject_);
            ptr2H->setIsdaEngineParameters(QuantLib::IsdaCdsEngine::NumericalFix::Taylor,
                QuantLib::IsdaCdsEngine::AccrualBias::NoBias, 
                QuantLib::IsdaCdsEngine::ForwardsInCouponPeriod::Flat);
        }
    }

    UpfrontCdsHelper::UpfrontCdsHelper(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            const QuantLib::Handle<QuantLib::Quote>& quote,
            QuantLib::Rate runningSpread,
            const QuantLib::Period& period,
            QuantLib::Natural settlementDays,
            const QuantLib::Calendar& calendar,
            QuantLib::Frequency frequency,
            QuantLib::BusinessDayConvention paymentConvention,
            QuantLib::DateGeneration::Rule rule,
            const QuantLib::DayCounter& dayCounter,
            QuantLib::Real recoveryRate,
            const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldTS,
            QuantLib::Natural upfrontSettlementDays,
            bool settlesAccrual,
            bool paysAtDefaultTime,
            bool includeLastDay,
            bool rebatesAccrual,
            bool useIsdaEngine,
            bool permanent) : DefaultProbabilityHelper(properties, permanent) {


        QuantLib::DayCounter lastPeriodDC;
        if(includeLastDay) 
            lastPeriodDC = QuantLib::Actual360(true);

        libraryObject_ = boost::shared_ptr<QuantLib::DefaultProbabilityHelper>(new
		       QuantLib::UpfrontCdsHelper(quote,
                                          runningSpread,
                                          period,
                                          settlementDays,
                                          calendar,
                                          frequency,
                                          paymentConvention,
                                          rule,
                                          dayCounter,
                                          recoveryRate,
                                          yieldTS,
                                          upfrontSettlementDays,
                                          settlesAccrual,
                                          paysAtDefaultTime,
                                          lastPeriodDC,
                                          rebatesAccrual,
                                          useIsdaEngine));
    }

    HazardRateCurve::HazardRateCurve(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            const std::vector<QuantLib::Date>& dates,
            const std::vector<QuantLib::Rate>& hazardRates,
            const QuantLib::DayCounter& dayCounter,
            bool permanent) 
        : DefaultProbabilityTermStructure(properties, permanent) {
        QL_REQUIRE(!dates.empty(), "no input dates given");
        QL_REQUIRE(dates.size() == hazardRates.size(), 
                   "vector sizes differ");
        libraryObject_ = boost::shared_ptr<QuantLib::Extrapolator>(
        new QuantLib::InterpolatedHazardRateCurve<QuantLib::BackwardFlat>(
 				 dates, hazardRates, dayCounter));
    }

    PiecewiseHazardRateCurve::PiecewiseHazardRateCurve(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            const std::vector<boost::shared_ptr<QuantLib::DefaultProbabilityHelper> >& helpers,
            const QuantLib::DayCounter& dayCounter,
            const QuantLib::Calendar& calendar,
            const std::string& interpolator,
            QuantLib::Real accuracy,
            bool permanent) 
        : DefaultProbabilityTermStructure(properties, permanent) {

        if(interpolator == std::string("LINEAR")){
            libraryObject_ = boost::shared_ptr<QuantLib::Extrapolator>(new
                   QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate,
                        QuantLib::Linear>(
                            0, 
                            calendar,
                            helpers, 
                            dayCounter));
        }else if(interpolator == std::string("BACKWARDFLAT")) {
            libraryObject_ = boost::shared_ptr<QuantLib::Extrapolator>(new
                   QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate,
                        QuantLib::BackwardFlat>(
                            0, 
                            calendar,
                            helpers, 
                            dayCounter));
        }else{
            QL_FAIL("Unrecognised interpolator");
        }

        libraryObject_->enableExtrapolation();
    }

    // ptr type check here would correspond to template spez in the  
    //   subscribers factory solution.
    const std::vector<QuantLib::Date>& PiecewiseHazardRateCurve::dates() const {
        typedef QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate, QuantLib::BackwardFlat> flat_curve;
        typedef QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate, QuantLib::Linear> lin_curve;
        boost::shared_ptr<flat_curve> ptrBF =
            boost::dynamic_pointer_cast<flat_curve>(libraryObject_);
        if(ptrBF) return ptrBF->dates();
        boost::shared_ptr<lin_curve> ptrLIN =
            boost::dynamic_pointer_cast<lin_curve>(libraryObject_);
        if(ptrLIN) return ptrLIN->dates();
        QL_FAIL("Unable to cast default probability term structure.");
    }

    const std::vector<QuantLib::Real>& PiecewiseHazardRateCurve::data() const {
        typedef QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate, QuantLib::BackwardFlat> flat_curve;
        typedef QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate, QuantLib::Linear> lin_curve;
        boost::shared_ptr<flat_curve> ptrBF =
            boost::dynamic_pointer_cast<flat_curve>(libraryObject_);
        if(ptrBF) return ptrBF->data();
        boost::shared_ptr<lin_curve> ptrLIN =
            boost::dynamic_pointer_cast<lin_curve>(libraryObject_);
        if(ptrLIN) return ptrLIN->data();
        QL_FAIL("Unable to cast default probability term structure.");
        }        


    PiecewiseFlatForwardCurve::PiecewiseFlatForwardCurve(
            const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
            const QuantLib::Date& referenceDate,
            const std::vector<boost::shared_ptr<QuantLib::RateHelper> >& helpers,
            const QuantLib::DayCounter& dayCounter,
            QuantLib::Real accuracy,
            bool permanent)
        : YieldTermStructure(properties, permanent) {
        libraryObject_ = boost::shared_ptr<QuantLib::Extrapolator>(new
               QuantLib::PiecewiseYieldCurve<QuantLib::Discount,QuantLib::LogLinear>(referenceDate, helpers, dayCounter));
    }



    RiskyFixedBond::RiskyFixedBond(
        const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
        std::string name,
        QuantLib::Currency ccy,
        QuantLib::Real recoveryRate,
        QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure> defaultTS,
        const boost::shared_ptr<QuantLib::Schedule>& schedule,
        QuantLib::Real rate,
        QuantLib::DayCounter dayCounter,
        QuantLib::BusinessDayConvention paymentConvention,
        QuantLib::Real notional,
        QuantLib::Handle<QuantLib::YieldTermStructure> yieldTS,
        QuantLib::Date npvDate,
        bool permanent)
    : Instrument(properties, permanent) {

        std::vector<QuantLib::Real> notionals(1,notional);

        libraryObject_ = boost::shared_ptr<QuantLib::RiskyFixedBond>(
            new QuantLib::RiskyFixedBond(
                    name,ccy,recoveryRate,defaultTS,*schedule,rate,dayCounter,
                    paymentConvention,notionals,yieldTS, npvDate
                                       ));

    }


}
