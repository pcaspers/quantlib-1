/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2012 Peter Caspers

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

/*! \file fdmzabrop.hpp
    \brief Zabr linear pricing operator
*/

#ifndef quantlib_fdm_zabr_op_hpp
#define quantlib_fdm_zabr_op_hpp

#include <ql/experimental/finitedifferences/firstderivativeop.hpp>
#include <ql/experimental/finitedifferences/triplebandlinearop.hpp>
#include <ql/experimental/finitedifferences/ninepointlinearop.hpp>
#include <ql/experimental/finitedifferences/fdmlinearopcomposite.hpp>

namespace QuantLib {

	class FdmZabrUnderlyingPart {
      public:
        FdmZabrUnderlyingPart(const boost::shared_ptr<FdmMesher>& mesher, const Real beta, const Real nu, const Real rho, const Real gamma);

        void setTime(Time t1, Time t2);
        const TripleBandLinearOp& getMap() const;

      protected:
        const Array volatilityValues_;
		const Array forwardValues_;
        TripleBandLinearOp mapT_;

        const boost::shared_ptr<FdmMesher> mesher_;
    };

	class FdmZabrVolatilityPart {
      public:
        FdmZabrVolatilityPart(
            const boost::shared_ptr<FdmMesher>& mesher,
			const Real beta, const Real nu, const Real rho, const Real gamma);

        void setTime(Time t1, Time t2);
        const TripleBandLinearOp& getMap() const;

      protected:
  	    const Array volatilityValues_;
		const Array forwardValues_;
        TripleBandLinearOp mapT_;

        const boost::shared_ptr<FdmMesher> mesher_;
	};


    class FdmZabrOp : public FdmLinearOpComposite {
      public:
        FdmZabrOp(const boost::shared_ptr<FdmMesher>& mesher,
			const Real beta,
			const Real nu,
			const Real rho,
			const Real gamma = 1.0); // gamma=1.0 recovers the classic sabr model

        Size size() const;
        void setTime(Time t1, Time t2);

        Disposable<Array> apply(const Array& r) const;
        Disposable<Array> apply_mixed(const Array& r) const;

        Disposable<Array> apply_direction(Size direction,
                                          const Array& r) const;
        Disposable<Array> solve_splitting(Size direction,
                                          const Array& r, Real s) const;
        Disposable<Array> preconditioner(const Array& r, Real s) const;

      private:
	    const Array volatilityValues_;
		const Array forwardValues_;
        NinePointLinearOp dxyMap_;
        FdmZabrUnderlyingPart dxMap_;
        FdmZabrVolatilityPart dyMap_;
    };
}

#endif
