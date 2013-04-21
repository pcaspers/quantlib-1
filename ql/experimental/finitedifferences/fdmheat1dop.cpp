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

#include <ql/experimental/finitedifferences/fdmheat1dop.hpp>
#include <ql/methods/finitedifferences/operators/secondderivativeop.hpp>

#if !defined(QL_NO_UBLAS_SUPPORT)
#include <boost/numeric/ublas/matrix.hpp>
#endif

namespace QuantLib {


    FdmHeat1dOp::FdmHeat1dOp(const boost::shared_ptr<FdmMesher>& mesher,
			const Real alpha)
    : mesher_(mesher), alpha_(alpha),
	  xValues_(mesher->locations(0)),
	  one_(mesher->locations(0).size(),1.0),
      mapT_(SecondDerivativeOp(0,mesher).mult(alpha*one_)) {

    }

    void FdmHeat1dOp::setTime(Time t1, Time t2) {
    }

    Size FdmHeat1dOp::size() const {
        return 1;
    }

    Disposable<Array> FdmHeat1dOp::apply(const Array& u) const {
        return mapT_.apply(u);
    }

    Disposable<Array> FdmHeat1dOp::apply_direction(Size direction,
                                                   const Array& r) const {
        if (direction == 0)
            return mapT_.apply(r);
        QL_FAIL("direction too large");
    }

    Disposable<Array> FdmHeat1dOp::apply_mixed(const Array& r) const {
		Array s(r);
        return s;
    }

    Disposable<Array>
        FdmHeat1dOp::solve_splitting(Size direction,
                                     const Array& r, Real a) const {

        if (direction == 0) {
            return mapT_.solve_splitting(r, a, 1.0);
        }
        QL_FAIL("direction too large");
    }

    Disposable<Array>
        FdmHeat1dOp::preconditioner(const Array& r, Real dt) const {

        return solve_splitting(0, r, dt);
    }

#if !defined(QL_NO_UBLAS_SUPPORT)
    Disposable<std::vector<SparseMatrix> >
    FdmHeat1dOp::toMatrixDecomp() const {
        std::vector<SparseMatrix> retVal(1);
        retVal[0] = mapT_.toMatrix();
        return retVal;
    }
#endif

}
