/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Peter Caspers

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

/*! \file qladolc.hpp
    \brief defines for ADOLC usage
*/

#ifndef ql_adolc_hpp
#define ql_adolc_hpp

// currently errors are thrown, to be investigated
#define QL_NO_UBLAS_SUPPORT

#include <adolc/adolc.h>

namespace QLFCT {

	inline const adouble CondExpLt(const adouble& x, const adouble& y, const adouble& a, const adouble& b) { return x < y ? a : b; }
	inline const adouble CondExpLe(const adouble& x, const adouble& y, const adouble& a, const adouble& b) { return x <= y ? a : b; }
	inline const adouble CondExpGt(const adouble& x, const adouble& y, const adouble& a, const adouble& b) { return x > y ? a : b; }
	inline const adouble CondExpGe(const adouble& x, const adouble& y, const adouble& a, const adouble& b) { return x >= y ? a : b; }
	inline const adouble CondExpEq(const adouble& x, const adouble& y, const adouble& a, const adouble& b) { return x == y ? a : b; }
	inline const adouble max(const adouble& x, const adouble& y) { return QLFCT::CondExpGt(x, y, x, y); }
	inline const adouble min(const adouble& x, const adouble& y) { return QLFCT::CondExpLt(x, y, x, y); }
    inline const adouble abs(const adouble& x) { return x > adouble(0.0) ? x : -x; }

    inline const adouble sqrt(const adouble& x) { return ::sqrt(x); }
    inline const adouble exp(const adouble& x) { return ::exp(x); }
    inline const adouble log(const adouble& x) { return ::log(x); }
    inline const adouble pow(const adouble& x, const adouble& y) { return ::pow(x,y); }
    inline const adouble sin(const adouble& x) { return ::sin(x); }
    inline const adouble cos(const adouble& x) { return ::cos(x); }
    inline const adouble tan(const adouble& x) { return ::tan(x); }
    inline const adouble sinh(const adouble& x) { return ::sinh(x); }
    inline const adouble cosh(const adouble& x) { return ::cosh(x); }
    inline const adouble tanh(const adouble& x) { return ::tanh(x); }
    inline const adouble asin(const adouble& x) { return ::asin(x); }
    inline const adouble acos(const adouble& x) { return ::acos(x); }
    inline const adouble atan(const adouble& x) { return ::atan(x); }
    inline const adouble erf(const adouble& x) { return ::erf(x); }

}

#endif
