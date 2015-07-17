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

//#define QL_NO_UBLAS_SUPPORT  

#include <adolc/adolc.h>

namespace QLFCT {

	template<class T> inline const T CondExpLt(const T& x, const T& y, const T& a, const T& b) { return x < y ? a : b; }
	template<class T> inline const T CondExpLe(const T& x, const T& y, const T& a, const T& b) { return x <= y ? a : b; }
	template<class T> inline const T CondExpGt(const T& x, const T& y, const T& a, const T& b) { return x > y ? a : b; }
	template<class T> inline const T CondExpGe(const T& x, const T& y, const T& a, const T& b) { return x >= y ? a : b; }
	template<class T> inline const T CondExpEq(const T& x, const T& y, const T& a, const T& b) { return x == y ? a : b; }
	template<class T> inline const T max(const T& x, const T& y) { return QLFCT::CondExpGt(x, y, x, y); }
	template<class T> inline const T min(const T& x, const T& y) { return QLFCT::CondExpLt(x, y, x, y); }
    template<class T> inline const T abs(const T& x) { return x > T(0.0) ? x : -x; }
    template<class T> inline const T sqrt(const T& x) { return sqrt(x); }
    template<class T> inline const T exp(const T& x) { return exp(x); }
    template<class T> inline const T log(const T& x) { return log(x); }
    template<class T> inline const T pow(const T& x, const T& y) { return pow(x,y); }
    template<class T> inline const T sin(const T& x) { return sin(x); }
    template<class T> inline const T cos(const T& x) { return cos(x); }
    template<class T> inline const T tan(const T& x) { return tan(x); }
    template<class T> inline const T sinh(const T& x) { return sinh(x); }
    template<class T> inline const T cosh(const T& x) { return cosh(x); }
    template<class T> inline const T tanh(const T& x) { return tanh(x); }
    template<class T> inline const T asin(const T& x) { return asin(x); }
    template<class T> inline const T acos(const T& x) { return acos(x); }
    template<class T> inline const T atan(const T& x) { return atan(x); }
    template<class T> inline const T erf(const T& x) { return erf(x); }

}

#endif
