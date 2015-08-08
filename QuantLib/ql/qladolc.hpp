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
	template<class T> inline T CondExpLt(const T& x, const T& y, const T& a, const T& b) { return x < y ? a : b; }
	template<> inline adouble CondExpLt(const adouble& x, const adouble& y, const adouble& a, const adouble& b) { adouble tmp; condassign(tmp, y - x , a , b); return tmp; }
	template<class T> inline T CondExpGt(const T& x, const T& y, const T& a, const T& b) { return x > y ? a : b; }
	template<> inline adouble CondExpGt(const adouble& x, const adouble& y, const adouble& a, const adouble& b) { adouble tmp; condassign(tmp, x - y , a , b); return tmp; }
	template<class T> inline T CondExpEq(const T& x, const T& y, const T& a, const T& b) { return x == y ? a : b; }
	template<> inline adouble CondExpEq(const adouble& x, const adouble& y, const adouble& a, const adouble& b) { adouble tmp1, tmp2; condassign(tmp1, x - y , b , a); condassign(tmp2, y - x, b, tmp1); return tmp2; }

}

using namespace QLFCT;
#endif
