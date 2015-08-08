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

/*! \file qlcppad.hpp
    \brief defines for CppAD usage
*/

#ifndef ql_cppad_hpp
#define ql_cppad_hpp

// currently errors are thrown, to be investigated
#define QL_NO_UBLAS_SUPPORT  

#include <cppad/cppad.hpp>

namespace CppAD {

template <class T> inline
const T fmax(const  T& x, const T& y) {
    return CondExpGt(x, y, x, y);
}

template <class T> inline
const T fmin(const T& x, const T& y) {
    return CondExpLt(x, y, x, y);
}

}

using namespace CppAD;

#endif
