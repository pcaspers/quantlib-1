
/*
 * Copyright (C) 2000-2001 QuantLib Group
 *
 * This file is part of QuantLib.
 * QuantLib is a C++ open source library for financial quantitative
 * analysts and developers --- http://quantlib.org/
 *
 * QuantLib is free software and you are allowed to use, copy, modify, merge,
 * publish, distribute, and/or sell copies of it under the conditions stated
 * in the QuantLib License.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the license for more details.
 *
 * You should have received a copy of the license along with this file;
 * if not, please email quantlib-users@lists.sourceforge.net
 * The license is also available at http://quantlib.org/LICENSE.TXT
 *
 * The members of the QuantLib Group are listed in the Authors.txt file, also
 * available at http://quantlib.org/group.html
*/

/*! \file actualactualhistorical.hpp
    \brief act/act historical day count convention

    \fullpath
    ql/DayCounters/%actualactualhistorical.hpp
*/

// $Id$

#ifndef quantlib_actualactualhistorical_day_counter_h
#define quantlib_actualactualhistorical_day_counter_h

#include "ql/daycounter.hpp"

namespace QuantLib {

    namespace DayCounters {

        //! Actual/Actual day count according to ISDA
        /*! Actual/Actual day count according to ISDA. Also known as
            "Actual/Actual (Historical)".
            For more detail please refer to:
            http://www.isda.org/c_and_a/pdf/mktc1198.pdf
        */
        class ActualActualHistorical : public DayCounter {
          public:
            std::string name() const { return std::string("act/act(h)"); }
            int dayCount(const Date& d1,
                         const Date& d2) const {
                            return (d2-d1); }
            Time yearFraction(const Date& d1,
                              const Date& d2,
                              const Date& refPeriodStart = Date(),
                              const Date& refPeriodEnd = Date()) const;
        };

    }

}


#endif
