
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

/*! \file processingiterator.hpp
    \brief Iterator mapping a unary function to an underlying sequence

    \fullpath
    ql/Utilities/%processingiterator.hpp
*/

// $Id$

#ifndef quantlib_processing_iterator_h
#define quantlib_processing_iterator_h

#include <iterator>

namespace QuantLib {

    namespace Utilities {

        //! Iterator mapping a unary function to an underlying sequence
        /*! This iterator advances an underlying iterator and returns the 
            values obtained by applying a unary function to the values such 
            iterator points to.

            This class was implemented based on Christopher Baus and Thomas
            Becker, <i>Custom Iterators for the STL</i>, included in the
            proceedings of the First Workshop on C++ Template Programming,
            Erfurt, Germany, 2000 (http://www.oonumerics.org/tmpw00/)
        */
        template <class Iterator, class UnaryFunction>
        class processing_iterator : public QL_ITERATOR<
            typename QL_ITERATOR_TRAITS<Iterator>::iterator_category,
            typename UnaryFunction::result_type,
            typename QL_ITERATOR_TRAITS<Iterator>::difference_type,
            const typename UnaryFunction::result_type*,
            const typename UnaryFunction::result_type&>
        {
          public:
            /* These typedefs are needed even though inherited from 
               QL_ITERATOR (see 14.6.2.3 of the standard).  */
            typedef typename UnaryFunction::result_type value_type;
            typedef typename QL_ITERATOR_TRAITS<Iterator>::difference_type
                difference_type;
            typedef const value_type* pointer;
            typedef const value_type& reference;
            // constructor
            processing_iterator(const Iterator&, const UnaryFunction&);
            //! \name Dereferencing
            //@{
            reference operator*()  const;
            pointer   operator->() const;
            //@}
            //! \name Random access
            //@{
            value_type operator[](int ) const;
            //@}
            //! \name Increment and decrement
            //@{
            processing_iterator& operator++();
            processing_iterator  operator++(int );
            processing_iterator& operator--();
            processing_iterator  operator--(int );
            processing_iterator& operator+=(difference_type);
            processing_iterator& operator-=(difference_type);
            processing_iterator  operator+(difference_type);
            processing_iterator  operator-(difference_type);
            //@}
            //! \name Difference
            //@{
            difference_type operator-(
                const processing_iterator<Iterator,UnaryFunction>&);
            //@}
            //! \name Comparisons
            //@{
            bool operator==(const processing_iterator<Iterator,UnaryFunction>&);
            bool operator!=(const processing_iterator<Iterator,UnaryFunction>&);
            bool operator< (const processing_iterator<Iterator,UnaryFunction>&);
            bool operator> (const processing_iterator<Iterator,UnaryFunction>&);
            bool operator<=(const processing_iterator<Iterator,UnaryFunction>&);
            bool operator>=(const processing_iterator<Iterator,UnaryFunction>&);
            //@}
          private:
            Iterator it_;
            UnaryFunction f_;
            mutable value_type x_;
        };

        //! helper function to create processing iterators
        /*! \relates processing_iterator */
        template <class Iterator, class UnaryFunction>
        processing_iterator<Iterator,UnaryFunction>
        make_processing_iterator(Iterator it, UnaryFunction p);



        // inline definitions

        template <class Iterator, class UnaryFunction>
        inline processing_iterator<Iterator,UnaryFunction>::processing_iterator(
          const Iterator& it, const UnaryFunction& f)
        : it_(it), f_(f) {}

        template <class Iterator, class UnaryFunction>
        inline processing_iterator<Iterator,UnaryFunction>&
        processing_iterator<Iterator,UnaryFunction>::operator++() {
            ++it_;
            return *this;
        }

        template <class Iterator, class UnaryFunction>
        inline processing_iterator<Iterator,UnaryFunction>
        processing_iterator<Iterator,UnaryFunction>::operator++(int ) {
            processing_iterator<Iterator,UnaryFunction> temp = *this;
            ++it_;
            return temp;
        }

        template <class Iterator, class UnaryFunction>
        inline processing_iterator<Iterator,UnaryFunction>&
        processing_iterator<Iterator,UnaryFunction>::operator--() {
            --it_;
            return *this;
        }

        template <class Iterator, class UnaryFunction>
        inline processing_iterator<Iterator,UnaryFunction>
        processing_iterator<Iterator,UnaryFunction>::operator--(int ) {
            processing_iterator<Iterator,UnaryFunction> temp = *this;
            --it_;
            return temp;
        }

        template <class Iterator, class UnaryFunction>
        inline processing_iterator<Iterator,UnaryFunction>&
        processing_iterator<Iterator,UnaryFunction>::operator+=(
          processing_iterator<Iterator,UnaryFunction>::difference_type i) {
            it_+=i;
            return *this;
        }

        template <class Iterator, class UnaryFunction>
        inline processing_iterator<Iterator,UnaryFunction>&
        processing_iterator<Iterator,UnaryFunction>::operator-=(
          processing_iterator<Iterator,UnaryFunction>::difference_type i) {
            it_-=i;
            return *this;
        }

        template <class Iterator, class UnaryFunction>
        inline processing_iterator<Iterator,UnaryFunction>::reference
        processing_iterator<Iterator,UnaryFunction>::operator*() const {
            x_ = f_(*it_);
            return x_;
        }

        template <class Iterator, class UnaryFunction>
        inline processing_iterator<Iterator,UnaryFunction>::pointer
        processing_iterator<Iterator,UnaryFunction>::operator->() const {
            x_ = f_(*it_);
            return &x_;
        }

        template <class Iterator, class UnaryFunction>
        inline processing_iterator<Iterator,UnaryFunction>::value_type
        processing_iterator<Iterator,UnaryFunction>::operator[](int i) const {
            return *(*this+i);
        }

        template <class Iterator, class UnaryFunction>
        inline processing_iterator<Iterator,UnaryFunction>
        processing_iterator<Iterator,UnaryFunction>::operator+(
          processing_iterator<Iterator,UnaryFunction>::difference_type i) {
            return processing_iterator<Iterator,UnaryFunction>(it_+i,f_);
        }

        template <class Iterator, class UnaryFunction>
        inline processing_iterator<Iterator,UnaryFunction>
        processing_iterator<Iterator,UnaryFunction>::operator-(
          processing_iterator<Iterator,UnaryFunction>::difference_type i) {
            return processing_iterator<Iterator,UnaryFunction>(it_-i,f_);
        }

        template <class Iterator, class UnaryFunction>
        inline processing_iterator<Iterator,UnaryFunction>::difference_type
        processing_iterator<Iterator,UnaryFunction>::operator-(
          const processing_iterator<Iterator,UnaryFunction>& i) {
            return (it_-i.it_);
        }

        template <class Iterator, class UnaryFunction>
        inline bool processing_iterator<Iterator,UnaryFunction>::operator==(
          const processing_iterator<Iterator,UnaryFunction>& i) {
            return (it_ == i.it_);
        }

        template <class Iterator, class UnaryFunction>
        inline bool processing_iterator<Iterator,UnaryFunction>::operator!=(
          const processing_iterator<Iterator,UnaryFunction>& i) {
            return (it_ != i.it_);
        }

        template <class Iterator, class UnaryFunction>
        inline bool processing_iterator<Iterator,UnaryFunction>::operator<(
          const processing_iterator<Iterator,UnaryFunction>& i) {
            return (it_ < i.it_);
        }

        template <class Iterator, class UnaryFunction>
        inline bool processing_iterator<Iterator,UnaryFunction>::operator>(
          const processing_iterator<Iterator,UnaryFunction>& i) {
            return (it_ > i.it_);
        }

        template <class Iterator, class UnaryFunction>
        inline bool processing_iterator<Iterator,UnaryFunction>::operator<=(
          const processing_iterator<Iterator,UnaryFunction>& i) {
            return (it_ <= i.it_);
        }

        template <class Iterator, class UnaryFunction>
        inline bool processing_iterator<Iterator,UnaryFunction>::operator>=(
          const processing_iterator<Iterator,UnaryFunction>& i) {
            return (it_ >= i.it_);
        }

        template <class Iterator, class UnaryFunction>
        inline processing_iterator<Iterator,UnaryFunction>
        make_processing_iterator(Iterator it, UnaryFunction p) {
            return processing_iterator<Iterator,UnaryFunction>(it,p);
        }

    }

}


#endif
