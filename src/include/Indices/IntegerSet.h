/*
	Indices tools.
    Copyright (C) 2012 Martin Vincent

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
*/


#ifndef INTEGERSET_H
#define	INTEGERSET_H

#include "boost/random.hpp"

#include <string>
#include <iostream>
using namespace std;

#include "../armadillo.hpp"
using namespace arma;

template <typename E>
class IntegerSetExpression {
public:

    s32 lower_bound() const {
        return static_cast<E const &> (*this).lower_bound();
    }

    s32 upper_bound() const {
        return static_cast<E const &> (*this).upper_bound();
    }

    bool operator [](s32 element) const {
        return static_cast<E const&> (*this)[element];
    }

    operator E&() {
        return static_cast<E&> (*this);
    }

    operator E const&() const {
        return static_cast<const E&> (*this);
    }

};

class IntegerSet : public IntegerSetExpression<IntegerSet> {

    friend ostream & operator<<(ostream & output, const IntegerSet & set) {
        if (set.is_empty()) {
            output << "empty" << endl;
        } else {
            output << trans(set.getElements());
        }
        return output;
    }

private:

    /* _elements is a vector of length upper_bound-lower_bound + 1
     * a non zero value at an index means that the element represented by the index is not in the set.
     */
    s32 _lower_bound;
    s32 _upper_bound;
    uvec _elements;

    template <class RandomNumberGenerator>
    static ivec _randomSubset(ivec set, u32 subset_length, RandomNumberGenerator & gen);

    template <typename E>
    void objectSet(IntegerSetExpression<E> const & set) {
        E const & s = set;
        const s32 n = s.upper_bound() - s.lower_bound() + 1;

        if (n <= 0) {
            _lower_bound = 0;
            _upper_bound = 0;
            _elements.zeros(1);
            return;
        }

        _elements.zeros(n);
        _lower_bound = s.lower_bound();
        _upper_bound = s.upper_bound();

        for (s32 i = 0; i != n; i++) {
            _elements(i) = set[i + _lower_bound] ? 1 : 0;
        }
    }

public:

    IntegerSet(const ivec & elements) : _lower_bound(min(elements)), _upper_bound(max(elements)), _elements(zeros<uvec>(_upper_bound - _lower_bound + 1)) {

        for (u32 i = 0; i != elements.n_elem; i++) {
            _elements(elements(i) - _lower_bound) = 1;
        }

    }

    IntegerSet(const uvec & elements) : _lower_bound(min(elements)), _upper_bound(max(elements)), _elements(zeros<uvec>(_upper_bound - _lower_bound + 1)) {

        for (u32 i = 0; i != elements.n_elem; i++) {
            _elements(elements(i) - _lower_bound) = 1;
        }
    }

    /*
     * The empty set
     */
    IntegerSet() : _lower_bound(0), _upper_bound(0), _elements(zeros<uvec>(1)) {
    }

    static IntegerSet EmptySet() {
        return IntegerSet();
    }

    /*
     * The set {lowe, ..., upper}
     */
    IntegerSet(s32 lower, s32 upper) : _lower_bound(lower), _upper_bound(upper), _elements(ones<uvec>(upper - lower + 1)) {
        ASSERT(lower < upper, "Lower >= upper.");
    }

    /*
     * The set {0, .., size-1}
     */
    IntegerSet(u32 size) : _lower_bound(0), _upper_bound(size - 1), _elements(ones<uvec>(size)) {
    }

    /*
     * The singleton set {element}
     */
    static IntegerSet Singleton(s32 the_element) {

        ivec elm(1);
        elm(0) = the_element;

        return IntegerSet(elm);
    }

    template <typename E>
    IntegerSet(IntegerSetExpression<E> const & set) {
        objectSet<E > (set);
    }

    //Assignment operator

    template <typename E>
            IntegerSet & operator=(IntegerSetExpression<E> const & set) {
        IntegerSet & this_set = *this;
        this_set = IntegerSet(set);
        return (this_set);
    }

    IntegerSet & operator=(IntegerSet const & set) {
        _lower_bound = set.lower_bound();
        _upper_bound = set.upper_bound();
        _elements = set._elements;

        return (* this);
    }

    s32 lower_bound() const {
        return _lower_bound;
    }

    s32 upper_bound() const {
        return _upper_bound;
    }

    bool is_empty() const {
        return size() == 0;
    }

    u32 size() const {
        return accu(_elements);
    }

    u32 pos(s32 element) const {
        return element - _lower_bound;
    }

    ivec getElements() const {

        if (is_empty()) {
            return zeros<ivec > (0);
        }

        ivec elem = zeros<ivec > (size());

        u32 j = 0;
        for (s32 i = _lower_bound; i <= _upper_bound; i++) {
            if ((*this)[i]) {
                elem(j) = i;
                j++;
            }
        }

        return elem;
    }

    bool operator [](s32 element) const {
        return element >= _lower_bound && element <= _upper_bound && _elements(pos(element)) != 0;
    }

    template <class RandomNumberGenerator>
    IntegerSet randomSubset(u32 size, RandomNumberGenerator & gen) const {
        if (this->size() == 0) {
            return IntegerSet();
        } else {
            return IntegerSet(_randomSubset(getElements(), size, gen));
        }
    }

    /*
     * DisjointSubsets
     */

    template <class RandomNumberGenerator>
    field<IntegerSet> disjointSubsets(const uvec & sizes, RandomNumberGenerator & gen) const;

    template <class RandomNumberGenerator>
    field<IntegerSet> disjointSubsets(u32 fold, RandomNumberGenerator & gen) const {

        ASSERT(fold <= size(), "Fold > size of set.")

        uvec sizes = zeros<uvec > (fold);

        sizes.fill(this->size() / fold);

        for (u32 i = 0; i < this->size() % fold; i++) {
            sizes(i)++;
        }

        return disjointSubsets(sizes, gen);
    }

};

template <class RandomNumberGenerator>
ivec IntegerSet::_randomSubset(ivec set, u32 subset_length, RandomNumberGenerator & gen) {

    ASSERT(set.n_elem >= subset_length, "_randomSubset - To few elements.");

    if (set.n_elem == 0) {
        return set;
    }

    boost::uniform_int<> uni(0, set.n_elem - 1);
    boost::variate_generator<RandomNumberGenerator&, boost::uniform_int<> > die(gen, uni);

    for (u32 i = 0; i < subset_length; i++) {
        set.swap_rows(i, die() % (set.n_elem - i) + i);
    }

    return set.rows(0, subset_length - 1);
}

template <typename E1, typename E2>
class SetUnion : public IntegerSetExpression<SetUnion<E1, E2> > {
    E1 const & _u;
    E2 const & _v;

public:

    SetUnion(IntegerSetExpression<E1> const& u, IntegerSetExpression<E2> const& v) : _u(u), _v(v) {
    }

    s32 lower_bound() const {
        return _u.lower_bound() < _v.lower_bound() ? _u.lower_bound() : _v.lower_bound();
    }

    s32 upper_bound() const {
        return _u.upper_bound() > _v.upper_bound() ? _u.upper_bound() : _v.upper_bound();
    }

    bool operator [](s32 element) const {
        return _u[element] || _v[element];
    }

};

template <typename E1, typename E2>
class SetDifference : public IntegerSetExpression<SetDifference<E1, E2> > {
    E1 const & _u;
    E2 const & _v;

public:

    SetDifference(IntegerSetExpression<E1> const& u, IntegerSetExpression<E2> const& v) : _u(u), _v(v) {
    }

    s32 lower_bound() const {
        return _u.lower_bound();
    }

    s32 upper_bound() const {
        return _u.upper_bound();
    }

    bool operator [](s32 element) const {
        return _u[element] && (!_v[element]);
    }

};

template <typename E1, typename E2>
class SetIntersection : public IntegerSetExpression<SetIntersection<E1, E2> > {
    E1 const & _u;
    E2 const & _v;

public:

    SetIntersection(IntegerSetExpression<E1> const& u, IntegerSetExpression<E2> const& v) : _u(u), _v(v) {
    }

    s32 lower_bound() const {
        return _u.lower_bound() > _v.lower_bound() ? _u.lower_bound() : _v.lower_bound();
    }

    s32 upper_bound() const {
        return _u.upper_bound() < _v.upper_bound() ? _u.upper_bound() : _v.upper_bound();
    }

    bool operator [](s32 element) const {
        return _u[element] && _v[element];
    }

};

/*
 * Overloaded operators
 */

template <typename E1, typename E2>
SetUnion<E1, E2> const operator+(IntegerSetExpression<E1> const& u, IntegerSetExpression<E2> const& v) {
    return SetUnion<E1, E2 > (u, v);
}

template <typename E1, typename E2>
SetDifference<E1, E2> const operator-(IntegerSetExpression<E1> const& u, IntegerSetExpression<E2> const& v) {
    return SetDifference<E1, E2 > (u, v);
}

template <typename E1, typename E2>
SetIntersection<E1, E2> const operator*(IntegerSetExpression<E1> const& u, IntegerSetExpression<E2> const& v) {
    return SetIntersection<E1, E2 > (u, v);
}

template <typename E1, typename E2>
bool const operator<=(IntegerSetExpression<E1> const& u, IntegerSetExpression<E2> const& v) {
    return IntegerSet(u - v).is_empty();
}

template <typename E1, typename E2>
bool const operator>=(IntegerSetExpression<E1> const& u, IntegerSetExpression<E2> const& v) {
    return v <= u;
}

template <typename E1, typename E2>
bool const operator==(IntegerSetExpression<E1> const& u, IntegerSetExpression<E2> const& v) {
    return v <= u && u <= v;
}

template <typename E1, typename E2>
bool const operator!=(IntegerSetExpression<E1> const& u, IntegerSetExpression<E2> const& v) {
    return !(v <= u) || !(u <= v);
}

template <class RandomNumberGenerator>
field<IntegerSet> IntegerSet::disjointSubsets(const uvec & sizes, RandomNumberGenerator & gen) const {

    ASSERT(accu(sizes) <= size(), "disjointSubsets - To few elements");

    field<IntegerSet> subsets(sizes.n_elem);
    IntegerSet remaining = *this;

    for (u32 i = 0; i != sizes.n_elem; i++) {
        subsets(i) = remaining.randomSubset(sizes(i), gen);
        remaining = (remaining - subsets(i));
    }

    return subsets;
}


#endif	/* INTEGERSET_H */

