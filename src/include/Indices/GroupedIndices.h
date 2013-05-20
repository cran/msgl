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

#ifndef GROUPEDINDICES_H
#define	GROUPEDINDICES_H

#include "Indices.h"

class GroupedIndices : public Indices {
public:

    friend ostream & operator<<(ostream & output, const GroupedIndices & indices) {

        for (u32 i = 0; i < indices._groupIndices.n_elem; i++) {
            output << "[Group " << i << " ] " << indices._groupIndices(i) << endl;
        }
        return output;
    }

private:

    field<Indices> _groupIndices;

    static Indices collectAll(field<Indices> set) {

        Indices indices;

        for (u32 i = 0; i < set.n_elem; i++) {
            indices = indices + set(i);
        }

        return indices;
    }

public:

    GroupedIndices() : Indices(), _groupIndices(field<Indices>(0)) {
    };

    GroupedIndices(uvec indices, uvec grouping) : Indices(indices), _groupIndices(field<Indices>(max(grouping) + 1)) {

        ASSERT(indices.n_elem == grouping.n_elem, "size of indices and grouping must be equal");

        u32 numberOfGroups = max(grouping) + 1;

        field<uvec> groupedIndices(numberOfGroups);

        for (u32 i = 0; i < grouping.n_elem; i++) {
            uvec elm(1);
            elm = indices(i);
            groupedIndices(grouping(i)) = join_cols(groupedIndices(grouping(i)), elm);
        }

        for (u32 i = 0; i < numberOfGroups; i++) {
            _groupIndices(i) = Indices(groupedIndices(i));
        }
    }

    GroupedIndices(u32 from, u32 to, uvec grouping) : Indices(from, to), _groupIndices(field<Indices>(max(grouping) + 1)) {

        ASSERT(to - from + 1 == grouping.n_elem, "size of indices and grouping must be equal");

        u32 numberOfGroups = max(grouping) + 1;

        field<uvec> groupedIndices(numberOfGroups);

        for (u32 i = 0; i < grouping.n_elem; i++) {
            uvec elm(1);
            elm = from + i;
            groupedIndices(grouping(i)) = join_cols(groupedIndices(grouping(i)), elm);
        }

        for (u32 i = 0; i < numberOfGroups; i++) {
            if (groupedIndices(i).n_elem == 0) {
                _groupIndices(i) = Indices();
            } else {
                _groupIndices(i) = Indices(groupedIndices(i));
            }
        }
    }

    GroupedIndices(field<Indices> groupIndices) : Indices(collectAll(groupIndices)), _groupIndices(groupIndices) {
    };

    /**
     * Retrive random subset of indices
     *
     * @param sizes vector of length equal to the number of groups, each element giving the number of indices to pick from the given group
     * @param gen random number genrator
     * @return the subset of indices
     */
    template <typename RandomNumberGenerator > GroupedIndices randomSubset(uvec const & sizes, RandomNumberGenerator & gen) const {

        ASSERT(sizes.n_elem == numberOfGroups(), "length of sizes vector must be equal to the number of groups");

        field<Indices> subsets(sizes.n_elem);

        for (u32 i = 0; i < sizes.n_elem; i++) {
            subsets(i) = _groupIndices(i).randomSubset(sizes(i), gen);
        }

        return GroupedIndices(subsets);
    }

    template <typename RandomNumberGenerator > GroupedIndices randomSubset(vec const & fractions, RandomNumberGenerator & gen) const {

        ASSERT((max(fractions) <= 1) && (min(fractions) >= 0), "invalid fraction vector");
        ASSERT(fractions.n_elem == numberOfGroups(), "length of fractions vector must be equal to the number of groups");

        uvec sizes(fractions.n_elem);

        for (u32 i = 0; i < fractions.n_elem; i++) {
            sizes(i) = floor(fractions(i) * _groupIndices(i).size());
        }

        return randomSubset(sizes, gen);
    }

    //TODO template <typename RandomNumberGenrator> GroupedIndices  blancedRandomSubset(u32 size_of_subset, RandomNumberGenerator & gen) const {}

    template <typename RandomNumberGenerator > GroupedIndices blancedRandomSubset(double fraction, RandomNumberGenerator & gen) const {

        ASSERT(fraction <= 1 && fraction >= 0, "invalid fraction");

        vec fractions(numberOfGroups());
        fractions.fill(fraction);

        return randomSubset(fractions, gen);
    }

    template <typename RandomNumberGenerator> field<GroupedIndices> blancedRandomSubsets(double fraction, u32 number_of_subsets, RandomNumberGenerator & gen) const {
        field<GroupedIndices> subsets(number_of_subsets);

        for (u32 i = 0; i < number_of_subsets; i++) {
            subsets(i) = blancedRandomSubset(fraction, gen);
        }

        return subsets;
    }
    // with grouping sizes - rows = subsets, cols = groups

    template <typename RandomNumberGenerator > field<GroupedIndices> groupedDisjointSubsets(umat const & sizes, RandomNumberGenerator & gen) const {

        ASSERT(sizes.n_cols == numberOfGroups(), "dimension mismatch")

        field<Indices> subsets(sizes.n_rows, sizes.n_cols);
        for (u32 i = 0; i < numberOfGroups(); i++) {
            subsets.col(i) = _groupIndices(i).disjointSubsets(sizes.col(i), gen);
        }

        field<GroupedIndices> groupedIndices(sizes.n_cols);
        for (u32 i = 0; i < sizes.n_cols; i++) {
            groupedIndices(i) = GroupedIndices(subsets.row(i));
        }

        return groupedIndices;
    }

    template <typename RandomNumberGenerator > field<GroupedIndices> groupedDisjointSubsets(mat const & fractions, RandomNumberGenerator & gen) const {

        ASSERT((max(max(fractions)) <= 1) && (min(min(fractions)) >= 0), "invalid fraction matrix");
        ASSERT(fractions.n_cols == numberOfGroups(), "dimension mismatch")

        field<Indices> subsets(fractions.n_rows, fractions.n_cols);
        for (u32 i = 0; i < numberOfGroups(); i++) {
            subsets.col(i) = _groupIndices(i).disjointSubsets(fractions.col(i), gen);
        }

        field<GroupedIndices> groupedIndices(fractions.n_cols);
        for (u32 i = 0; i < fractions.n_cols; i++) {
            groupedIndices(i) = GroupedIndices(subsets.row(i));
        }

        return groupedIndices;

    }

    //TODO template <typename RandomNumberGenrator> field<GroupedIndices> groupedDisjointSubsets(uvec const & sizes, RandomNumberGenerator & gen) const;

    template <typename RandomNumberGenerator > field<GroupedIndices> groupedDisjointSubsets(vec const & fractions, RandomNumberGenerator & gen) const {

        ASSERT(max(fractions) <= 1 && min(fractions) >= 0, "invalid fraction matrix");

        field<Indices> subsets(fractions.n_elem, numberOfGroups());
        for (u32 i = 0; i < numberOfGroups(); i++) {
            subsets.col(i) = _groupIndices(i).disjointSubsets(fractions, gen);
        }

        field<GroupedIndices> groupedIndices(fractions.n_elem);
        for (u32 i = 0; i < fractions.n_elem; i++) {
            groupedIndices(i) = GroupedIndices(subsets.row(i));
        }

        return groupedIndices;

    }

    template <typename RandomNumberGenerator > field<GroupedIndices> groupedDisjointSubsets(u32 fold, RandomNumberGenerator & gen) const {
        ASSERT(fold <= size(), "fold > size of indices set.")
		ASSERT(fold <= max(groupSizes()), "fold > max group size")

        vec fractions(fold);
        fractions.fill(1 / (double) fold);

        return groupedDisjointSubsets(fractions, gen);
    }

    u32 numberOfGroups() const {
        return _groupIndices.n_elem;
    }

    uvec groupSizes() const {
    	uvec sizes(numberOfGroups());
    	for(u32 i = 0; i < numberOfGroups(); ++i) {
    		sizes(i) = _groupIndices(i).size();
    	}

    	return sizes;
    }

    uvec grouping() const {
        uvec grouping(size());

        for (u32 i = 0; i < numberOfGroups(); i++) {

            if (_groupIndices(i).isEmpty()) {
                continue;
            }

            uvec elements = _groupIndices(i).getElements();
            for (u32 j = 0; j < elements.n_elem; j++) {
                grouping(elements(j)) = i;
            }
        }

        return grouping;
    }

};

#endif	/* GROUPEDINDICES_H */

