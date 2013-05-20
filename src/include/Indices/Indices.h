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

#ifndef INDICES_H
#define	INDICES_H

#include <string>
#include <iostream>
using namespace std;

#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

#include <armadillo>
using namespace arma;

#include <exception>
#include "IntegerSet.h"
#include "SubMatrix.h"

class Indices
{

  // Operators
public:

  friend ostream &
  operator<<(ostream & output, const Indices & indices)
  {

    output << indices._indices;
    return output;
  }

  friend Indices
  operator-(const Indices & a, const Indices & b)
  {
    return Indices(a._indices - b._indices);
  }

  friend Indices
  operator+(const Indices & a, const Indices & b)
  {
    return Indices(a._indices + b._indices);
  }

  friend Indices
  operator*(const Indices & a, const Indices & b)
  {
    return Indices(a._indices * b._indices);
  }

private:
  IntegerSet _indices;

public:

  Indices() :
      _indices()
  {
  }

  Indices(u32 from, u32 to) :
      _indices(from, to)
  {
  }

  Indices(u32 size) :
      _indices(size)
  {
  }

  Indices(uvec const& elements) :
      _indices(elements)
  {
    //TODO assert that every element of indices is unique
  }

  //Prevent automatic conversion
  explicit
  Indices(IntegerSet const& indices) :
      _indices(indices)
  {
    //TODO assert that every element of indices is unique
  }

  //Copy constructor
  Indices(Indices const& indices) :
      _indices(indices.getElements())
  {

  }


  template<typename T>
    SpMat<T> const
    select_rows(SpMat<T> const& matrix_object) const
    {
      //TODO empty check

      uvec elements = getElements();

      /*TODO is it somewhat unsatisfactory that we have to convert the matrix to a dense representation before subsetting
       * an efficient solution using sparse representation would be nice
       */
      arma::Mat<T> m(matrix_object);
      SpMat<T> sub_object(m.rows(elements));

      return sub_object;
    }

  template<typename T>
    Mat<T> const
    select_rows(Mat<T> const& matrix_object) const
    {

      uvec elements = getElements();

      return matrix_object.rows(elements);
    }

  template<typename T>
    field<T> const
    select_rows(field<T> const& matrix_object) const
    {

      uvec elements = getElements();

      field<T> sub_object(elements.n_elem, matrix_object.n_cols);

      for (u32 j = 0; j < matrix_object.n_cols; ++j)
        {
          for (u32 i = 0; i < elements.n_elem; ++i)
            {
              sub_object(i, j) = matrix_object(elements(i), j);
            }
        }

      return sub_object;
    }

  template<typename T>
    Col<typename T::elem_type> const
    select_indices(T const& v) const
    {
      //TODO empty check

      uvec elements = getElements();

      return v.elem(elements);

    }

   template<typename T> SubMatrixRow<T> select_rows_view(T & M) const {
       //TODO empty check
        return SubMatrixRow<T>(M, getElements());
    }

  Indices const
  select(Indices const& indices) const
  {
    return (*this) * indices;
  }

  /**
   * Retrive random subset of indices
   *
   * @param size size of subset
   * @param gen random number genrator
   * @return the subset of indices
   */
  template<typename RandomNumberGenerator>
    Indices
    randomSubset(u32 size_of_subset, RandomNumberGenerator & gen) const
    {

      ASSERT(size() >= size_of_subset, "_randomSubset - To few elements.");

      Indices subset(_indices.randomSubset(size_of_subset, gen));
      return subset;
    }

  template<typename RandomNumberGenerator>
    field<Indices>
    randomSubsets(u32 size_of_subset, u32 number_of_subsets,
        RandomNumberGenerator & gen) const
    {
      field<Indices> subsets(number_of_subsets);

      for (u32 i = 0; i < number_of_subsets; i++)
        {
          subsets(i) = randomSubset(size_of_subset, gen);
        }

      return subsets;
    }

  template<typename RandomNumberGenerator>
    Indices
    randomSubset(double fraction, u32 number_of_subsets,
        RandomNumberGenerator & gen) const
    {
      return randomSubsets((u32) floor(fraction * size()), number_of_subsets,
          gen);
    }

  template<typename RandomNumberGenerator>
    Indices
    randomSubset(double fraction, RandomNumberGenerator & gen) const
    {
      return randomSubset<RandomNumberGenerator>((u32) floor(fraction * size()),
          gen);
    }

  template<typename RandomNumberGenerator>
    field<Indices>
    disjointSubsets(uvec const & sizes, RandomNumberGenerator & gen) const
    {

      ASSERT(sum(sizes) <= size(),
          "sum of sizes must be less or equal to the number of indices");

      field<IntegerSet> sets = _indices.disjointSubsets(sizes, gen);
      field<Indices> subsets(sizes.n_elem);

      for (u32 i = 0; i < sizes.n_elem; i++)
        {
          subsets(i) = Indices(sets(i));
        }

      return (subsets);
    }

  template<typename RandomNumberGenerator>
    field<Indices>
    disjointSubsets(vec const & fractions, RandomNumberGenerator & gen) const
    {

      ASSERT(sum(fractions) <= 1 + std::numeric_limits<double>::infinity(),
          "sum of fractions must be less or equal to 1");

      uvec sizes(fractions.n_elem);

      vec sorted_fractions = sort(fractions, 1);

      for (u32 i = 0; i < sizes.n_elem; i++)
        {
          sizes(i) = floor(sorted_fractions(i) * size());
        }

      u32 r = size() - sum(sizes);

      for (u32 i = 0; i < r; i++)
        {
          sizes(i)++;
        }

      return disjointSubsets(sizes, gen);
    }

  template<typename RandomNumberGenerator>
    field<Indices>
    disjointSubsets(u32 fold, RandomNumberGenerator & gen) const
    {

      ASSERT(fold <= size(), "Fold > size of indices set.");

      vec fractions(fold);
      fractions.fill(1 / (double) fold);

      return disjointSubsets(fractions, gen);
    }

  u32
  size() const
  {
    return _indices.size();
  }

  uvec
  getElements() const
  {

    if (_indices.is_empty())
      {
        return NULL;
      }

    return conv_to<uvec>::from(_indices.getElements());
  }

  bool
  isEmpty() const
  {
    return _indices.is_empty();
  }
};

#endif	/* INDICES_H */

