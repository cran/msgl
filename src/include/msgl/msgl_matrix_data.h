/* Routines for multinomial sparse group lasso regression.
 Intended for use with R.
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

#ifndef MSGL_MATRIX_DATA_H_
#define MSGL_MATRIX_DATA_H_

//template argument T = matrix type

template<typename T>
class GroupedMatrixData;

template<typename T>
class MatrixData {

public:

	sgl::natural const n_samples;
	T const data_matrix; // data matrix rows -> samples, cols -> features

	MatrixData() :
			n_samples(0), data_matrix(0, 0) {
	}

	MatrixData(T const& data_matrix, bool intercept) :
			n_samples(data_matrix.n_rows), data_matrix(create_data_matrix(intercept, data_matrix)) {
	}

	MatrixData(T const& data_matrix) :
			n_samples(data_matrix.n_rows), data_matrix(data_matrix) {
	}

	MatrixData(MatrixData<T> const& data) :
			n_samples(data.n_samples), data_matrix(data.data_matrix) {
		//TODO efficiency number of times copied
	}

	~MatrixData() {
	}

	//TODO remove
	//GroupedMatrixData<T> * create_new(sgl::integere_vector const& grouping) const;

	const MatrixData<T> operator()(Indices const& indices) const {
		return MatrixData<T>(indices.select_rows(data_matrix));
	}

	void set_matrix(T const& data_matrix) {
		const_cast<T&>(this->data_matrix) = data_matrix;
		const_cast<sgl::natural&>(this->n_samples) = data_matrix.n_rows;
	}

private:

	const T create_data_matrix(bool intercept, T const& data_matrix) {

		T m(data_matrix.n_rows, 1);

		if (intercept) {
			m.col(0).ones();
		} else {
			m.col(0).zeros();
		}

		return join_rows(m, data_matrix);
	}

};

template<typename T>
class GroupedMatrixData: public MatrixData<T> {

public:

	sgl::natural_vector const grouping; // grouping

	using MatrixData<T>::data_matrix;
	using MatrixData<T>::n_samples;

	GroupedMatrixData() :
			MatrixData<T>(), grouping(static_cast<sgl::natural>(0)) {
	}

	GroupedMatrixData(T const& data_matrix, sgl::natural_vector const& grouping, bool intercept) :
			MatrixData<T>(data_matrix, intercept), grouping(grouping) {

		this->validity();
	}

	GroupedMatrixData(T const& data_matrix, sgl::natural_vector const& grouping) :
			MatrixData<T>(data_matrix), grouping(grouping) {

		this->validity();

		//TODO clean up intercept add
	}

	GroupedMatrixData(GroupedMatrixData<T> const& data) :
			MatrixData<T>(data), grouping(data.grouping) {
	}
	;

	void set_grouping(sgl::natural_vector const& grouping) {
		const_cast<sgl::natural_vector&>(this->grouping) = grouping;

		this->validity();
	}

	GroupedMatrixData<T> & operator =(GroupedMatrixData<T> const& other) {

		if (this != &other) {

			set_matrix(other.data_matrix);
			set_grouping(other.grouping);
		}

		return *this;
	}

	const GroupedMatrixData<T> operator()(Indices const& indices) const {
		return GroupedMatrixData<T>(indices.select_rows(data_matrix), indices.select_indices(grouping));
	}

	GroupedMatrixData<T> * create_new(Indices const& indices) const {
		return new GroupedMatrixData<T>(indices.select_rows(data_matrix), indices.select_indices(grouping));
	}

private:

	void validity() {

		//TODO error msgs
		//TODO domain checks grouping - all groups represented from 0 to max??

		if (this->grouping.n_elem != n_samples) {
			throw std::domain_error("Dimension mismatch");
		}
	}
};

template<typename T>
class WeightedGroupedMatrixData: public MatrixData<T> {

public:

	sgl::vector const weights;
	sgl::natural_vector const grouping; // grouping

	using MatrixData<T>::data_matrix;
	using MatrixData<T>::n_samples;

	WeightedGroupedMatrixData() :
			MatrixData<T>(), weights(static_cast<sgl::natural>(0)), grouping(static_cast<sgl::natural>(0)) {
	}

	WeightedGroupedMatrixData(T const& data_matrix, sgl::natural_vector const& grouping, sgl::vector const& sample_weights, bool intercept) :
			MatrixData<T>(data_matrix, intercept), weights(sample_weights), grouping(grouping) {

		this->validity();
	}

	WeightedGroupedMatrixData(T const& data_matrix, sgl::natural_vector const& grouping, sgl::vector const& sample_weights) :
			MatrixData<T>(data_matrix), weights(sample_weights), grouping(grouping) {

		this->validity();

		//TODO clean up intercept add
	}

	WeightedGroupedMatrixData(WeightedGroupedMatrixData const& s) :
			MatrixData<T>(s.data_matrix), weights(s.weights), grouping(s.grouping) {
		this->validity();
	}

	void set_grouping(sgl::natural_vector const& grouping) {
		const_cast<sgl::natural_vector&>(this->grouping) = grouping;

		this->validity();
	}

	void set_weights(sgl::vector const& weights) {
		const_cast<sgl::vector&>(this->weights) = weights;

		this->validity();
	}

	WeightedGroupedMatrixData<T> & operator =(WeightedGroupedMatrixData<T> const& other) {

		if (this != &other) {

			this->set_matrix(other.data_matrix);
			const_cast<sgl::vector&>(this->weights) = other.weights;
			const_cast<sgl::natural_vector&>(this->grouping) = other.grouping;

			this->validity();
		}

		return *this;
	}

	const WeightedGroupedMatrixData<T> operator()(Indices const& indices) const {
		return WeightedGroupedMatrixData<T>(indices.select_rows(data_matrix), indices.select_indices(grouping),
				indices.select_indices(weights));
	}

	WeightedGroupedMatrixData<T> * create_new(Indices const& indices) const {
		return new WeightedGroupedMatrixData<T>(indices.select_rows(data_matrix), indices.select_indices(grouping),
				indices.select_indices(weights));
	}

private:

	void validity() {

		//TODO error msgs
		//TODO domain checks grouping - all groups represented from 0 to max??

		if (grouping.n_elem != n_samples || weights.n_elem != n_samples) {
			throw std::domain_error("WeightedGroupedMatrixData: dimension mismatch");
		}
	}
};

template<typename T>
class GroupFilter {

private:

	sgl::natural_vector groups;
	bool execulde;

public:

	GroupFilter(sgl::natural const& group, bool execulde = false) :
			groups(1), execulde(execulde) {
		groups(0) = group;
	}

	GroupFilter(sgl::natural_vector const groups, bool execulde = false) :
			groups(groups), execulde(execulde) {
	}

	const GroupedMatrixData<T> * create_filtered_data(GroupedMatrixData<T> const& data) const {

		sgl::natural_vector tmp(data.n_samples);
		tmp.zeros();

		for (sgl::natural i = 0; i < groups.n_elem; ++i) {
			if (execulde) {
				tmp += (data.grouping != groups(i));
			} else {
				tmp += (data.grouping == groups(i));
			}
		}

		//TODO what if tmp == zero

		return boost::shared_ptr<GroupedMatrixData<T> >(data.create_new(Indices(find(tmp))));
	}

};

#endif /* MSGL_MATRIX_DATA_H_ */
