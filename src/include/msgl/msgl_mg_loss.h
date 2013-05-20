/*
 Routines for multinomial sparse group lasso regression.
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

#ifndef MSGL_MG_LOSS_H_
#define MSGL_MG_LOSS_H_

class MultinomialSampleLoss {

public:

	const sgl::natural n_samples;
	const sgl::natural n_classes;

private:

	sgl::natural_vector const& Y;
	sgl::vector const& W; //weights - vector of length n_samples

	sgl::matrix prob; //probabilities: n_samples x n_classes

public:

	mutable sgl::matrix_field hessian_matrices;
	mutable bool hessians_computed;


	MultinomialSampleLoss() :
			n_samples(0), n_classes(0), Y(static_cast<u32>(0)), W(static_cast<u32>(0)), prob(n_samples, n_classes), hessian_matrices(static_cast<u32>(0)), hessians_computed(false) {
	}

	template<typename E>
	MultinomialSampleLoss(WeightedGroupedMatrixData<E> const& data) :
			n_samples(data.grouping.n_elem), n_classes(max(data.grouping) + 1), Y(data.grouping), W(data.weights), prob(
					n_samples, n_classes), hessian_matrices(n_samples), hessians_computed(false) {

		set_lp_zero();
	}

	//template<typename T>
	void set_lp(sgl::matrix const& link);
	void set_lp_zero();

	const sgl::matrix gradients() const {

		sgl::matrix grad = trans(prob);

		for (sgl::natural i = 0; i < n_samples; ++i) {
			grad(Y(i), i) -= 1;
			grad.col(i) *= W(i);
		}

		return grad;
	}

	void compute_hessians() const {

		if (hessians_computed) {
			return;
		}

		for(sgl::natural i = 0; i < n_samples; ++i) {
			hessian_matrices(i) = W(i) * (diagmat(prob.row(i)) - trans(prob.row(i)) * prob.row(i));
		}

		hessians_computed = true;
	}

	const sgl::matrix& hessians(u32 i) const {
		return hessian_matrices(i);
	}

	const sgl::matrix& probabilities() const {
		return prob;
	}

	const sgl::numeric sum_values() const {

		sgl::numeric val = 0;

		for (sgl::natural i = 0; i < n_samples; ++i) {
			val += -W(i) * log(prob(i, Y(i)));
		}

		return (val);

	}

};

// linp n_samples x n_classes the linear predictors
//template<typename T>
void MultinomialSampleLoss::set_lp(sgl::matrix const& linp) {

	TIMER_START;

	prob = exp(linp);

	for (sgl::natural i = 0; i < n_samples; ++i) {
		prob.row(i) *= 1 / as_scalar(sum(prob.row(i), 1));
	}

	ASSERT_IS_FINITE(prob);
	hessians_computed = false;
}

void MultinomialSampleLoss::set_lp_zero() {

	prob.fill(1 / static_cast<sgl::numeric>(n_classes));
	hessians_computed = false;
}

template<typename T, typename E>
class GenralizedLinearLossBase: public T {

public:

	typedef WeightedGroupedMatrixData<E> data_type;

	sgl::DimConfig const& dim_config;

protected:

	E const& X; //design matrix - n_features x n_samples

	sgl::natural const n_samples;
	sgl::natural const n_classes;
	sgl::natural const n_features;

	mutable sgl::matrix partial_hessian;
	mutable sgl::natural_vector hessian_diag_mat_computed;
	mutable sgl::matrix_field hessian_diag_mat;

	sgl::parameter current_parameters;

	sgl::vector x_norm;
	sgl::numeric x_norm_max;
	mutable sgl::numeric partial_hessian_norm;
	mutable sgl::numeric level0_bound;

public:

	GenralizedLinearLossBase(data_type const& data, sgl::DimConfig const& dim_config);

	~GenralizedLinearLossBase() {
	}

	//Change evaluation point

	void at(sgl::parameter const& parameters);
	void at_zero();

	//Evaluation

	sgl::numeric evaluate() const;

	//Gradient

	sgl::vector const gradient() const;
	sgl::block_vector const gradient(sgl::natural_vector const& indices) const;

	sgl::vector const compute_block_gradient(sgl::natural feature_index) const;

	sgl::numeric hessian_bound_level0() const;
	sgl::numeric hessian_bound_level1(sgl::natural block_index) const;

	static sgl::natural compute_unit_size(data_type const& data) {
		return max(data.grouping) + 1;
	}

	static sgl::natural compute_number_of_units(data_type const& data) {
		return data.data_matrix.n_cols;
	}

protected:

	void compute_hessian_norm() const;
};

template<typename T, typename E>
GenralizedLinearLossBase<T, E>::GenralizedLinearLossBase(data_type const& data, sgl::DimConfig const& dim_config) :
		T(data), dim_config(dim_config), X(data.data_matrix), n_samples(data.n_samples), n_classes(max(data.grouping) + 1), n_features(
				X.n_cols), partial_hessian(n_classes, n_samples), hessian_diag_mat_computed(dim_config.n_blocks), hessian_diag_mat(
				dim_config.n_blocks), current_parameters(dim_config), x_norm(dim_config.n_blocks) {

	TIMER_START;

	//Initialize outer and x_norm
	for (sgl::natural j = 0; j < dim_config.n_blocks; ++j) {
		x_norm(j) = as_scalar(max(sqrt(sum(square(X.cols(dim_config.block_start_index(j) / n_classes, dim_config.block_end_index(j) / n_classes)))), 1));
	}

	x_norm_max = x_norm.max();
}

template<typename T, typename E>
void GenralizedLinearLossBase<T, E>::at(const sgl::parameter & parameters) {

	TIMER_START;

	current_parameters = parameters;

	sgl::matrix lp(X * trans(parameters.as_matrix()));

	T::set_lp(lp);

	partial_hessian.zeros();
	hessian_diag_mat_computed.zeros();
	compute_hessian_norm();
}

template<typename T, typename E>
void GenralizedLinearLossBase<T, E>::at_zero() {

	current_parameters.zeros();
	T::set_lp_zero();

	partial_hessian.zeros();
	hessian_diag_mat_computed.zeros();
	compute_hessian_norm();
}

template<typename T, typename E>
sgl::vector const GenralizedLinearLossBase<T, E>::gradient() const {

	TIMER_START;

	return reshape(T::gradients() * X, n_features * n_classes, 1);
}

template<typename T, typename E>
sgl::block_vector const GenralizedLinearLossBase<T, E>::gradient(sgl::natural_vector const& indices) const {

	throw std::runtime_error("GenralizedLinearLoss::gradient : Not implemented");

	sgl::block_vector gradient(n_features, n_classes);

	return gradient;
}

//note this function may return inf
template<typename T, typename E>
sgl::numeric GenralizedLinearLossBase<T, E>::evaluate() const {
	return T::sum_values();
}

template<typename T, typename E>
inline sgl::vector const GenralizedLinearLossBase<T, E>::compute_block_gradient(sgl::natural block_index) const {

	TIMER_START;

	return reshape(
			partial_hessian
					* X.cols(dim_config.block_start_index(block_index) / n_classes, dim_config.block_end_index(block_index) / n_classes),
			dim_config.block_dim(block_index), 1);
}

//TODO rename
template<typename T, typename E>
inline void GenralizedLinearLossBase<T, E>::compute_hessian_norm() const {

	TIMER_START;

	partial_hessian_norm = sqrt(as_scalar(max(sum(square(partial_hessian), 1))));
	level0_bound = partial_hessian_norm * x_norm_max;
}

template<typename T, typename E>
inline sgl::numeric GenralizedLinearLossBase<T, E>::hessian_bound_level0() const {
	return level0_bound;
}

template<typename T, typename E>
inline sgl::numeric GenralizedLinearLossBase<T, E>::hessian_bound_level1(sgl::natural block_index) const {

	return partial_hessian_norm * x_norm(block_index);
}

template<typename T, typename E>
class GenralizedLinearLoss: public GenralizedLinearLossBase<T, E> {

public:

	typedef typename GenralizedLinearLossBase<T, E>::data_type data_type;

	using GenralizedLinearLossBase<T, E>::dim_config;

private:

	using GenralizedLinearLossBase<T, E>::X;

	using GenralizedLinearLossBase<T, E>::n_samples;
	using GenralizedLinearLossBase<T, E>::n_classes;
	using GenralizedLinearLossBase<T, E>::n_features;

	using GenralizedLinearLossBase<T, E>::partial_hessian;
	using GenralizedLinearLossBase<T, E>::hessian_diag_mat_computed;
	using GenralizedLinearLossBase<T, E>::hessian_diag_mat;
	using GenralizedLinearLossBase<T, E>::current_parameters;

	using GenralizedLinearLossBase<T, E>::x_norm;
	using GenralizedLinearLossBase<T, E>::x_norm_max;
	using GenralizedLinearLossBase<T, E>::partial_hessian_norm;
	using GenralizedLinearLossBase<T, E>::level0_bound;

public:

	GenralizedLinearLoss(data_type const& data, sgl::DimConfig const& dim_config);

	//Hessian
	sgl::matrix const hessian_diag(sgl::natural block_index) const;

	void hessian_update(sgl::natural block_index, sgl::parameter_block_vector const& z);

};

template<typename T, typename E>
GenralizedLinearLoss<T, E>::GenralizedLinearLoss(data_type const& data, sgl::DimConfig const& dim_config) :
		GenralizedLinearLossBase<T, E>(data, dim_config) {
}

template<typename T, typename E>
inline sgl::matrix const GenralizedLinearLoss<T, E>::hessian_diag(sgl::natural block_index) const {

	TIMER_START;

	if (hessian_diag_mat_computed(block_index) != 0) {
		return hessian_diag_mat(block_index);
	}

	T::compute_hessians();

	hessian_diag_mat(block_index).zeros(dim_config.block_dim(block_index), dim_config.block_dim(block_index));

	sgl::matrix tmp(X.cols(dim_config.block_start_index(block_index) / n_classes, dim_config.block_end_index(block_index) / n_classes));

	for (sgl::natural i = 0; i < n_samples; ++i) {

		sgl::matrix H = T::hessians(i);

		for (sgl::natural j = 0; j < tmp.n_cols; ++j) {
			for (sgl::natural k = j; k < tmp.n_cols; ++k) {
				hessian_diag_mat(block_index).submat(j * n_classes, k * n_classes, (j + 1) * n_classes - 1, (k + 1) * n_classes - 1) += tmp(
						i, j) * tmp(i, k) * H;

			}
		}
	}

	hessian_diag_mat(block_index) = symmatu(hessian_diag_mat(block_index));

	hessian_diag_mat_computed(block_index) = 1;

	return hessian_diag_mat(block_index);
}

template<typename T, typename E>
void GenralizedLinearLoss<T, E>::hessian_update(sgl::natural block_index, sgl::parameter_block_vector const& z) {

	TIMER_START;

	//Update
	T::compute_hessians();

	E tmp1 = X.cols(dim_config.block_start_index(block_index) / n_classes, dim_config.block_end_index(block_index) / n_classes);

	sgl::parameter_block_vector tmp2(z - current_parameters.block(block_index));
	tmp2.reshape(n_classes, dim_config.block_dim(block_index) / n_classes);

	for (sgl::natural i = 0; i < n_samples; ++i) {
		partial_hessian.col(i) += T::hessians(i) * tmp2 * trans(tmp1.row(i));
	}

	this->compute_hessian_norm();

	//Update cureent x
	current_parameters.set_block(block_index, z);

}

// Sparse matrix specializations

template<>
void GenralizedLinearLoss<MultinomialSampleLoss, sgl::sparse_matrix>::hessian_update(sgl::natural block_index,
		sgl::parameter_block_vector const& z) {

	TIMER_START;

	MultinomialSampleLoss::compute_hessians();

	sgl::matrix tmp(z - current_parameters.block(block_index));
	tmp.reshape(n_classes, dim_config.block_dim(block_index) / n_classes);

	sgl::vector tmp2(n_classes);

	for (sgl::natural i = dim_config.block_start_index(block_index) / n_classes; i < dim_config.block_end_index(block_index) / n_classes+1; ++i) {

		tmp2 = tmp.col(i-dim_config.block_start_index(block_index) / n_classes);

		for (sgl::natural j = X.col_ptrs[i]; j < X.col_ptrs[i+1]; ++j) {

			sgl::natural row = X.row_indices[j];

			partial_hessian.col(row) += MultinomialSampleLoss::hessians(row) * tmp2 * X.values[j];
		}
	}


	this->compute_hessian_norm();

	//Update current x
	current_parameters.set_block(block_index, z);
}

template<>
inline sgl::matrix const GenralizedLinearLoss<MultinomialSampleLoss, sgl::sparse_matrix>::hessian_diag(sgl::natural block_index) const {

	TIMER_START;

	if (hessian_diag_mat_computed(block_index) != 0) {
		return hessian_diag_mat(block_index);
	}

	MultinomialSampleLoss::compute_hessians();

	hessian_diag_mat(block_index).zeros(dim_config.block_dim(block_index), dim_config.block_dim(block_index));

	sgl::sparse_matrix tmp(
			X.cols(dim_config.block_start_index(block_index) / n_classes, dim_config.block_end_index(block_index) / n_classes));

	for (sgl::natural j = 0; j < tmp.n_cols; ++j) {
		for (sgl::natural k = j; k < tmp.n_cols; ++k) {
			for (sgl::natural i1 = tmp.col_ptrs[j]; i1 < tmp.col_ptrs[j + 1]; ++i1) {

				sgl::natural row1 = tmp.row_indices[i1];

				sgl::numeric vi2 = 0;

				for (sgl::natural i2 = tmp.col_ptrs[k]; i2 < tmp.col_ptrs[k + 1]; ++i2) {

					sgl::natural row2 = tmp.row_indices[i2];

					if (row1 != row2) {
						continue;
					}

					vi2 += tmp.values[i2];
				}

				if(vi2 != 0) {
					hessian_diag_mat(block_index).submat(j * n_classes, k * n_classes, (j + 1) * n_classes - 1, (k + 1) * n_classes - 1) += vi2 * tmp.values[i1] * MultinomialSampleLoss::hessians(row1);
				}
			}
		}
	}
	
	hessian_diag_mat(block_index) = symmatu(hessian_diag_mat(block_index));

	hessian_diag_mat_computed(block_index) = 1;

	return hessian_diag_mat(block_index);
}

typedef sgl::ObjectiveFunctionType<GenralizedLinearLoss<MultinomialSampleLoss, sgl::matrix> , WeightedGroupedMatrixData<sgl::matrix> > gl_weighted_multinomial;
typedef sgl::ObjectiveFunctionType<GenralizedLinearLoss<MultinomialSampleLoss, sgl::sparse_matrix>
		, WeightedGroupedMatrixData<sgl::sparse_matrix> > sp_gl_weighted_multinomial;

#endif /* MSGL_MG_LOSS_H_ */
