/*
 Sgl template library for optimizing sparse group lasso penalized objectives.
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

#ifndef DIMCONFIG_H_
#define DIMCONFIG_H_

class DimConfig {

private:

	sgl::natural_vector index;

public:

	sgl::vector const L1_penalty_weights;

	sgl::vector const L2_penalty_weights;

	/*
	 * Dimension of each block (i.e this is a vector of length n_blocks)
	 */
	const sgl::natural_vector block_dim;

	/*
	 * Number of blocks
	 */
	const sgl::natural n_blocks;

	/*
	 * Total dimension i.e. sum block_dim
	 */
	const sgl::natural dim;

	/*
	 * Block unit dimension - all blocks has size = 0 mod block_unit_dim
	 */
	const sgl::natural block_unit_dim;

	DimConfig() :
			index(static_cast<sgl::natural>(0)), L1_penalty_weights(static_cast<sgl::natural>(0)), L2_penalty_weights(
					static_cast<sgl::natural>(0)), block_dim(static_cast<sgl::natural>(0)), n_blocks(0), dim(0), block_unit_dim(0) {
	}

	DimConfig(sgl::natural_vector const& block_dim, sgl::natural unit_dim, sgl::vector const& L1_penalty_weights, sgl::vector const& L2_penalty_weights) :
			index(block_dim.n_elem + 1), L1_penalty_weights(L1_penalty_weights), L2_penalty_weights(L2_penalty_weights), block_dim(
					block_dim), n_blocks(block_dim.n_elem), dim(sum(block_dim)), block_unit_dim(unit_dim) {

		//Domain checks
		if (!sgl::is_non_negative(L1_penalty_weights)) {
			throw std::domain_error("L1 weights are negative");
		}
		if (!sgl::is_non_negative(L2_penalty_weights)) {
			throw std::domain_error("L2 weights are negative");
		}
		if (sum(block_dim) != L1_penalty_weights.n_elem) {
			throw std::logic_error("DimConfig : Dimension mismatch");
		}
		if (block_dim.n_elem != L2_penalty_weights.n_elem) {
			throw std::logic_error("DimConfig : Dimension mismatch");
		}

		//Initialise index
		index(0) = 0;
		for (sgl::natural i = 1; i < index.n_elem; i++) {
			index(i) = index(i - 1) + block_dim(i - 1);
		}

	}

	sgl::natural block_start_index(sgl::natural block_index) const;

	sgl::natural block_end_index(sgl::natural block_index) const;

	sgl::vector L1_penalty_weight(sgl::natural block_index) const;

	sgl::numeric L2_penalty_weight(sgl::natural block_index) const;

};

inline sgl::natural DimConfig::block_start_index(sgl::natural block_index) const {
	return index(block_index);
}

inline sgl::natural DimConfig::block_end_index(sgl::natural block_index) const {
	return index(block_index + 1) - 1;
}

inline sgl::vector DimConfig::L1_penalty_weight(sgl::natural block_index) const {
	return L1_penalty_weights.subvec(block_start_index(block_index), block_end_index(block_index));
}

inline sgl::numeric DimConfig::L2_penalty_weight(sgl::natural block_index) const {
	return L2_penalty_weights(block_index);
}

//FIXME remove
//DimConfig createDimConfig(sgl::vector const& feature_weights, sgl::vector const& class_weights, bool block_zero_no_penalty = true, double feature0 = 0) {
//
//	sgl::natural number_of_features = feature_weights.n_elem + (block_zero_no_penalty ? 1 : 0);
//	sgl::natural number_of_classes = class_weights.n_elem;
//
//	//Create block_dim
//	sgl::natural_vector block_dim(number_of_features);
//	block_dim.fill(number_of_classes);
//
//	//Create L1 weights
//	sgl::vector penalty_weight_L1(zeros(number_of_features * number_of_classes));
//
//	for (sgl::natural i = block_zero_no_penalty ? 1 : 0; i < number_of_features; ++i) {
//		penalty_weight_L1.subvec(i * number_of_classes, (i + 1) * number_of_classes - 1) = class_weights;
//	}
//
//	//Create L2 weights
//	sgl::vector penalty_weight_L2(number_of_features);
//	penalty_weight_L2(0) = feature0;
//	penalty_weight_L2.subvec(block_zero_no_penalty ? 1 : 0, number_of_features - 1) = feature_weights;
//
//	return DimConfig(block_dim, number_of_classes, penalty_weight_L1, penalty_weight_L2);
//}

//FIXME remove
//template<bool add_block_zero = true>
//DimConfig createDimConfig(sgl::vector const& feature_weights, sgl::matrix const& parameter_weights) {
//
//	sgl::natural number_of_features = feature_weights.n_elem + (add_block_zero ? 1 : 0);
//	sgl::natural number_of_classes = parameter_weights.n_rows;
//
//	//Create block_dim
//	sgl::natural_vector block_dim(number_of_features);
//	block_dim.fill(number_of_classes);
//
//	//Create L1 weights
//	sgl::vector penalty_weight_L1(zeros(number_of_features * number_of_classes));
//
//	for (sgl::natural i = add_block_zero ? 0 : 1; i < number_of_features; ++i) {
//		penalty_weight_L1.subvec(i * number_of_classes, (i + 1) * number_of_classes - 1) = parameter_weights.col(
//				i - add_block_zero ? 1 : 0);
//	}
//
//	//Create L2 weights
//	sgl::vector penalty_weight_L2(number_of_features);
//
//	if (add_block_zero) {
//		penalty_weight_L2(0) = 0;
//		penalty_weight_L2.subvec(1, number_of_features - 1) = feature_weights;
//	}
//
//	else {
//		penalty_weight_L2.subvec(0, number_of_features - 1) = feature_weights;
//	}
//
//	return DimConfig(block_dim, number_of_classes, penalty_weight_L1, penalty_weight_L2);
//}

DimConfig createDimConfig(sgl::natural_vector const& block_dim, sgl::vector const& blocks_weights, sgl::matrix const& parameter_weights) {

	//Create L1 weights
	sgl::vector penalty_weight_L1 = reshape(parameter_weights, parameter_weights.n_elem, 1);

	if(block_dim.n_elem != blocks_weights.n_elem) {
		throw std::logic_error("createDimConfig : Dimension mismatch");
	}

	if(sum(block_dim) != penalty_weight_L1.n_elem) {
		throw std::logic_error("createDimConfig : Dimension mismatch");
	}

	return DimConfig(block_dim, parameter_weights.n_rows, penalty_weight_L1, blocks_weights);
}

DimConfig createDimConfig(sgl::natural number_of_features, sgl::natural number_of_classes, bool block_zero_no_penalty = true) {

	//Create block_dim
	sgl::natural_vector block_dim(number_of_features);
	block_dim.fill(number_of_classes);

	//Create L1 weights
	sgl::vector penalty_weight_L1(ones(number_of_features * number_of_classes));

	//Create L2 weights
	sgl::vector penalty_weight_L2(ones(number_of_features));

	if (block_zero_no_penalty) {
		penalty_weight_L1.subvec(0, number_of_classes - 1).zeros();
		penalty_weight_L2(0) = 0;
	}

	return DimConfig(block_dim, number_of_classes, penalty_weight_L1, penalty_weight_L2);
}

#endif /* DIMCONFIG_H_ */
