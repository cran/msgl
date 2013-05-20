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

#ifndef MSGL_MULTINOMIAL_WEIGHTED_SPARSE_H_
#define MSGL_MULTINOMIAL_WEIGHTED_SPARSE_H_

extern "C" {

R::SEXP r_msgl_wb_sp_basic(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_sample_weights, R::SEXP r_block_dim, R::SEXP r_blockWeights,
		R::SEXP r_parameterWeights, R::SEXP r_alpha, R::SEXP r_lambda, R::SEXP r_needed_solutions, R::SEXP r_config);

R::SEXP r_msgl_wb_sp_lambda_seq(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_sample_weights, R::SEXP r_block_dim, R::SEXP r_blockWeights,
		R::SEXP r_parameterWeights, R::SEXP r_alpha, R::SEXP r_numberOfModels, R::SEXP r_lambdaMin, R::SEXP r_config);

R::SEXP r_msgl_wb_sp_cv(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_sample_weights, R::SEXP r_block_dim, R::SEXP r_blockWeights,
		R::SEXP r_parameterWeights, R::SEXP r_alpha, R::SEXP r_lambda_seq, R::SEXP r_fold, R::SEXP r_cv_indices, R::SEXP r_use_cv_indices,
		R::SEXP r_number_of_threads, R::SEXP r_seed, R::SEXP r_config);

R::SEXP r_msgl_wb_sp_subsampling(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_sample_weights, R::SEXP r_block_dim, R::SEXP r_blockWeights,
		R::SEXP r_parameterWeights, R::SEXP r_alpha, R::SEXP r_lambda_seq, R::SEXP r_training_samples, R::SEXP r_test_samples,
		R::SEXP r_number_of_threads, R::SEXP r_config);

R::SEXP r_msgl_sparse_predict(R::SEXP r_x, R::SEXP r_beta);

}

R::SEXP msgl_wb_sp_lambda_seq(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_sample_weights, R::SEXP r_block_dim, R::SEXP r_blockWeights,
		R::SEXP r_parameterWeights, R::SEXP r_alpha, R::SEXP r_numberOfModels, R::SEXP r_lambdaMin, R::SEXP r_config) {

	TIMER_START;

	//Map data
	const sgl::sparse_matrix X = get_value < sgl::sparse_matrix > (r_x);
	const sgl::natural_vector Y = get_value < sgl::natural_vector > (r_classes);
	const sgl::vector W = get_value < sgl::vector > (r_sample_weights);
	const sgl::natural_vector block_dim = get_value < sgl::natural_vector > (r_block_dim);
	const sgl::vector blockWeights = get_value < sgl::vector > (r_blockWeights);
	const sgl::matrix parameterWeights = get_value < sgl::matrix > (r_parameterWeights);
	const sgl::numeric alpha = get_value < sgl::numeric > (r_alpha);

	// Create optimiser
	rList rlist_config(r_config);
	const msgl::AlgorithmConfiguration config(rlist_config);

	sgl::DimConfig dim_config = sgl::createDimConfig(block_dim, blockWeights, parameterWeights);
	msgl::WeightedGroupedMatrixData < sgl::sparse_matrix > data(X, Y, W, true);
	msgl::sp_gl_weighted_multinomial obj_type(data);

	sgl::Interface < msgl::AlgorithmConfiguration, msgl::sp_gl_weighted_multinomial > sgl_optimizer(obj_type, dim_config, alpha, config);

	sgl::vector result = sgl_optimizer.lambda_sequence(sgl_optimizer.lambda_max(), get_value < sgl::numeric > (r_lambdaMin),
			get_value < sgl::natural > (r_numberOfModels));

	return (rObject(result));
}

R::SEXP r_msgl_wb_sp_lambda_seq(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_sample_weights, R::SEXP r_block_dim, R::SEXP r_blockWeights,
		R::SEXP r_parameterWeights, R::SEXP r_alpha, R::SEXP r_numberOfModels, R::SEXP r_lambdaMin, R::SEXP r_config) {

	try {

		return msgl_wb_sp_lambda_seq(r_x, r_classes, r_sample_weights, r_block_dim, r_blockWeights, r_parameterWeights, r_alpha,
				r_numberOfModels, r_lambdaMin, r_config);

		//Catch unhandled exceptions

	} catch (std::exception & e) {

		if(e.what() != NULL) {
			SGL_ERROR(e.what());
		}

		else {
			SGL_ERROR("Unknown error");
		}

	} catch (...) {
		SGL_ERROR("Unknown error");
	}

	return R::R_NilValue; //Avoid compiler warnings
}

/* msgl_basic
 *
 */

R::SEXP msgl_wb_sp_basic(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_sample_weights, R::SEXP r_block_dim, R::SEXP r_blockWeights,
		R::SEXP r_parameterWeights, R::SEXP r_alpha, R::SEXP r_lambda_seq, R::SEXP r_needed_solutions, R::SEXP r_config) {

	//Start scope timer, note will only be activated if SGL_TIMING is defined
	TIMER_START;

	//Load data
	const sgl::sparse_matrix X = get_value < sgl::sparse_matrix > (r_x);
	const sgl::natural_vector Y = get_value < sgl::natural_vector > (r_classes);
	const sgl::vector W = get_value < sgl::vector > (r_sample_weights);
	const sgl::natural_vector block_dim = get_value < sgl::natural_vector > (r_block_dim);
	const sgl::vector blockWeights = get_value < sgl::vector > (r_blockWeights);
	const sgl::matrix parameterWeights = get_value < sgl::matrix > (r_parameterWeights);
	const sgl::numeric alpha = get_value < sgl::numeric > (r_alpha);
	const sgl::vector lambda_seq = get_value < sgl::vector > (r_lambda_seq);
	const sgl::natural_vector needed_solutions = get_value < sgl::natural_vector > (r_needed_solutions);

	// Create optimiser
	rList rlist_config(r_config);
	const msgl::AlgorithmConfiguration config(rlist_config);
	sgl::DimConfig dim_config = sgl::createDimConfig(block_dim, blockWeights, parameterWeights);

	if (config.verbose) {
		rout << "Msgl, sparse" << endl;
		rout << "Number of blocks : " << dim_config.n_blocks << " - total dimension : " << dim_config.dim << " - L2 penalty for block 0 : "
				<< dim_config.L2_penalty_weight(0) << endl;
	}

	msgl::WeightedGroupedMatrixData < sgl::sparse_matrix > data(X, Y, W, true);
	msgl::sp_gl_weighted_multinomial obj_type(data);

	sgl::Interface < msgl::AlgorithmConfiguration, msgl::sp_gl_weighted_multinomial > sgl_optimizer(obj_type, dim_config, alpha, config);

	// Solve problem
	sgl::block_vector_field x_field(needed_solutions.n_elem);
	sgl::vector object_value(needed_solutions.n_elem);
	sgl::vector function_value(needed_solutions.n_elem);

	sgl_optimizer.optimize(x_field, needed_solutions, object_value, function_value, lambda_seq);

	boost::shared_ptr < rList > res_ptr;

	//Build result R list
	res_ptr = boost::shared_ptr < rList > (new rList(4));
	rList & res = *res_ptr.get();

	//TODO better solution
	sgl::sparse_matrix_field beta(x_field.n_elem);
	for (sgl::natural i = 0; i < beta.n_elem; ++i) {
		beta(i) = x_field(i).as_matrix();
	}

	res.attach(rObject(beta), "beta");
	res.attach(rObject(object_value), "loss");
	res.attach(rObject(function_value), "objective");
	res.attach(r_lambda_seq, "lambda");

	return (res);
}

R::SEXP r_msgl_wb_sp_basic(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_sample_weights, R::SEXP r_block_dim, R::SEXP r_blockWeights,
		R::SEXP r_parameterWeights, R::SEXP r_alpha, R::SEXP r_lambda_seq, R::SEXP r_needed_solutions, R::SEXP r_config) {

	try {

		return msgl_wb_sp_basic(r_x, r_classes, r_sample_weights, r_block_dim, r_blockWeights, r_parameterWeights, r_alpha, r_lambda_seq,
				r_needed_solutions, r_config);

		//Catch unhandled exceptions
	} catch (std::exception & e) {

		if(e.what() != NULL) {
			SGL_ERROR(e.what());
		}

		else {
			SGL_ERROR("Unknown error");
		}

	} catch (...) {
		SGL_ERROR("Unknown error");
	}

	return R::R_NilValue; //Avoid compiler warnings
}

/* msgl_cv
 *
 */

R::SEXP msgl_wb_sp_cv(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_sample_weights, R::SEXP r_block_dim, R::SEXP r_blockWeights,
		R::SEXP r_parameterWeights, R::SEXP r_alpha, R::SEXP r_lambda_seq, R::SEXP r_fold, R::SEXP r_cv_indices, R::SEXP r_use_cv_indices,
		R::SEXP r_number_of_threads, R::SEXP r_seed, R::SEXP r_config) {

	//Map data
	const sgl::sparse_matrix X = get_value < sgl::sparse_matrix > (r_x);
	const sgl::natural_vector Y = get_value < sgl::natural_vector > (r_classes);
	const sgl::vector W = get_value < sgl::vector > (r_sample_weights);
	const sgl::natural_vector block_dim = get_value < sgl::natural_vector > (r_block_dim);
	const sgl::vector blockWeights = get_value < sgl::vector > (r_blockWeights);
	const sgl::matrix parameterWeights = get_value < sgl::matrix > (r_parameterWeights);
	const sgl::numeric alpha = get_value < sgl::numeric > (r_alpha);
	const sgl::vector lambda_seq = get_value < sgl::vector > (r_lambda_seq);
	const sgl::natural fold = get_value < sgl::natural > (r_fold);
	const sgl::natural number_of_threads = get_value < sgl::natural > (r_number_of_threads);
	const unsigned int seed = get_value<unsigned int>(r_seed);
	const bool use_cv_indices = get_value<bool>(r_use_cv_indices);

	//Configuration
	rList rlist_config(r_config);
	const msgl::AlgorithmConfiguration config(rlist_config);

	sgl::DimConfig dim_config = sgl::createDimConfig(block_dim, blockWeights, parameterWeights);

	if (config.verbose) {
		rout << "Msgl, cross validation, sparse" << endl;
		rout << "Number of blocks : " << dim_config.n_blocks << " - total dimension : " << dim_config.dim << " - L2 penalty for block 0 : "
				<< dim_config.L2_penalty_weight(0) << endl;
	}

	msgl::WeightedGroupedMatrixData < sgl::sparse_matrix > data(X, Y, W, true);
	msgl::sp_gl_weighted_multinomial obj_type(data);

	sgl::Interface < msgl::AlgorithmConfiguration, msgl::sp_gl_weighted_multinomial > sgl_optimizer(obj_type, dim_config, alpha, config);

	//Cv split
	field < Indices > cvgroups;

	GroupedIndices const indices(0, data.n_samples - 1, data.grouping);

	if (use_cv_indices) {
		cvgroups = get_field < Indices > (r_cv_indices);
	}

	else {

		boost::mt19937 gen;
		gen.seed(seed);

		cvgroups = conv < Indices, GroupedIndices > (indices.groupedDisjointSubsets(fold, gen));
	}

	msgl::MultinomialPredictor < sgl::sparse_matrix > predictor;

	//Do cv
	boost::tuple < field<msgl::MultinomialResponse>, sgl::vector, sgl::vector > response_field = sgl_optimizer.regular_cv(predictor,
			lambda_seq, cvgroups, indices, number_of_threads);

	//Build R list
	boost::tuple < sgl::matrix_field, sgl::matrix_field, sgl::natural_matrix > result = convert(response_field.get<0>());

	rList res(6);

	res.attach(rObject(result.get<0>()), "link");
	res.attach(rObject(result.get<1>()), "response");
	res.attach(rObject(result.get<2>()), "classes");

	res.attach(rObject(cvgroups), "cv.indices");

	res.attach(rObject(response_field.get<1>()), "features");
	res.attach(rObject(response_field.get<2>()), "parameters");

	return res;
}

R::SEXP r_msgl_wb_sp_cv(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_sample_weights, R::SEXP r_block_dim, R::SEXP r_blockWeights,
		R::SEXP r_parameterWeights, R::SEXP r_alpha, R::SEXP r_lambda_seq, R::SEXP r_fold, R::SEXP r_cv_indices, R::SEXP r_use_cv_indices,
		R::SEXP r_number_of_threads, R::SEXP r_seed, R::SEXP r_config) {

	try {

		return msgl_wb_sp_cv(r_x, r_classes, r_sample_weights, r_block_dim, r_blockWeights, r_parameterWeights, r_alpha, r_lambda_seq,
				r_fold, r_cv_indices, r_use_cv_indices, r_number_of_threads, r_seed, r_config);

		//Catch unhandled exceptions

	} catch (std::exception & e) {

		if(e.what() != NULL) {
			SGL_ERROR(e.what());
		}

		else {
			SGL_ERROR("Unknown error");
		}

	} catch (...) {
		SGL_ERROR("Unknown error");
	}

	return R::R_NilValue; //Avoid compiler warnings
}

boost::tuple<field<field<msgl::MultinomialResponse> >, sgl::natural_matrix, sgl::natural_matrix> subsampling(sgl::sparse_matrix const& X,
		sgl::natural_vector const& Y, sgl::vector const& W, sgl::natural_vector const& block_dim, sgl::vector const& blockWeights,
		sgl::matrix const& parameterWeights, sgl::numeric alpha, sgl::vector const& lambda_seq, sgl::natural number_of_threads,
		field<Indices> const& training_samples, field<Indices> const& test_samples,  msgl::AlgorithmConfiguration const& config) {

	sgl::DimConfig dim_config = sgl::createDimConfig(block_dim, blockWeights, parameterWeights);
	msgl::WeightedGroupedMatrixData < sgl::sparse_matrix > data(X, Y, W, true);
	msgl::sp_gl_weighted_multinomial obj_type(data);

	if (config.verbose) {
		rout << "Msgl, subsampling, sparse" << endl;
		rout << "Number of blocks : " << dim_config.n_blocks << " - total dimension : " << dim_config.dim << " - L2 penalty for block 0 : "
				<< dim_config.L2_penalty_weight(0) << endl;
	}

	sgl::Interface < msgl::AlgorithmConfiguration, msgl::sp_gl_weighted_multinomial > sgl_optimizer(obj_type, dim_config, alpha, config);

	msgl::MultinomialPredictor < sgl::sparse_matrix > predictor;

	//Do subsampling
	return sgl_optimizer.subsampling(predictor, lambda_seq, training_samples, test_samples, number_of_threads);
}

R::SEXP msgl_wb_sp_subsampling(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_sample_weights, R::SEXP r_block_dim, R::SEXP r_blockWeights,
		R::SEXP r_parameterWeights, R::SEXP r_alpha, R::SEXP r_lambda_seq, R::SEXP r_training_samples, R::SEXP r_test_samples,
		R::SEXP r_number_of_threads, R::SEXP r_config) {

	//Map data
	const sgl::sparse_matrix X = get_value < sgl::sparse_matrix > (r_x);
	const sgl::natural_vector Y = get_value < sgl::natural_vector > (r_classes);
	const sgl::vector W = get_value < sgl::vector > (r_sample_weights);
	const sgl::natural_vector block_dim = get_value < sgl::natural_vector > (r_block_dim);
	const sgl::vector blockWeights = get_value < sgl::vector > (r_blockWeights);
	const sgl::matrix parameterWeights = get_value < sgl::matrix > (r_parameterWeights);
	const sgl::numeric alpha = get_value < sgl::numeric > (r_alpha);
	const sgl::vector lambda_seq = get_value < sgl::vector > (r_lambda_seq);
	const sgl::natural number_of_threads = get_value < sgl::natural > (r_number_of_threads);
	const field<Indices> training_samples = get_field < Indices > (r_training_samples);
	const field<Indices> test_samples = get_field < Indices > (r_test_samples);

	//Configuration
	rList rlist_config(r_config);
	const msgl::AlgorithmConfiguration config(rlist_config);

	boost::tuple < field<field<msgl::MultinomialResponse> >, sgl::natural_matrix, sgl::natural_matrix > response_field = subsampling(X, Y,
			W, block_dim, blockWeights, parameterWeights, alpha, lambda_seq, number_of_threads, training_samples, test_samples, config);

	sgl::natural n_subsamples = training_samples.n_elem;

	field < sgl::matrix_field > link(n_subsamples);
	field < sgl::matrix_field > response(n_subsamples);
	field < sgl::natural_matrix > classes(n_subsamples);

	for (sgl::natural i = 0; i < n_subsamples; ++i) {
		boost::tuple < sgl::matrix_field, sgl::matrix_field, sgl::natural_matrix > result = convert(response_field.get<0>()(i));
		link(i) = result.get<0>();
		response(i) = result.get<1>();
		classes(i) = result.get<2>();
	}

	rList res(6);

	res.attach(rObject(link), "link");
	res.attach(rObject(response), "response");
	res.attach(rObject(classes), "classes");

	res.attach(rObject(response_field.get<1>()), "features");
	res.attach(rObject(response_field.get<2>()), "parameters");

	return res;

}

R::SEXP r_msgl_wb_sp_subsampling(R::SEXP r_x, R::SEXP r_classes, R::SEXP r_sample_weights, R::SEXP r_block_dim, R::SEXP r_blockWeights,
		R::SEXP r_parameterWeights, R::SEXP r_alpha, R::SEXP r_lambda_seq, R::SEXP r_training_samples, R::SEXP r_test_samples,
		R::SEXP r_number_of_threads, R::SEXP r_config) {

	try {

		return msgl_wb_sp_subsampling(r_x, r_classes, r_sample_weights, r_block_dim, r_blockWeights, r_parameterWeights, r_alpha,
				r_lambda_seq, r_training_samples, r_test_samples, r_number_of_threads, r_config);

		//Catch unhandled exceptions

	} catch (std::exception & e) {

		if(e.what() != NULL) {
			SGL_ERROR(e.what());
		}

		else {
			SGL_ERROR("Unknown error");
		}

	} catch (...) {
		SGL_ERROR("Unknown error");
	}

	return R::R_NilValue; //Avoid compiler warnings
}

/* msgl_predict_classes
 *
 */

R::SEXP msgl_sparse_predict(R::SEXP r_x, R::SEXP r_beta) {

	//TODO domain checks

	//Map data
	const sgl::sparse_matrix X = get_value < sgl::sparse_matrix > (r_x);
	const sgl::sparse_matrix_field beta = get_field < sgl::sparse_matrix > (r_beta);

	msgl::MultinomialPredictor < sgl::sparse_matrix > predictor;

	boost::tuple < sgl::matrix_field, sgl::matrix_field, sgl::natural_matrix > result = convert(
			predictor.predict(msgl::MatrixData < sgl::sparse_matrix > (X, true), beta));

	rList res(3);
	res.attach(rObject(result.get<0>()), "link");
	res.attach(rObject(result.get<1>()), "response");
	res.attach(rObject(result.get<2>()), "classes");

	return res;
}

R::SEXP r_msgl_sparse_predict(R::SEXP r_x, R::SEXP r_beta) {

	try {

		return msgl_sparse_predict(r_x, r_beta);

		//Catch unhandled exceptions
	} catch (std::exception & e) {

		if(e.what() != NULL) {
			SGL_ERROR(e.what());
		}

		else {
			SGL_ERROR("Unknown error");
		}

	} catch (...) {
		SGL_ERROR("Unknown error");
	}

	return R::R_NilValue; //Avoid compiler warnings
}

#endif /* MSGL_MULTINOMIAL_WEIGHTED_SPARSE_H_ */
