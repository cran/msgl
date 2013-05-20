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
#ifndef INTERFACE_BASIC_H_
#define INTERFACE_BASIC_H_

template<typename CONFIG, typename ObjectiveFunctionType>
class Interface {

public:

	const sgl::numeric alpha;

	sgl::SglProblem<CONFIG> const sgl;
	sgl::SglOptimizer<sgl::SglProblem<CONFIG> > const optimizer;
	ObjectiveFunctionType const& objective_type;

	/**
	 * Construct a Interface.
	 *
	 * @param objective_type objective function type
	 * @param dim_config dimension configuration for the sgl optimizer
	 * @param alpha alpha values of the sgl optimizer
	 * @param config algorithm configuration
	 */
	Interface(ObjectiveFunctionType const& objective_type,
			sgl::DimConfig const& dim_config, sgl::numeric alpha,
			CONFIG const& config) :
			alpha(alpha), sgl(dim_config, config), optimizer(sgl, alpha), objective_type(
					objective_type) {

		//TODO obj_func -> check dim match with dim_config
	}

	/**
	 * Optimize the given objective.
	 * Returns results for all lambda values.
	 *
	 * @param lambda_sequence lambda sequence to optimize over
	 * @return respectively fitted parameter values, objective function values and penalised objective function_value.
	 */
	boost::tuple<sgl::block_vector_field, sgl::vector, sgl::vector>
	optimize(sgl::vector const& lambda_sequence) const;

	/**
	 * Optimize the given objective.
	 * Returns results only for the lambda values with indices given in needed_solutions.
	 * This reduces memory use.
	 *
	 * @param lambda_sequence the lambda sequence to optimize over
	 * @param needed_solutions indices of the lambda values for the needed solutions
	 * @return respectively fitted parameter values, objective function values and penalised objective function_value.
	 */
	boost::tuple<sgl::block_vector_field, sgl::vector, sgl::vector>
	optimize(sgl::vector const& lambda_sequence,
			sgl::natural_vector needed_solutions) const;

	/**
	 * Optimize the given objective.
	 *
	 * @param x_field will after the call returns contain the fitted parameters
	 * @param needed_solutions indices of the lambda values for the needed solutions
	 * @param object_value will after the call returns contain the value of the objective at the give lambda values
	 * @param function_value will after the call returns contain the value of the penalty + objective at the give lambda values
	 * @param lambda_sequence the lambda sequence to optimize over
	 */
	sgl::natural
	optimize(sgl::parameter_field & x_field,
			sgl::natural_vector needed_solutions, sgl::vector & object_value,
			sgl::vector & function_value,
			const sgl::vector & lambda_sequence) const;

	/**
	 *
	 * @param predictor
	 * @param lambda
	 * @param cv_indices
	 * @param indices_all
	 * @param number_of_threads
	 * @param do_refit
	 * @return
	 */
	template<typename Predictor>
	boost::tuple<field<typename Predictor::response_type>, sgl::vector,
			sgl::vector>
	regular_cv(Predictor const& predictor, sgl::vector const& lambda,
			field<Indices> const& cv_indices, Indices const& indices_all,
			sgl::natural number_of_threads) const;

	template<typename Predictor>
	boost::tuple<field<field<typename Predictor::response_type> >,
			sgl::natural_matrix, sgl::natural_matrix>
	subsampling(Predictor const& predictor, sgl::vector const& lambda_sequence,
			field<Indices> const& training_samples,
			field<Indices> const& test_samples,
			sgl::natural const number_of_threads) const;
	//Lambda

	/**
	 *
	 * @return
	 */
	sgl::numeric
	lambda_max() const;

	/**
	 *
	 * @param lambda_max
	 * @param lambda_min
	 * @param n
	 * @return
	 */
	sgl::vector
	lambda_sequence(sgl::numeric lambda_max, sgl::numeric lambda_min,
			sgl::natural n) const;
};

//TODO this interface is dangerous as a call changes the state of the objective
template<typename CONFIG, typename ObjectiveFunctionType>
sgl::numeric Interface<CONFIG, ObjectiveFunctionType>::lambda_max() const {

	sgl::vector gradient(sgl.setup.dim);

	typename ObjectiveFunctionType::instance_type objective_function =
			objective_type.create_instance(sgl.setup);

	objective_function.at_zero();
	gradient = objective_function.gradient();

	return sgl.compute_critical_lambda(gradient, alpha);

}

template<typename CONFIG, typename ObjectiveFunctionType>
sgl::vector Interface<CONFIG, ObjectiveFunctionType>::lambda_sequence(
		sgl::numeric lambda_max, sgl::numeric lambda_min,
		sgl::natural n) const {
	sgl::vector lambda_sequence(n);

	lambda_sequence(n - 1) = lambda_min;

	sgl::numeric const a = exp((log(lambda_max) - log(lambda_min)) / (n - 1));

	for (sgl::natural i = 1; i < n; i++) {
		lambda_sequence(n - 1 - i) = a * lambda_sequence(n - i);
	}

	return lambda_sequence;
}

template<typename CONFIG, typename ObjectiveFunctionType>
inline boost::tuple<sgl::block_vector_field, sgl::vector, sgl::vector> Interface<
		CONFIG, ObjectiveFunctionType>::optimize(
		const sgl::vector & lambda_sequence) const {

	//Domain checks
	if (!sgl::is_decreasing(lambda_sequence)
			|| !sgl::is_positive(lambda_sequence)) {
		throw std::domain_error(
				"the lambda sequence must be decreasing and positive");
	}

	sgl::natural_vector needed_solutions(lambda_sequence.n_elem);
	sgl::seq(needed_solutions, 0, 1);

	typename ObjectiveFunctionType::instance_type objective =
			objective_type.create_instance(sgl.setup);

	return optimizer.optimize(objective, lambda_sequence, needed_solutions,
			true);
}

template<typename CONFIG, typename ObjectiveFunctionType>
inline boost::tuple<sgl::block_vector_field, sgl::vector, sgl::vector> Interface<
		CONFIG, ObjectiveFunctionType>::optimize(
		const sgl::vector & lambda_sequence,
		sgl::natural_vector needed_solutions) const {

	//Domain checks
	if (!sgl::is_decreasing(lambda_sequence)
			|| !sgl::is_positive(lambda_sequence)) {
		throw std::domain_error(
				"the lambda sequence must be decreasing and positive");
	}

	//TODO check that all elements of needed_solutions are unique and less than the length of lambda_sequence

	typename ObjectiveFunctionType::instance_type objective =
			objective_type.create_instance();

	return optimizer.optimize(objective, lambda_sequence, needed_solutions,
			true);
}

template<typename CONFIG, typename ObjectiveFunctionType>
inline sgl::natural Interface<CONFIG, ObjectiveFunctionType>::optimize(
		sgl::parameter_field & x_field, sgl::natural_vector needed_solutions,
		sgl::vector & object_value, sgl::vector & function_value,
		const sgl::vector & lambda_sequence) const {

	//Domain checks
	if (!sgl::is_decreasing(lambda_sequence)
			|| !sgl::is_positive(lambda_sequence)) {
		throw std::domain_error(
				"the lambda sequence must be decreasing and positive");
	}

	//TODO check that all elements of needed_solutions are unique and less than the length of lambda_sequence

	typename ObjectiveFunctionType::instance_type objective =
			objective_type.create_instance(sgl.setup);

	return optimizer.optimize(x_field, needed_solutions, object_value, function_value,
			objective, lambda_sequence, true);
}

template<typename CONFIG, typename ObjectiveFunctionType>
template<typename Predictor>
inline boost::tuple<field<typename Predictor::response_type>, sgl::vector,
		sgl::vector> Interface<CONFIG, ObjectiveFunctionType>::regular_cv(
		Predictor const& predictor, sgl::vector const& lambda_sequence,
		field<Indices> const& cv_indices, Indices const& indices_all,
		sgl::natural const number_of_threads) const {

	//Domain checks
	if (!sgl::is_decreasing(lambda_sequence)
			|| !sgl::is_positive(lambda_sequence)) {
		throw std::domain_error(
				"the lambda sequence must be decreasing and positive");
	}

	//TODO indices_all

	//Result matrix
	field<typename Predictor::response_type> response(
			objective_type.data.n_samples, lambda_sequence.n_elem);

	sgl::vector average_number_of_features = zeros<sgl::vector>(lambda_sequence.n_elem);
	sgl::vector average_number_of_parameters = zeros<sgl::vector>(lambda_sequence.n_elem);

	//Training indices
	field<Indices> training_indices(cv_indices.n_elem);
	for (u32 i = 0; i < cv_indices.n_elem; ++i) {
		training_indices(i) = indices_all - cv_indices(i);
	}

	const int n_indices = cv_indices.n_elem;

	bool exception_caught = false;
	string exception_msg;

	// create progress monitor
	Progress p(lambda_sequence.n_elem * cv_indices.n_elem, sgl.config.verbose);

#ifdef SGL_USE_OPENMP
	omp_set_num_threads(number_of_threads);

#pragma omp parallel for schedule(dynamic)
#endif
	for (int i = 0; i < n_indices; i++) {

		if ( ! exception_caught || ! p.is_aborted() ) {

			try {

				ObjectiveFunctionType traning_objective; //Note traning_objective stores the X matrix
				typename ObjectiveFunctionType::data_type test_data;

#ifdef SGL_USE_OPENMP
#pragma omp critical
#endif
				{
					traning_objective.set_data(
							objective_type.data(training_indices(i)));
					test_data = objective_type.data(cv_indices(i));
				}

				typename ObjectiveFunctionType::instance_type objective(
						traning_objective.create_instance(sgl.setup));

				//Fit
				//sgl::block_vector_field x_field(lambda.n_elem);
				sgl::parameter x(sgl.setup);
				sgl::parameter x0(sgl.setup);
				sgl::vector gradient(sgl.setup.dim);

				//Start at zero
				x.zeros();
				objective.at_zero();
				gradient = objective.gradient();

				//Lambda loop
				sgl::natural lambda_index = 0;

				while (true) {

					sgl::numeric const lambda = lambda_sequence(lambda_index);

					optimizer.optimize_single(x, x0, gradient, objective,
							lambda);

					field<typename Predictor::response_type> r =
							predictor.predict(test_data, x);

					//Update average number of features / parameters
#ifdef SGL_USE_OPENMP
#pragma omp critical
#endif
					{
						average_number_of_features(lambda_index) +=
								x.n_nonzero_blocks;
						average_number_of_parameters(lambda_index) +=
								x.n_nonzero;
						//Predict fold

						cv_indices(i).select_rows_view(response).col(
								lambda_index) = r;

					}

					//next lambda
					++lambda_index;
					//Increas progress monitor
					p.increment();

					if (lambda_index >= lambda_sequence.n_elem) {
						//No more lambda values - exit
						break;
					}

					//Go one step back, (avoid computing the gradient) - hence start at x0
					x = x0;
					objective.at(x0);

				}
			} catch (SGL_EXCEPTIONS & ex) {

#ifdef SGL_USE_OPENMP
#pragma omp critical //Needed in the case when tow or more threads throws an exception at the same time
#endif
				{
					if (!exception_caught) {

						//Mark exception caught
						exception_caught = true;

						//Copy msg
						if (ex.what() != NULL) {
							exception_msg = ex.what();
						}

						else {
							exception_msg = "Unknown error";
						}

						//Interrupt all threads
						SGL_INTERRUPT;

					}
				}
			}
		}
	}

	if (exception_caught) {
		SGL_INTERRUPT_RESET;

		//handle exception
		throw std::runtime_error(exception_msg.c_str());
	}

	if(p.is_aborted()) {
		throw std::runtime_error("Aborted by user");
	}

	return boost::make_tuple(response,
			average_number_of_features / static_cast<sgl::numeric>(n_indices),
			average_number_of_parameters / static_cast<sgl::numeric>(n_indices));
}

template<typename CONFIG, typename ObjectiveFunctionType>
template<typename Predictor>
inline boost::tuple<field<field<typename Predictor::response_type> >,
		sgl::natural_matrix, sgl::natural_matrix> Interface<CONFIG,
		ObjectiveFunctionType>::subsampling(Predictor const& predictor,
		sgl::vector const& lambda_sequence,
		field<Indices> const& training_samples,
		field<Indices> const& test_samples,
		sgl::natural const number_of_threads) const {

	//Domain checks
	if (!sgl::is_decreasing(lambda_sequence)
			|| !sgl::is_positive(lambda_sequence)) {
		throw std::domain_error(
				"subsampling : the lambda sequence must be decreasing and positive");
	}

	if (training_samples.n_elem != test_samples.n_elem) {
		throw std::domain_error(
				"subsampling : number of training and test subsamples do not match");
	}

	//TODO domain checks

	sgl::natural n_subsamples = training_samples.n_elem;

	//Result matrix
	field < field<typename Predictor::response_type>
			> response_field_subsamples(n_subsamples);

	sgl::natural_matrix number_of_features(n_subsamples,
			lambda_sequence.n_elem);
	sgl::natural_matrix number_of_parameters(n_subsamples,
			lambda_sequence.n_elem);

	bool exception_caught = false;
	string exception_msg;

	// create progress monitor
	Progress p(lambda_sequence.n_elem * n_subsamples, sgl.config.verbose);


#ifdef SGL_USE_OPENMP
	omp_set_num_threads(number_of_threads);

#pragma omp parallel for schedule(dynamic)
#endif
	for (int i = 0; i < static_cast<int>(n_subsamples); i++) {

		if ( ! exception_caught || ! p.is_aborted() ) {

			try {

				ObjectiveFunctionType traning_objective; //Note traning_objective stores the X matrix
				typename ObjectiveFunctionType::data_type test_data;

#ifdef SGL_USE_OPENMP
#pragma omp critical
#endif
				{
					traning_objective.set_data(
							objective_type.data(training_samples(i)));
					test_data = objective_type.data(test_samples(i));
				}

				typename ObjectiveFunctionType::instance_type objective(
						traning_objective.create_instance(sgl.setup));

				//Response field
				field<typename Predictor::response_type> response_field(
						test_samples(i).size(), lambda_sequence.n_elem);

				//Fit
				sgl::parameter x(sgl.setup);
				sgl::parameter x0(sgl.setup);
				sgl::vector gradient(sgl.setup.dim);

				//Start at zero
				x.zeros();
				x0.zeros();
				objective.at_zero();
				gradient = objective.gradient();

				//Lambda loop
				sgl::natural lambda_index = 0;

				while (true) {

					sgl::numeric const lambda = lambda_sequence(lambda_index);

					optimizer.optimize_single(x, x0, gradient, objective,
							lambda);

					//set number of features / parameters
					number_of_features(i, lambda_index) = x.n_nonzero_blocks;
					number_of_parameters(i, lambda_index) = x.n_nonzero;

					//Predict fold
					response_field.col(lambda_index) = predictor.predict(
							test_data, x);

					//next lambda
					++lambda_index;
					//Increas progress monitor
					p.increment();

					if (lambda_index >= lambda_sequence.n_elem) {
						//No more lambda values - exit
						break;
					}

					//Go one step back, (avoid computing the gradient) - hence start at x0
					x = x0;
					objective.at(x0);

				}

#ifdef SGL_USE_OPENMP
#pragma omp critical
#endif
				{
					response_field_subsamples(i) = response_field;
				}
			} catch (SGL_EXCEPTIONS & ex) {

#ifdef SGL_USE_OPENMP
#pragma omp critical //Needed in the case when tow or more threads throws an exception at the same time
#endif
				{
					if (!exception_caught) {

						//Mark exception caught
						exception_caught = true;

						//Copy msg
						if (ex.what() != NULL) {
							exception_msg = ex.what();
						}

						else {
							exception_msg = "Unknown error";
						}

						//Interrupt all threads
						SGL_INTERRUPT;

					}
				}
			}
		}
	}

	if (exception_caught) {

		SGL_INTERRUPT_RESET;

		//handle exception

		throw std::runtime_error(exception_msg.c_str());
	}

	if(p.is_aborted()) {
		throw std::runtime_error("Aborted by user");
	}

	return boost::make_tuple(response_field_subsamples, number_of_features,
			number_of_parameters);
}

#endif /* INTERFACE_BASIC_H_ */
