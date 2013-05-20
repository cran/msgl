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

#ifndef MSGL_MULTINOMIAL_PREDICTOR_H_
#define MSGL_MULTINOMIAL_PREDICTOR_H_


template<typename T, typename R, typename S>
class Predictor {

public:

	typedef MatrixData<T> data_type;
	typedef S response_type;

	//TODO remove
//	inline const field<typename Predictor::response_type> predict(const data_type & sample_data, Indices const& samples,
//			const sgl::parameter & parameter) const {
//
//		TIMER_START
//
//		return do_predict(rows(sample_data.data_matrix, samples.getElements()), parameter);
//	}
//
//	inline const field<typename Predictor::response_type> predict(const data_type & sample_data, Indices const& samples,
//			const sgl::sparse_matrix_field & parameters) const {
//
//		TIMER_START
//
//		return do_predict(row_subview(sample_data.data_matrix, samples.getElements()), parameters);
//	}

	inline const field<typename Predictor::response_type> predict(const data_type & sample_data,
			const sgl::sparse_matrix_field & parameters) const {

		TIMER_START

		return do_predict(sample_data.data_matrix, parameters);
	}

	inline const field<typename Predictor::response_type> predict(const data_type & sample_data,
			const sgl::parameter & parameters) const {

		TIMER_START

		return do_predict(sample_data.data_matrix, parameters);
	}

private:

	template<typename E>
	field<response_type> const do_predict(E const& sample_data_matrix, const sgl::sparse_matrix_field & parameters) const {

		field<response_type> response(sample_data_matrix.n_rows, parameters.n_elem);

		for (sgl::natural j = 0; j < parameters.n_elem; ++j) {

			response.col(j) = do_predict(sample_data_matrix, parameters(j));
		}

		return response;

	}

	template<typename E>
	field<response_type> const do_predict(E const& sample_data_matrix, const sgl::sparse_matrix & parameter) const {

		sgl::natural n_samples = sample_data_matrix.n_rows;

		field<response_type> response(n_samples);

		sgl::matrix lp(sample_data_matrix);
		lp = parameter * trans(lp);

		for (sgl::natural i = 0; i < n_samples; ++i) {
			response(i) = static_cast<R const&>(*this).createResponse(lp.col(i));
		}

		return response;

	}

};

template<typename T>
class MultinomialPredictor : public Predictor<T, MultinomialPredictor<T>, MultinomialResponse> {

public:

	typedef typename Predictor<T, MultinomialPredictor<T>, MultinomialResponse>::response_type response_type;

	template<typename E>
	response_type createResponse(E const& linear_predictors) const {
		return MultinomialResponse(static_cast<sgl::vector>(linear_predictors));
	}

};

template<typename R>
static boost::tuple<sgl::matrix_field, sgl::matrix_field, sgl::natural_matrix> convert(field<R> const& field_obj) {

	sgl::natural number_of_samples = field_obj.n_rows;
	sgl::natural length_of_lambda = field_obj.n_cols;
	sgl::natural number_of_groups = field_obj(0, 0).number_of_classes();

	sgl::matrix_field link(length_of_lambda);
	sgl::matrix_field response(length_of_lambda);
	sgl::natural_matrix classes(number_of_samples, length_of_lambda);

	for (sgl::natural i = 0; i < length_of_lambda; ++i) {

		link(i).set_size(field_obj(0, i).linear_predictor().n_elem, number_of_samples);
		response(i).set_size(number_of_groups, number_of_samples);

		for (sgl::natural j = 0; j < number_of_samples; ++j) {

			link(i).col(j) = field_obj(j, i).linear_predictor();
			response(i).col(j) = field_obj(j, i).response();
			classes(j, i) = field_obj(j, i).predicted_class();

		}

	}

	return boost::make_tuple(link, response, classes);
}

#endif /* MSGL_MULTINOMIAL_PREDICTOR_H_ */
