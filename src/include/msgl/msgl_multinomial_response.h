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

#ifndef MSGL_MULTINOMIAL_RESPONSE_H_
#define MSGL_MULTINOMIAL_RESPONSE_H_


class MultinomialResponse {

private:

	//TODO make public const members
	sgl::natural n_classes;

	sgl::natural p_predicted_class;
	sgl::vector p_linear_predictors;
	sgl::vector p_probabilities;

public:

	MultinomialResponse(sgl::vector const& linear_predictors) :
			n_classes(linear_predictors.n_elem), p_predicted_class(argmax(linear_predictors)), p_linear_predictors(
					linear_predictors), p_probabilities(exp(linear_predictors) * (1/ sum(exp(linear_predictors)))) {
	}

	//Needed so that we can use fields
	MultinomialResponse() :
			n_classes(0), p_predicted_class(0), p_linear_predictors(), p_probabilities() {
	}

    sgl::natural number_of_classes() const;

    sgl::vector linear_predictor() const;

    sgl::natural predicted_class() const;

    sgl::vector response() const;

};

sgl::natural MultinomialResponse::number_of_classes() const
{
    return n_classes;
}

sgl::vector MultinomialResponse::linear_predictor() const
{
    return p_linear_predictors;
}

sgl::natural MultinomialResponse::predicted_class() const
{
    return p_predicted_class;
}

sgl::vector MultinomialResponse::response() const
{
    return p_probabilities;
}

#endif /* MSGL_MULTINOMIAL_RESPONSE_H_ */
