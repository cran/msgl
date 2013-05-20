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

#ifndef MSGL_ALGORITHM_CONFIG_H_
#define MSGL_ALGORITHM_CONFIG_H_

template<typename type>
static type getConfigAttribute(rList const& config, std::string const& name) {

	int index;
	if (index = config.getIndex(name), index >= 0) {

		return get_value < type > (config.get(index));

	} else {

		string msg = "Missing configuration parameter : ";
		throw std::domain_error(msg.append(name).c_str());
		return type(); //avoid compiler warnings
	}
}

class AlgorithmConfiguration {

public:

	sgl::numeric const tolerance_penalized_main_equation_loop;

	sgl::numeric const tolerance_penalized_inner_loop_alpha;
	sgl::numeric const tolerance_penalized_inner_loop_beta;

	sgl::numeric const tolerance_penalized_middel_loop_alpha;

	sgl::numeric const tolerance_penalized_outer_loop_alpha;
	sgl::numeric const tolerance_penalized_outer_loop_beta;
	sgl::numeric const tolerance_penalized_outer_loop_gamma;

	bool const use_bound_optimization;

	bool const use_stepsize_optimization_in_penalizeed_loop;
	sgl::numeric const stepsize_opt_penalized_initial_t;
	sgl::numeric const stepsize_opt_penalized_a;
	sgl::numeric const stepsize_opt_penalized_b;

	bool const verbose;

	AlgorithmConfiguration(rList const& config) :

			tolerance_penalized_main_equation_loop(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_main_equation_loop")),

			tolerance_penalized_inner_loop_alpha(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_inner_loop_alpha")),

			tolerance_penalized_inner_loop_beta(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_inner_loop_beta")),

			tolerance_penalized_middel_loop_alpha(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_middel_loop_alpha")),

			tolerance_penalized_outer_loop_alpha(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_outer_loop_alpha")),

			tolerance_penalized_outer_loop_beta(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_outer_loop_beta")),

			tolerance_penalized_outer_loop_gamma(
					getConfigAttribute<sgl::numeric>(config,
							"tolerance_penalized_outer_loop_gamma")),

			use_bound_optimization(
					getConfigAttribute<bool>(config,
							"use_bound_optimization")),

			use_stepsize_optimization_in_penalizeed_loop(
					getConfigAttribute<bool>(config,
							"use_stepsize_optimization_in_penalizeed_loop")),

			stepsize_opt_penalized_initial_t(
					getConfigAttribute<sgl::numeric>(config,
							"stepsize_opt_penalized_initial_t")),

			stepsize_opt_penalized_a(
					getConfigAttribute<sgl::numeric>(config,
							"stepsize_opt_penalized_a")),

			stepsize_opt_penalized_b(
					getConfigAttribute<sgl::numeric>(config,
							"stepsize_opt_penalized_b")),

			verbose(getConfigAttribute<bool>(config, "verbose")) {
	}
};

#endif /* MSGL_ALGORITHM_CONFIG_H_ */
