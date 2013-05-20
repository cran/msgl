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
#ifndef ALGORITHMCONFIGURATIONDEFAULT_H_
#define ALGORITHMCONFIGURATIONDEFAULT_H_

class AlgorithmConfigurationDefault {

public:

	sgl::numeric tolerance_penalized_main_equation_loop;

	sgl::numeric tolerance_penalized_inner_loop_alpha;
	sgl::numeric tolerance_penalized_inner_loop_beta;

	sgl::numeric tolerance_penalized_middel_loop_alpha;

	sgl::numeric tolerance_penalized_outer_loop_alpha;
	sgl::numeric tolerance_penalized_outer_loop_beta;
	sgl::numeric tolerance_penalized_outer_loop_gamma;

	bool use_bound_optimization;

	bool use_stepsize_optimization_in_penalizeed_loop;
	sgl::numeric stepsize_opt_penalized_initial_t;
	sgl::numeric stepsize_opt_penalized_a;
	sgl::numeric stepsize_opt_penalized_b;

	bool verbose;

	AlgorithmConfigurationDefault() :
		tolerance_penalized_main_equation_loop(1e-10),

		tolerance_penalized_inner_loop_alpha(1e-4),
		tolerance_penalized_inner_loop_beta(0),

		tolerance_penalized_middel_loop_alpha(0.01),

		tolerance_penalized_outer_loop_alpha(0.01),
		tolerance_penalized_outer_loop_beta(0),
		tolerance_penalized_outer_loop_gamma(5e-4),

		use_bound_optimization(true),

		use_stepsize_optimization_in_penalizeed_loop(true),
		stepsize_opt_penalized_initial_t(1),
		stepsize_opt_penalized_a(0.1),
		stepsize_opt_penalized_b(0.5),

		verbose(true) {
	}

};

#endif /* ALGORITHMCONFIGURATIONDEFAULT_H_ */
