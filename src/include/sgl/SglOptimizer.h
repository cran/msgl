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

#ifndef SGLOPTIMIZER_H_
#define SGLOPTIMIZER_H_

template<typename SGL>
class SglOptimizer {

private:

	SGL const& sgl;

public:

	const sgl::numeric alpha;

	SglOptimizer(SGL const& sgl, const sgl::numeric alpha);

	//Will only return the solutions (i.e x's and objective values ect) for the lambda indices values given by solution_index
	template<typename T>
	boost::tuple<sgl::block_vector_field, sgl::vector, sgl::vector> optimize(T & objective, sgl::vector const& lambda_sequence,
			sgl::natural_vector const& solution_index, bool handle_exceptions = false) const;

	template<typename T>
	void optimize(sgl::parameter_field & x_field, sgl::natural_vector const& x_field_index, sgl::vector & object_value,
			sgl::vector & function_value, T & objective, sgl::vector const& lambda_sequence, bool handle_exceptions = false) const;

	template<typename T>
	sgl::parameter optimize_single(sgl::parameter & x, sgl::parameter & x_old, sgl::vector & gradient, T & objective, sgl::numeric const& lambda) const;

private:

	template<typename T>
	void optimize_quadratic(T & objective, sgl::parameter & x, sgl::vector const& gradient, sgl::vector const& critical_bounds,
			sgl::numeric const alpha, sgl::numeric const lambda) const;

	void optimize_inner(sgl::vector const& gradient_at_x0, sgl::matrix const& second_order_term, sgl::numeric penalty_constant_L2,
			sgl::vector const& penalty_constant_L1, sgl::parameter_block_vector & x, sgl::parameter_block_vector const& x0) const;

	//Approximately Solve c + h*x + penalty*x/sqrt(x*x+r) = 0, returns the absolute value of the result
	sgl::numeric solve_main_equation(sgl::numeric const c, sgl::numeric const h, sgl::numeric const p, sgl::numeric const r,
			sgl::numeric const x) const;

	sgl::numeric update_x(sgl::numeric const g, sgl::numeric const h, sgl::numeric const penalty_constant_L2,
			sgl::vector const& penalty_constant_L1, sgl::numeric const& x, sgl::numeric const r, sgl::natural const i) const;

	//TODO move to SglProblem
	void argmin_subgradient(sgl::parameter_block_vector & x, sgl::vector const& v, sgl::vector const& p) const;

	//TODO move
	sgl::numeric function_value(sgl::parameter_block_vector const& x, sgl::vector const& gradient_at_x, sgl::matrix const& second_order_term,
			sgl::numeric penalty_constant_L2, sgl::vector const& penalty_constant_L1) const;

	void locate_safe_point(sgl::parameter_block_vector & safe_point, sgl::parameter_block_vector const& x, sgl::vector const& gradient_at_x,
			sgl::matrix const& second_order_term, sgl::numeric penalty_constant_L2, sgl::vector const& penalty_constant_L1) const;

	template<typename T>
	sgl::numeric
	stepsize_optimize_penalized(T & objective, sgl::parameter const& x1, sgl::parameter const& x0, sgl::vector const& gradient_at_x0,
			sgl::numeric const likelihood_at_x0, sgl::numeric const lambda) const;

	//TODO move
	bool is_stopping_criteria_fulfilled(sgl::parameter const& x, sgl::parameter const& x_old, sgl::numeric const f,
			sgl::numeric const f_old) const;

	void interrupt_check() const;

	void init_interrupt_handler();
};

template<typename SGL>
inline SglOptimizer<SGL>::SglOptimizer(const SGL & sgl, const sgl::numeric alpha) : sgl(sgl), alpha(alpha) {
  if (0 > alpha || alpha > 1) {
      throw std::domain_error("alpha is not in the range 0 to 1");
  }
}

template<typename SGL>
template<typename T>
inline boost::tuple<sgl::block_vector_field,sgl::vector,sgl::vector>
SglOptimizer<SGL>::optimize(T & objective, const sgl::vector & lambda_sequence, const sgl::natural_vector & solution_index, bool handle_exceptions) const
{
    sgl::block_vector_field x_field(solution_index.n_elem);
    sgl::vector object_value(solution_index.n_elem);
    sgl::vector function_value(solution_index.n_elem);

    optimize(x_field, solution_index, object_value, function_value, objective, lambda_sequence, handle_exceptions);

    return boost::make_tuple(x_field, object_value, function_value);
}

//TODO remove gradinet_old
//TODO stopping criteria an object of it own
template<typename SGL>
bool SglOptimizer<SGL>::is_stopping_criteria_fulfilled(sgl::parameter const& x, sgl::parameter const& x_old, sgl::numeric const f,
		sgl::numeric const f_old) const {

	TIMER_START;

	//Check stopping criteria

	//TODO Switch on/off layers in config

	//Layer : Objective function value stopping rule

	if (abs(f_old - f) / abs(f_old) > sgl.config.tolerance_penalized_outer_loop_gamma) {

		return false;
	}

	//Layer : Distance stopping rule
	sgl::numeric d = sgl.dist(x_old, x);

	if (d > sgl.config.tolerance_penalized_outer_loop_alpha) {
		return false;
	}

	//Layer : Model stopping rule

	if (sgl.discrete_dist(x_old, x) > sgl.config.tolerance_penalized_outer_loop_beta) {
		return false;
	}

	return true;

}

//Given x and gradient should correspond when called
// x and gradient will not correspond after the call (for the k'th x, gradient will be the k-1'th)
//Returns the k-1'th parameter - corresponding to last updated gradient
//x and x_old must have same size
template<typename SGL>
template<typename T>
sgl::parameter SglOptimizer<SGL>::optimize_single(sgl::parameter & x, sgl::parameter & x_old, sgl::vector & gradient, T & objective,
		sgl::numeric const& lambda) const {

	TIMER_START;

	//Quadratic approximation loop

	sgl::numeric f;
	sgl::numeric f_old;

	//Initialise converges checker
	CONVERGENCE_CHECK(1e3); //TODO configable convergens checker

	sgl::vector critical_bounds(sgl.setup.n_blocks); //TODO this will only be used when use_bound_optimization = true

	do {

		CONVERGENCE_CHECK_INCREASE;

		x_old = x;
		f_old = objective.evaluate() + sgl.penalty(x, alpha, lambda);

		//Compute critical bounds
		if (sgl.config.use_bound_optimization) {
			critical_bounds = sgl.compute_bounds(gradient, x, alpha, lambda);
		}

		optimize_quadratic(objective, x, gradient, critical_bounds, alpha, lambda);

		objective.at(x);
		f = objective.evaluate() + sgl.penalty(x, alpha, lambda);

		// Stepsize optimization
		//TODO use_stepsize_optimization_in_penalizeed_loop levels - level 0 - never use step size opt - level 1 only when f > f_old - level 2 always
		if (f > f_old && sgl.config.use_stepsize_optimization_in_penalizeed_loop && !is_stopping_criteria_fulfilled(x, x_old, f, f_old)) {
			//Start scope timer, note will only be activated if SGL_TIMING is defined
			TIMER_START;

			objective.at(x_old);

			sgl::numeric t = stepsize_optimize_penalized(objective, x, x_old, gradient, objective.evaluate(), lambda);

			if (t != 1) {
				x = (1 - t) * x_old + t * x;
			}

#ifdef SGL_DEBUG_INFO_STEPSIZE
			rout << "stepsize = " << t << std::endl;
#endif

			objective.at(x);
			f = objective.evaluate() + sgl.penalty(x, alpha, lambda);

		}

#ifdef SGL_DEBUG_INFO_QUADRATIC
		rout << " parameter distance = " << sgl.dist(x_old, x) << std::endl;
		rout << " parameter discrete distance = " << sgl.discrete_dist(x_old, x) << std::endl;
		rout << " function value = " << f << std::endl;
#endif

		if (!is_stopping_criteria_fulfilled(x, x_old, f, f_old)) {

			//Continue quadratic loop

			ASSERT_IS_NON_NEGATIVE(f_old - f + std::numeric_limits<sgl::numeric>::epsilon());

			gradient = objective.gradient();

		} else {

			//Exit quadratic loop
			return x_old;
		}

	} while (true);

}

template<typename SGL>
template<typename T>
void SglOptimizer<SGL>::optimize(sgl::parameter_field & x_field, sgl::natural_vector const& needed_solutions, sgl::vector & object_value,
		sgl::vector & function_value, T & objective, sgl::vector const& lambda_sequence, bool handle_exceptions) const {

	//Start scope timer, note will only be activated if SGL_TIMING is defined
	TIMER_START;

	//Ensure x_field_index is ordered
	sgl::natural_vector x_field_order = sort(needed_solutions, 0);

	sgl::vector gradient(sgl.setup.dim);
	sgl::parameter x(sgl.setup);
	sgl::parameter x0(sgl.setup);

	//Start at zero
	x.zeros();
	x0.zeros();
	objective.at_zero();
	gradient = objective.gradient();

	//Lambda loop
	sgl::natural lambda_index = 0;

	//Needed solutions index
	sgl::natural x_field_index = 0;

	try {
		while (true) {

			sgl::numeric const lambda = lambda_sequence(lambda_index);

			if (sgl.config.verbose) {
				SGL_STD_OUT << "At index " << lambda_index << " - lambda = " << lambda << " - obj. fun. value = " << objective.evaluate()
								<< " - non zero blocks = " << x.n_nonzero_blocks << " - non zero parameters "
							<< x.n_nonzero << endl;
			}

			 optimize_single(x, x0, gradient, objective, lambda);

			//Check if we must save the solution
			if (lambda_index == x_field_order(x_field_index)) {

				//Save solution
				//TODO remove x_field(x_field_index).copy_size(x);
				x_field(x_field_index) = x;

				//Save objective function values
				object_value(x_field_index) = objective.evaluate();
				function_value(x_field_index) = object_value(x_field_index) + sgl.penalty(x, alpha, lambda);

				//Next x field
				++x_field_index;
			}

			//next lambda
			++lambda_index;

			if (lambda_index >= lambda_sequence.n_elem || x_field_index >= x_field_order.n_elem) {
				//No more lambda values or no more solutions needed - exit
				break;
			}

			//Go one step back, (avoid computing the gradient) - hence start at x0
			x = x0;
			objective.at(x0);

		}

	} catch (SGL_EXCEPTIONS & e) {

		//Reset intrrupt flag
		SGL_INTERRUPT_RESET;

		throw;

	} catch (SGL_EXCEPTIONS_GENRAL) {

		//Reset intrrupt flag
		SGL_INTERRUPT_RESET;

		throw;
	}
}

template<typename SGL>
template<typename T>
sgl::numeric SglOptimizer<SGL>::stepsize_optimize_penalized(T & objective, sgl::parameter const& x1, sgl::parameter const& x0,
		sgl::vector const& gradient_at_x0, sgl::numeric const likelihood_at_x0, sgl::numeric const lambda) const {

	TIMER_START;

	sgl::numeric t = sgl.config.stepsize_opt_penalized_initial_t;

	sgl::numeric penalty_at_x0 = sgl.penalty(x0, alpha, lambda);

	sgl::numeric delta = sgl.config.stepsize_opt_penalized_a
			* (dot(gradient_at_x0, (x1 - x0)) + sgl.penalty(x1, alpha, lambda) - penalty_at_x0);

	sgl::numeric value_at_x0 = likelihood_at_x0 + penalty_at_x0;

	sgl::numeric value;
	while (objective.at((1 - t) * x0 + t * x1), value = objective.evaluate() + sgl.penalty((1 - t) * x0 + t * x1, alpha, lambda), value
			> value_at_x0 + t * delta && t > 0) {

		t *= sgl.config.stepsize_opt_penalized_b;

		if(t - std::numeric_limits<sgl::numeric>::epsilon() < 0) {
		    //TODO what should we do in this case
		    return(0);
		}
	}

	ASSERT_IS_POSITIVE(t);

	return t;
}

//gradient should be the gradient at x
template<typename SGL>
template<typename T>
inline void SglOptimizer<SGL>::optimize_quadratic(T & objective, sgl::parameter & x, sgl::vector const& gradient,
		sgl::vector const& critical_bounds, sgl::numeric const alpha, sgl::numeric const lambda) const {

	//Initialize converges checker
	CONVERGENCE_CHECK(1e4); //TODO configable converges checker

	sgl::vector block_gradient;
	sgl::parameter_block_vector x_new;

        sgl::numeric dist;

        do {

		TIMER_START;

		CONVERGENCE_CHECK_INCREASE;

#ifdef SGL_DEBUG_INFO_GB_OPT
		sgl::natural computed_gbs = 0;
#endif

		// **** Regular coordinate decent loop

		// Reset distance
		dist = 0;

		for (sgl::natural block_index = 0; block_index < sgl.setup.n_blocks; ++block_index) {

			block_gradient.set_size(sgl.setup.block_dim(block_index));
			x_new.set_size(sgl.setup.block_dim(block_index));

			// **** Block optimization

			sgl::natural block_start = sgl.setup.block_start_index(block_index);
			sgl::natural block_end = sgl.setup.block_end_index(block_index);

			if (!sgl.config.use_bound_optimization
					|| (!x.is_block_zero(block_index)
					    || (critical_bounds(block_index) <= objective.hessian_bound_level0()
									&& critical_bounds(block_index) <= objective.hessian_bound_level1(block_index)))) {

				//Block could be active, check needed

#ifdef SGL_DEBUG_INFO_GB_OPT
				++computed_gbs;
#endif

				block_gradient = gradient.subvec(block_start, block_end) + objective.compute_block_gradient(block_index);

				//Check if block is active

				bool block_is_active;
				if (x.is_block_zero(block_index)) {

					block_is_active = sgl.is_block_active(block_gradient, block_index, alpha, lambda);

				} else {

			            block_is_active = sgl.is_block_active(block_gradient - objective.hessian_diag(block_index)*x.block(block_index), block_index, alpha, lambda);

				}

				if (block_is_active) {

				        sgl::parameter_block_vector x_block(x.block(block_index));

				        //Block active (non zero), Optimise
					optimize_inner(block_gradient, objective.hessian_diag(block_index),
							lambda * (1 - alpha) * sgl.setup.L2_penalty_weight(block_index),
							lambda * alpha * sgl.setup.L1_penalty_weight(block_index), x_new, x_block);

					// Update

					double const dist_block = sgl.max_dist(x_block, x_new);

					// Max dist
					if (dist < dist_block) {
						dist = dist_block;
					}

					//TODO remove
                                        //double const discrete_dist_block = sgl::discrete_dist(x_block, x_new);
//					if (discrete_dist < discrete_dist_block) {
//						discrete_dist = discrete_dist_block;
//					}

					objective.hessian_update(block_index, x_new);

					x.set_block(block_index, x_new);

					continue;

				}

			}

			//Block is inactive

#ifdef SGL_DEBUG_BLOCK_ACTIVE

			// This is a complex debug section,
			// if compiled with this section active performance will be considerably lower than optimal performance

			block_gradient = gradient.subvec(block_start, block_end) + objective.compute_block_gradient(block_index);

			if (!x.is_block_zero(block_index)) {
				block_gradient -= objective.hessian_diag(block_index) * x.block(block_index);
			}

			//Check if block is active
			if (sgl.is_block_active(block_gradient, block_index, alpha, lambda)) {

				rout << "block = " << block_index << " gab = " << sgl.compute_K(abs(block_gradient) - lambda * alpha * sgl.setup.L1_penalty_weight(block_index), 0) - sgl::square(lambda * (1 - alpha) * sgl.setup.L2_penalty_weight(block_index)) << endl;

				rout << "critical bound = " << critical_bounds(block_index) << " hessian level 0 bound = " << objective.hessian_bound_level0()
								<< " hessian level 1 bound = " << objective.hessian_bound_level1(block_index) << endl;

				//throw std::runtime_error("error - hessian bound");
			}
#endif

			if (x.is_block_zero(block_index)) {
				continue;
			}

			// Update
			x_new.zeros();

                        sgl::parameter_block_vector x_block(x.block(block_index));
			double const dist_block = sgl.max_dist(x_block, x_new);

			// Max dist
			if (dist < dist_block) {
				dist = dist_block;
			}

			objective.hessian_update(block_index, x_new);
			x.set_block_zero(block_index);
		}

#ifdef SGL_DEBUG_INFO_GB_OPT
		rout << "Computed block gradients " << computed_gbs << std::endl;
#endif

#ifdef SGL_DEBUG_QUADRATIC_STOPPING
		rout << "Quadratic loop - dist = " << dist << std::endl;
#endif

	} while (dist > sgl.config.tolerance_penalized_middel_loop_alpha);

}

template<typename SGL>
void SglOptimizer<SGL>::optimize_inner(sgl::vector const& gradient_at_x0, sgl::matrix const& second_order_term,
		sgl::numeric penalty_constant_L2, sgl::vector const& penalty_constant_L1, sgl::parameter_block_vector & x,
		sgl::parameter_block_vector const& x0) const {

	//Start timer, only active if SGL_TIMING is defined
	TIMER_START;

	//Initialise converges checker
	CONVERGENCE_CHECK(1e5); //TODO configable

	const sgl::natural dim = x0.n_elem;

	sgl::vector gradient = gradient_at_x0;
	x = x0; //start at x0

	sgl::numeric sumsq = as_scalar(sum(square(x)));
	sgl::parameter_block_vector x_old(dim);

	do {
		CONVERGENCE_CHECK_INCREASE;
		SGL_INTERRUPT_CHECK;

		x_old = x;

		for (sgl::natural i = 0; i < dim; i++) {

			sgl::numeric const xi = x(i);

			//Compute new x
			sgl::numeric x_new = update_x(gradient(i), second_order_term(i, i), penalty_constant_L2, penalty_constant_L1, xi,
					sgl::pos(sumsq - sgl::square(xi)), i);

			ASSERT_IS_NUMBER(x_new);
			ASSERT_IS_FINITE(x_new);

			//Update gradient and x
			if (x_new != xi) {
				gradient += (x_new - xi) * second_order_term.col(i);
				sumsq += sgl::square(x_new) - sgl::square(xi);

				x(i) = x_new;

			}

		}

		//Check if we ended up near zero
		if (sumsq < 1e-20 && function_value(x, gradient, second_order_term, penalty_constant_L2, penalty_constant_L1) >= 0) {
			//Find new x
			locate_safe_point(x, x0, gradient_at_x0, second_order_term, penalty_constant_L2, penalty_constant_L1);
			gradient = gradient_at_x0 + second_order_term * (x - x0);
		}

	} while (sgl.max_dist(x_old, x) > sgl.config.tolerance_penalized_inner_loop_alpha
			|| sgl::discrete_dist(x_old, x) > sgl.config.tolerance_penalized_inner_loop_beta);

	ASSERT_IS_NON_ZERO(x);
}

//TODO better name for function
template<typename SGL>
void SglOptimizer<SGL>::locate_safe_point(sgl::parameter_block_vector & safe_point, sgl::parameter_block_vector const& x,
		sgl::vector const& gradient_at_x, sgl::matrix const& second_order_term, sgl::numeric penalty_constant_L2,
		sgl::vector const& penalty_constant_L1) const {

	sgl::parameter_block_vector x_decent(safe_point.n_elem);
	argmin_subgradient(x_decent, gradient_at_x - second_order_term * x, penalty_constant_L1);

	safe_point = x_decent;

	sgl::numeric t = 1; //TODO configable

	sgl::numeric value;
	while (value = function_value(safe_point, gradient_at_x + second_order_term * (safe_point - x), second_order_term, penalty_constant_L2,
			penalty_constant_L1), value >= 0) {

		t = 0.9 * t;

		safe_point = t * x_decent;
	}
}

template<typename SGL>
sgl::numeric SglOptimizer<SGL>::function_value(sgl::parameter_block_vector const& x, sgl::vector const& gradient_at_x,
		sgl::matrix const& second_order_term, sgl::numeric penalty_constant_L2, sgl::vector const& penalty_constant_L1) const {

	TIMER_START;

	return arma::as_scalar(
			trans(gradient_at_x - 1 / 2 * second_order_term * x) * x + penalty_constant_L2 * sqrt(sum(square(x)))
					+ sum(penalty_constant_L1 % abs(x)));
}

template<typename SGL>
void SglOptimizer<SGL>::argmin_subgradient(sgl::parameter_block_vector & x, sgl::vector const& v, sgl::vector const& p) const {

	for (sgl::natural i = 0; i < v.n_elem; ++i) {
		if (abs(v(i)) > p(i)) {
			x(i) = p(i) * sgl::sign(v(i))-v(i);
		} else {
			x(i) = 0;
		}
	}
}

//r = sum_{j != i} x_j^2
template<typename SGL>
sgl::numeric SglOptimizer<SGL>::update_x(sgl::numeric g, sgl::numeric const h, sgl::numeric const penalty_constant_L2,
		sgl::vector const& penalty_constant_L1, sgl::numeric const& x, sgl::numeric const r, sgl::natural const i) const {

	if (h == 0) {
		return 0;
	}

	//Special case no penalty
	if (penalty_constant_L1(i) == 0 && penalty_constant_L2 == 0) {
		return x - g / h;
	}

	//Special case L2 penalty = 0
	if (penalty_constant_L2 == 0) {

		sgl::numeric const penalty = penalty_constant_L1(i);

		if (abs(g - h * x) <= penalty) {
			return 0;
		}

		if (g - h * x < -penalty) {
			return x - (g + penalty) / h;
		}

		return x + (penalty - g) / h;
	}

	//Compute r
//	sgl::vector x_temp = x_vector;
//	x_temp.shed_row(i);
//	sgl::numeric const r = arma::as_scalar(sum(square(x_temp)));
//
//	ASSERT_IS_FINITE(r);

	//Case r = 0
	if (r == 0) {

		sgl::numeric const penalty = penalty_constant_L1(i) + penalty_constant_L2;

		if (abs(g - h * x) <= penalty) {
			return 0;
		}

		if (g - h * x < -penalty) {
			return x - (penalty + g) / h;
		}

		return x + (penalty - g) / h;
	}

	//Case r non zero

	//Special case L1 penalty = 0
	if (penalty_constant_L1(i) == 0) {

		sgl::numeric const penalty = penalty_constant_L2;

		return -sgl::sign(g - h * x) * solve_main_equation(abs(g - h * x), h, penalty, r, x);
	}

	//L1 penalty and L2 penalty non zero
	sgl::numeric const penalty = penalty_constant_L1(i);

	if (abs(g - h * x) <= penalty) {
		return 0;
	}

	if (g - h * x < -penalty) {

		//x pos
		return solve_main_equation(abs(g - h * x + penalty_constant_L1(i)), h, penalty_constant_L2, r, x);
	}

	//x negative
	return -solve_main_equation(abs(g - h * x - penalty_constant_L1(i)), h, penalty_constant_L2, r, x);

}

// solver for c + h*x + p*x/sqrt(x^2+r) = 0
template<typename SGL>
sgl::numeric SglOptimizer<SGL>::solve_main_equation(sgl::numeric const c, sgl::numeric const h, sgl::numeric const p, sgl::numeric const r,
		sgl::numeric const x_initial) const {

	//Start scope timer, note will only be activated if SGL_TIMING is defined
	TIMER_START;

	//Initialise converges checker
	CONVERGENCE_CHECK(1e8); //TODO configable

	//Domain checks
	ASSERT_IS_POSITIVE(c);
	ASSERT_IS_POSITIVE(h);
	ASSERT_IS_POSITIVE(p);
	ASSERT_IS_POSITIVE(r);

	sgl::numeric x0 = 0;
	sgl::numeric x1 = -c / h;

	//use initial point

	sgl::numeric const x_init = -abs(x_initial);

	if (x_init > x1) {
		if (c + h * x_init + p * x_init / sqrt(x_init * x_init + r) > 0) {
			x0 = x_init;
		} else {
			x1 = x_init;
		}
	}

	do {
		CONVERGENCE_CHECK_INCREASE;

		sgl::numeric new_x = x1 + (x0 - x1) / 2;

		sgl::numeric value = c + h * new_x + p * new_x / sqrt(new_x * new_x + r);

		if (abs(value) < 1e-10) { //TODO configurable
			x0 = new_x;
			break;
		} else if (value > 0) {
			x0 = new_x;
		} else {
			x1 = new_x;
		}

	} while (abs(x0 - x1) > sgl.config.tolerance_penalized_main_equation_loop);

	ASSERT_IS_FINITE(x0);
	return abs(x0);
}



#endif /* SGLOPTIMIZER_H_ */
