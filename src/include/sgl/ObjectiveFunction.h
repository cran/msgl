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

#ifndef OBJECTIVEFUNCTION_H_
#define OBJECTIVEFUNCTION_H_

//TODO no assignment operators

template<typename E>
struct function_traits;

template<typename E>
class ObjectiveFunctionExpression {

public:

	typedef typename function_traits<E>::data_type data_type;

	sgl::natural const n_blocks;
	sgl::natural_vector const block_sizes;

	ObjectiveFunctionExpression(sgl::natural const& n_blocks, sgl::natural_vector const& block_sizes) :
			n_blocks(n_blocks), block_sizes(block_sizes) {

		//Domain check
		if (block_sizes.n_elem != n_blocks) {
			throw std::domain_error("ObjectiveFunctionExpression : Dimension mismatch");
		}

	}

	void at(sgl::parameter const& x) {
		static_cast<E&>(*this).at(x);
	}

	void at_zero() {
		static_cast<E&>(*this).at_zero();
	}

	sgl::numeric evaluate() const {
		return static_cast<E&>(*this).evaluate();
	}

	sgl::vector const gradient() const {
		return static_cast<E&>(*this).gradient();
	}

	sgl::block_vector const gradient(sgl::natural_vector const& indices) const {
		return static_cast<E&>(*this).gradient(indices);
	}

	sgl::matrix const hessian_diag(sgl::natural feature_index) const {
		return static_cast<E&>(*this).hessian_diag(feature_index);
	}

	void hessian_update(sgl::natural feature_index, sgl::parameter_block_vector const& z) {
		static_cast<E&>(*this).hessian_update(feature_index, z);
	}

	sgl::vector const compute_block_gradient(sgl::natural feature_index) const {
		return static_cast<E&>(*this).compute_block_gradient(feature_index);
	}

	sgl::numeric hessian_bound_level0() const {
		return static_cast<E&>(*this).hessian_bound_level0();
	}

	sgl::numeric hessian_bound_level1(sgl::natural block_index) const {
		return static_cast<E&>(*this).hessian_bound_level1(block_index);
	}

	operator E&() {
		return static_cast<E&>(*this);
	}

	operator E const&() const {
		return static_cast<const E&>(*this);
	}
};

template<typename E>
struct function_traits<ObjectiveFunctionExpression<E> > {
	typedef typename E::data_type data_type;
};

template<typename E>
struct ObjectiveFunction_members {

	boost::shared_ptr<E> fun_ptr;

	E & fun;

	ObjectiveFunction_members(boost::shared_ptr<E> fun_ptr) :
			fun_ptr(fun_ptr), fun(*fun_ptr) {
	}

};

template<typename E>
class ObjectiveFunction: private ObjectiveFunction_members<E>, public ObjectiveFunctionExpression<ObjectiveFunction<E> > {

public:

	using ObjectiveFunction_members<E>::fun;

	typedef typename function_traits<ObjectiveFunction>::data_type data_type;

	ObjectiveFunction(data_type const& data, DimConfig const& dim_config) :
			ObjectiveFunction_members<E>(boost::shared_ptr<E>(new E(data, dim_config))), ObjectiveFunctionExpression<ObjectiveFunction<E> >(
					dim_config.n_blocks, dim_config.block_dim) {
	}

	void at(sgl::parameter const& x) {
		fun.at(x);
	}

	void at_zero() {
		fun.at_zero();
	}

	sgl::numeric evaluate() const {
		return fun.evaluate();
	}

	sgl::vector const gradient() const {
		return fun.gradient();
	}

	sgl::block_vector const gradient(sgl::natural_vector const& indices) const {
		return fun.gradient(indices);
	}

	sgl::matrix const hessian_diag(sgl::natural feature_index) const {
		return fun.hessian_diag(feature_index);
	}

	void hessian_update(sgl::natural feature_index, sgl::parameter_block_vector const& z) {
		fun.hessian_update(feature_index, z);
	}

	sgl::vector const compute_block_gradient(sgl::natural feature_index) const {
		return fun.compute_block_gradient(feature_index);
	}

	sgl::numeric hessian_bound_level0() const {
		return fun.hessian_bound_level0();
	}

	sgl::numeric hessian_bound_level1(sgl::natural block_index) const {
		return fun.hessian_bound_level1(block_index);
	}

	static sgl::natural compute_unit_size(data_type const& data) {
		return E::compute_unit_size(data);
	}

	static sgl::natural compute_number_of_units(data_type const& data) {
		return E::compute_number_of_units(data);
	}
};

template<typename E>
struct function_traits<ObjectiveFunction<E> > {
	typedef typename E::data_type data_type;
};

template<typename E, typename Filter>
class DataFilter: public ObjectiveFunctionExpression<DataFilter<E, Filter> > {

public:

	typedef typename function_traits<DataFilter>::data_type data_type;

	Filter const& filter;

	boost::shared_ptr<const data_type> data_ptr;
	data_type const& data;

private:

	E fun;

public:

	template<typename F>
	DataFilter(F const& fun_type, data_type const& data, Filter const& filter, DimConfig const& dim_config) :
			filter(filter), data_ptr(filter.create_filtered_data(data)), data(*data_ptr), fun(
					fun_type.create_instance(this->data, dim_config)) {
	}

	void at(sgl::parameter const& x) {
		fun.at(x);
	}

	void at_zero() {
		fun.at_zero();
	}

	sgl::numeric evaluate() const {
		return fun.evaluate();
	}

	sgl::vector const gradient() const {
		return fun.gradient();
	}

	sgl::block_vector const gradient(sgl::natural_vector const& indices) const {
		return fun.gradient(indices);
	}

	sgl::matrix const hessian_diag(sgl::natural feature_index) const {
		return fun.hessian_diag(feature_index);
	}

	void hessian_update(sgl::natural feature_index, sgl::parameter_block_vector const& z) {
		fun.hessian_update(feature_index, z);
	}

	sgl::vector const compute_block_gradient(sgl::natural feature_index) const {
		return fun.compute_block_gradient(feature_index);
	}

	sgl::numeric hessian_bound_level0() const {
		return fun.hessian_bound_level0();
	}

	sgl::numeric hessian_bound_level1(sgl::natural block_index) const {
		return fun.hessian_bound_level1(block_index);
	}
};

template<typename E, typename Filter>
struct function_traits<DataFilter<E, Filter> > {
	typedef typename E::data_type data_type;
};

template<typename E>
class ObjectiveFunctionScalarMultiplication: public ObjectiveFunctionExpression<ObjectiveFunctionScalarMultiplication<E> > {

private:

	E const fun;
	sgl::numeric const lambda;

public:

	typedef typename function_traits<ObjectiveFunctionScalarMultiplication>::data_type data_type;

	ObjectiveFunctionScalarMultiplication(sgl::numeric lambda, ObjectiveFunctionExpression<E> const& fun) : ObjectiveFunctionExpression<ObjectiveFunctionScalarMultiplication<E> >(fun.n_blocks, fun.block_sizes),
			fun(fun), lambda(lambda) {
	}

	void at(sgl::parameter const& x) {
		const_cast<E&>(fun).at(x);
	}

	void at_zero() {
		const_cast<E&>(fun).at_zero();
	}

	sgl::numeric evaluate() const {
		return lambda * fun.evaluate();
	}

	sgl::vector const gradient() const {
		return lambda * fun.gradient();
	}

	sgl::block_vector const gradient(sgl::natural_vector const& indices) const {
		return lambda * fun.gradient(indices);
	}

	sgl::matrix const hessian_diag(sgl::natural feature_index) const {
		return lambda * fun.hessian_diag(feature_index);
	}

	void hessian_update(sgl::natural feature_index, sgl::parameter_block_vector const& z) {
		const_cast<E&>(fun).hessian_update(feature_index, z);
	}

	sgl::vector const compute_block_gradient(sgl::natural feature_index) const {
		return lambda * fun.compute_block_gradient(feature_index);
	}

	sgl::numeric hessian_bound_level0() const {
		return lambda * fun.hessian_bound_level0();
	}

	sgl::numeric hessian_bound_level1(sgl::natural block_index) const {
		return lambda * fun.hessian_bound_level1(block_index);
	}
};

template<typename E>
struct function_traits<ObjectiveFunctionScalarMultiplication<E> > {
	typedef typename E::data_type data_type;
};

template<typename A, typename B>
class ObjectiveFunctionAddition: public ObjectiveFunctionExpression<ObjectiveFunctionAddition<A, B> > {

private:

	A const fun1;
	B const fun2;

public:

	typedef typename function_traits<ObjectiveFunctionAddition>::data_type data_type;

	ObjectiveFunctionAddition(ObjectiveFunctionExpression<A> const& fun1, ObjectiveFunctionExpression<B> const& fun2) : ObjectiveFunctionExpression<ObjectiveFunctionAddition<A, B> >(fun1.n_blocks, fun1.block_sizes),
			fun1(fun1), fun2(fun2) {

		if(fun1.n_blocks != fun2.n_blocks) {
			throw std::runtime_error("ObjectiveFunctionAddition : dimension mismatch");
		}

		if(accu(fun1.block_sizes != fun2.block_sizes) > 0) {
			throw std::runtime_error("ObjectiveFunctionAddition : dimension mismatch");
		}
	}

	void at(sgl::parameter const& x) {
		const_cast<A&>(fun1).at(x);
		const_cast<B&>(fun2).at(x);
	}

	void at_zero() {
		const_cast<A&>(fun1).at_zero();
		const_cast<B&>(fun2).at_zero();
	}

	sgl::numeric evaluate() const {
		return fun1.evaluate() + fun2.evaluate();
	}

	sgl::vector const gradient() const {
		return fun1.gradient() + fun2.gradient();
	}

	sgl::block_vector const gradient(sgl::natural_vector const& indices) const {
		return fun1.gradient(indices) + fun2.gradient(indices);
	}

	sgl::matrix const hessian_diag(sgl::natural feature_index) const {
		return fun1.hessian_diag(feature_index) + fun2.hessian_diag(feature_index);
	}

	void hessian_update(sgl::natural feature_index, sgl::parameter_block_vector const& z) {
		const_cast<A&>(fun1).hessian_update(feature_index, z);
		const_cast<B&>(fun2).hessian_update(feature_index, z);
	}

	sgl::vector const compute_block_gradient(sgl::natural feature_index) const {
		return fun1.compute_block_gradient(feature_index) + fun2.compute_block_gradient(feature_index);
	}

	sgl::numeric hessian_bound_level0() const {
		return fun1.hessian_bound_level0() + fun2.hessian_bound_level0();
	}

	sgl::numeric hessian_bound_level1(sgl::natural block_index) const {
		return fun1.hessian_bound_level1(block_index) + fun2.hessian_bound_level1(block_index);
	}

};

template<typename A, typename B>
struct function_traits<ObjectiveFunctionAddition<A, B> > {
	typedef typename A::data_type data_type;
};

#endif /* OBJECTIVEFUNCTION_H_ */
