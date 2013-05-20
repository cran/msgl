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

#ifndef OBJECTIVEFUNCTIONEXPRESSIONTYPE_H_
#define OBJECTIVEFUNCTIONEXPRESSIONTYPE_H_

//TODO REMOVE
//class NoData {};
//
//static NoData noData;

template <typename E>
struct type_traits;

template< typename E>
class ObjectiveFunctionExpressionType {

public:

	typedef typename type_traits<E>::data_type data_type;
	typedef typename type_traits<E>::instance_type instance_type;

	data_type const& data;

	ObjectiveFunctionExpressionType(data_type const& data) : data(data) {}

	instance_type create_instance(sgl::DimConfig const& dim_config) const {
		return static_cast<E&> (*this).create_instance(dim_config);
	}

	ObjectiveFunctionExpressionType<E> operator()(Indices const& indices) const  {
		return static_cast<E&> (*this)(indices);
	}

	operator E&() {
		return static_cast<E&> (*this);
	}

	operator E const&() const {
		return static_cast<const E&> (*this);
	}
};

template <typename E>
struct type_traits<ObjectiveFunctionExpressionType<E> > {
	typedef typename E::data_type data_type;
	typedef typename E::instance_type instance_type;
};

//TODO do we need filters
template<typename E, typename Filter>
class DataFilterType : public ObjectiveFunctionExpressionType<DataFilterType<E, Filter> > {

public:

	typedef Filter filter_type;

	typedef typename type_traits<DataFilterType>::data_type data_type;
	typedef typename type_traits<DataFilterType>::instance_type instance_type;

private:

	Filter const& filter;

	E const fun_type;

public:

	using ObjectiveFunctionExpressionType<DataFilterType<E, Filter> >::data;

	DataFilterType(Filter const& filter, ObjectiveFunctionExpressionType<E> const& fun_type) : ObjectiveFunctionExpressionType<DataFilterType<E, Filter> >(fun_type.data), filter(filter), fun_type(fun_type) {

	}

	instance_type create_instance(sgl::DimConfig const& dim_config) const {
		return instance_type(fun_type, data, filter);
	}

	DataFilterType<E, Filter> operator()(Indices const& indices) const  {
			return DataFilterType<E, Filter>(filter(indices), fun_type(indices));
		}
};


template <typename E, typename Filter>
struct type_traits<DataFilterType<E, Filter> > {
	typedef typename E::data_type data_type;
	typedef DataFilter<typename E::instance_type, Filter> instance_type;
};

template<typename D>
class ObjectiveFunctionData {

protected:

	D const function_data;

public:

        ObjectiveFunctionData() : function_data() {}

	ObjectiveFunctionData(D const& data) : function_data(data) {}

        void set_data(D const& data) {
          const_cast<D&>(this->function_data) = data;
        }
};

template< typename E, typename D>
class ObjectiveFunctionType : public ObjectiveFunctionData<D>, public ObjectiveFunctionExpressionType<ObjectiveFunctionType<E, D > > {

private:
	ObjectiveFunctionType<E, D> const& operator=(ObjectiveFunctionType<E, D> const& s) {}

public:

	using ObjectiveFunctionExpressionType<ObjectiveFunctionType<E, D> >::data;

	typedef typename type_traits<ObjectiveFunctionType>::data_type data_type;
	typedef typename type_traits<ObjectiveFunctionType>::instance_type instance_type;

	ObjectiveFunctionType() : ObjectiveFunctionData<D>(), ObjectiveFunctionExpressionType<ObjectiveFunctionType<E, D> >(this->function_data) {}

	ObjectiveFunctionType(data_type const& data) : ObjectiveFunctionData<D>(data), ObjectiveFunctionExpressionType<ObjectiveFunctionType<E, D> >(this->function_data) {}

	ObjectiveFunctionType(ObjectiveFunctionType<E, D> const& s) : ObjectiveFunctionData<D>(s.function_data), ObjectiveFunctionExpressionType<ObjectiveFunctionType<E, D> >(this->function_data) {}

	instance_type create_instance(sgl::DimConfig const& dim_config) const {
		return instance_type(data, dim_config);
	}

	//TODO use constructor instead, this should save 2 constructor calls, hence 2 x data copy
	ObjectiveFunctionType<E, D> operator()(Indices const& indices) const  {
			return ObjectiveFunctionType<E, D>(data(indices));
	}

	//TODO do we need this ?
	template<typename Filter>
	DataFilterType<ObjectiveFunctionType<E, D>, Filter> operator [](Filter const& filter) const {
		return DataFilterType<ObjectiveFunctionType<E, D>, Filter>(filter, *this);
	}

};

template <typename E, typename D>
struct type_traits<ObjectiveFunctionType<E, D> > {
	typedef D data_type;
	typedef ObjectiveFunction<E> instance_type;
};


template<typename E>
class ObjectiveFunctionScalarMultiplicationType: public ObjectiveFunctionExpressionType<
		ObjectiveFunctionScalarMultiplicationType<E> > {

private:

	E const fun_type;
	sgl::numeric const lambda;

public:

	typedef typename type_traits<ObjectiveFunctionScalarMultiplicationType>::data_type data_type;
	typedef typename type_traits<ObjectiveFunctionScalarMultiplicationType>::instance_type instance_type;

	ObjectiveFunctionScalarMultiplicationType(sgl::numeric const& lambda, ObjectiveFunctionExpressionType<E> const& fun_type) : ObjectiveFunctionExpressionType<
			ObjectiveFunctionScalarMultiplicationType<E> >(fun_type.data), fun_type(fun_type), lambda(lambda) {
	}

	instance_type create_instance(sgl::DimConfig const& dim_config) const {
		return instance_type(lambda, fun_type.create_instance(dim_config));
	}

	ObjectiveFunctionScalarMultiplicationType<E> operator()(Indices const& indices) const  {
		return ObjectiveFunctionScalarMultiplicationType<E>(lambda, fun_type(indices));
	}

};


template <typename E>
struct type_traits<ObjectiveFunctionScalarMultiplicationType<E> > {
	typedef typename E::data_type data_type;
	typedef ObjectiveFunctionScalarMultiplication<typename E::instance_type> instance_type;
};

template<typename A, typename B>
class ObjectiveFunctionAdditionType : public ObjectiveFunctionExpressionType<ObjectiveFunctionAdditionType<A, B> > {

private:

	A const fun1_type;
	B const fun2_type;

public:

	typedef typename type_traits<ObjectiveFunctionAdditionType>::data_type data_type;
	typedef typename type_traits<ObjectiveFunctionAdditionType>::instance_type instance_type;


	ObjectiveFunctionAdditionType(ObjectiveFunctionExpressionType<A> const& fun1_type, ObjectiveFunctionExpressionType<B> const& fun2_type) : ObjectiveFunctionExpressionType<ObjectiveFunctionAdditionType<A, B> >(fun1_type.data), fun1_type(fun1_type), fun2_type(fun2_type) {
		//TODO assert that fun1 and fun2 have the same data
	}

	instance_type create_instance(sgl::DimConfig const& dim_config) const {
		return instance_type(fun1_type.create_instance(dim_config), fun2_type.create_instance(dim_config));
	}

	ObjectiveFunctionAdditionType<A, B> operator()(Indices const& indices) const  {
			return ObjectiveFunctionAdditionType<A, B>(fun1_type(indices), fun2_type(indices));
		}

};

template <typename A, typename B>
struct type_traits<ObjectiveFunctionAdditionType<A, B> > {
	typedef typename A::data_type data_type;
	typedef ObjectiveFunctionAddition<typename A::instance_type, typename B::instance_type> instance_type;
};

template<typename A, typename B>

ObjectiveFunctionAdditionType<A,B> operator+(ObjectiveFunctionExpressionType<A> const& fun1_type, ObjectiveFunctionExpressionType<B> const& fun2_type) {
  return ObjectiveFunctionAdditionType<A,B>(fun1_type, fun2_type);
}


#endif /* OBJECTIVEFUNCTIONEXPRESSIONTYPE_H_ */
