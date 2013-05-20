/*
	Lightweight tools for R and c++ integration.
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

#ifndef GET_VALUE_H_
#define GET_VALUE_H_

template<typename type>
type get_value(R::SEXP exp) {
	// no code should go here
	// Type unsupported => error here
	const type error_type_not_defined;
	&error_type_not_defined = 0;
}

template<>
double get_value(R::SEXP exp) {
	return static_cast<double>(*R::REAL(exp));
}

template<>
arma::u32 get_value(R::SEXP exp) {
	return static_cast<arma::u32>(*R::INTEGER(exp));
}

template<>
bool get_value(R::SEXP exp) {
	return static_cast<bool>(*R::LOGICAL(exp));
}

template<>
arma::Mat<double> get_value(R::SEXP exp) {

	double * ptr = R::REAL(exp);

	R::SEXP dim = R::getAttrib(exp, R::R_DimSymbol);

	unsigned int n_rows = R::INTEGER(dim)[0];
	unsigned int n_cols = R::INTEGER(dim)[1];

	return arma::conv_to< arma::Mat<double> >::from(arma::mat(ptr, n_rows, n_cols, false, true));
}

template<>
arma::Col<double> get_value(R::SEXP exp) {

	double *ptr = R::REAL(exp);

	return arma::conv_to< arma::Col<double> >::from(arma::vec(ptr, R::Rf_length(exp), false, true));
}

template<>
arma::Col<arma::u32> get_value(R::SEXP exp) {

	int *ptr = R::INTEGER(exp);

	return arma::conv_to< arma::Col<arma::u32> >::from(arma::Col<int>(ptr, R::Rf_length(exp), false, true));
}

template<>
arma::Col<arma::s32> get_value(R::SEXP exp) {

	int *ptr = R::INTEGER(exp);

	return arma::conv_to< arma::Col<arma::s32> >::from(arma::Col<int>(ptr, R::Rf_length(exp), false, true));
}

template<>
Indices get_value(R::SEXP exp) {
	return Indices(get_value<arma::Col<arma::u32> >(exp));
}

template<>
arma::sp_mat get_value(R::SEXP exp) {

	R::SEXP dim = R::VECTOR_ELT(exp, 0);
	unsigned int n_rows = R::INTEGER(dim)[0];
	unsigned int n_cols = R::INTEGER(dim)[1];

	R::SEXP col_ptrs = R::VECTOR_ELT(exp, 1);
	R::SEXP row_idx = R::VECTOR_ELT(exp, 2);
	R::SEXP values = R::VECTOR_ELT(exp, 3);

	unsigned int n_nonzero = R::length(values);

	arma::sp_mat m(n_rows, n_cols);

	if (n_nonzero == 0) {
		return m;
	}

	arma::uword* new_row_indices = arma::memory::acquire_chunked < uword > (n_nonzero + 1);
	double* new_values = arma::memory::acquire_chunked<double>(n_nonzero + 1);

	arma::arrayops::copy(new_values, R::REAL(values), n_nonzero);

	int * row_ptr = R::INTEGER(row_idx);
	for (unsigned int i = 0; i < n_nonzero; ++i) {
		new_row_indices[i] = static_cast<arma::uword>(row_ptr[i]);
	}

	new_row_indices[n_nonzero] = 0;

	int * col_ptr = R::INTEGER(col_ptrs);
	for (unsigned int i = 0; i < n_cols + 2; ++i) {
		arma::access::rwp(m.col_ptrs)[i] = static_cast<arma::uword>(col_ptr[i]);
	}

	arma::memory::release(m.values);
	arma::memory::release(m.row_indices);

	arma::access::rw(m.values) = new_values;
	arma::access::rw(m.row_indices) = new_row_indices;

	// Update counts and such.
	arma::access::rw(m.n_nonzero) = n_nonzero;

	return m;
}

template<typename type>
arma::field<type> get_field(R::SEXP exp) {

	arma::field<type> res(static_cast<arma::u32>(R::length(exp)));

	for (arma::u32 i = 0; i < static_cast<arma::u32>(R::length(exp)); ++i) {
		R::SEXP elm = R::VECTOR_ELT(exp, i);
		res(i) = get_value<type>(elm);
	}

	return res;
}

#endif /* GET_VALUE_H_ */
