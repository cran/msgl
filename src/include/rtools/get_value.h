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
type get_value(SEXP exp) {
	// no code should go here
	// Type unsupported => error here
	const type error_type_not_defined;
	&error_type_not_defined = 0;
}

template<>
double get_value(SEXP exp) {
	return static_cast<double>(*REAL(exp));
}

template<>
arma::u32 get_value(SEXP exp) {
	return static_cast<arma::u32>(*INTEGER(exp));
}

template<>
bool get_value(SEXP exp) {
	return static_cast<bool>(*LOGICAL(exp));
}

template<>
arma::Mat<double> get_value(SEXP exp) {

	double * ptr = REAL(exp);

	SEXP dim = Rf_getAttrib(exp, R_DimSymbol);

	unsigned int n_rows = INTEGER(dim)[0];
	unsigned int n_cols = INTEGER(dim)[1];

	return arma::conv_to< arma::Mat<double> >::from(arma::mat(ptr, n_rows, n_cols, false, true));
}

template<>
arma::Col<double> get_value(SEXP exp) {

	double *ptr = REAL(exp);

	return arma::conv_to< arma::Col<double> >::from(arma::vec(ptr, Rf_length(exp), false, true));
}

template<>
arma::Col<arma::u32> get_value(SEXP exp) {

	int *ptr = INTEGER(exp);

	return arma::conv_to< arma::Col<arma::u32> >::from(arma::Col<int>(ptr, Rf_length(exp), false, true));
}

template<>
arma::Col<arma::s32> get_value(SEXP exp) {

	int *ptr = INTEGER(exp);

	return arma::conv_to< arma::Col<arma::s32> >::from(arma::Col<int>(ptr, Rf_length(exp), false, true));
}

template<>
Indices get_value(SEXP exp) {
	return Indices(get_value<arma::Col<arma::u32> >(exp));
}

template<>
arma::sp_mat get_value(SEXP exp) {

	SEXP dim = VECTOR_ELT(exp, 0);
	unsigned int n_rows = INTEGER(dim)[0];
	unsigned int n_cols = INTEGER(dim)[1];

	SEXP col_ptrs = VECTOR_ELT(exp, 1);
	SEXP row_idx = VECTOR_ELT(exp, 2);
	SEXP values = VECTOR_ELT(exp, 3);

	unsigned int n_nonzero = Rf_length(values);

	arma::sp_mat m(n_rows, n_cols);

	if (n_nonzero == 0) {
		return m;
	}

	arma::uword* new_row_indices = arma::memory::acquire_chunked < uword > (n_nonzero + 1);
	double* new_values = arma::memory::acquire_chunked<double>(n_nonzero + 1);

	arma::arrayops::copy(new_values, REAL(values), n_nonzero);

	int * row_ptr = INTEGER(row_idx);
	for (unsigned int i = 0; i < n_nonzero; ++i) {
		new_row_indices[i] = static_cast<arma::uword>(row_ptr[i]);
	}

	new_row_indices[n_nonzero] = 0;

	int * col_ptr = INTEGER(col_ptrs);
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
arma::field<type> get_field(SEXP exp) {

	arma::field<type> res(static_cast<arma::u32>(Rf_length(exp)));

	for (arma::u32 i = 0; i < static_cast<arma::u32>(Rf_length(exp)); ++i) {
		SEXP elm = VECTOR_ELT(exp, i);
		res(i) = get_value<type>(elm);
	}

	return res;
}

#endif /* GET_VALUE_H_ */
