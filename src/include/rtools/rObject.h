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


#ifndef ROBJECT_H_
#define ROBJECT_H_

class rObject {

private:

	SEXP exp;

	arma::u32 number_of_protects;

	bool const unprotect_on_destruction;

public:

	rObject(arma::u32 value, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {
		PROTECT(exp = Rf_allocVector(INTSXP, 1));
		INTEGER(exp)[0] = value;
	}

	rObject(double value, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {
		PROTECT(exp = Rf_allocVector(REALSXP, 1));
		REAL(exp)[0] = value;
	}

	rObject(arma::Mat<double> const& m, bool unprotect_on_destruction = true) :
			number_of_protects(2), unprotect_on_destruction(
					unprotect_on_destruction) {

		SEXP matrixDim;
		PROTECT(matrixDim = Rf_allocVector(INTSXP, 2));
		INTEGER(matrixDim)[0] = m.n_rows;
		INTEGER(matrixDim)[1] = m.n_cols;

		PROTECT(exp = Rf_allocVector(REALSXP, m.n_elem));

		//Copy data
		arrayops::copy(REAL(exp), m.mem, m.n_elem);

		Rf_setAttrib(exp, R_DimSymbol, matrixDim);
	}

	rObject(arma::Mat<arma::u32> const& m, bool unprotect_on_destruction = true) :
			number_of_protects(2), unprotect_on_destruction(
					unprotect_on_destruction) {

		SEXP matrixDim;
		PROTECT(matrixDim = Rf_allocVector(INTSXP, 2));
		INTEGER(matrixDim)[0] = m.n_rows;
		INTEGER(matrixDim)[1] = m.n_cols;

		PROTECT(exp = Rf_allocVector(INTSXP, m.n_rows * m.n_cols));

		//Copy data
		copy_cast(INTEGER(exp), m.mem, m.n_elem);

		Rf_setAttrib(exp, R_DimSymbol, matrixDim);
	}

	rObject(Col<double> const& v, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {

		PROTECT(exp = Rf_allocVector(REALSXP, v.n_elem));

		//Copy data
		arrayops::copy(REAL(exp), v.mem, v.n_elem);

	}

	rObject(Col<arma::u32> const& v, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {

		PROTECT(exp = Rf_allocVector(INTSXP, v.n_elem));

		//Copy data
		copy_cast(INTEGER(exp), v.mem, v.n_elem);
	}

	rObject(Indices i, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {

		PROTECT(exp = Rf_allocVector(INTSXP, i.size()));

		//Copy data
		arma::uvec tmp = i.getElements();
		copy_cast(INTEGER(exp), tmp.mem, tmp.n_elem);

	}

	rObject(arma::sp_mat const& m, bool unprotect_on_destruction = true) :
			number_of_protects(0), unprotect_on_destruction(
					unprotect_on_destruction) {
		exp = create_sparse_matrix_object(m, number_of_protects);
	}

	SEXP create_sparse_matrix_object(arma::sp_mat const& spMat,
			arma::u32 & number_of_protects) {

		number_of_protects += 5;

		SEXP sexp_object;
		PROTECT(sexp_object = Rf_allocVector(VECSXP, 4)); // Creating a list with 4 elements

		//TODO names on list

		SEXP dim;
		PROTECT(dim = Rf_allocVector(INTSXP, 2));
		SET_VECTOR_ELT(sexp_object, 0, dim);
		INTEGER(dim)[0] = spMat.n_rows;
		INTEGER(dim)[1] = spMat.n_cols;

		SEXP col_ptrs;
		PROTECT(col_ptrs = Rf_allocVector(INTSXP, spMat.n_cols+1));
		SET_VECTOR_ELT(sexp_object, 1, col_ptrs);
		copy_cast(INTEGER(col_ptrs), spMat.col_ptrs, spMat.n_cols+1);

		SEXP row_indices;
		PROTECT(row_indices = Rf_allocVector(INTSXP, spMat.n_nonzero));
		SET_VECTOR_ELT(sexp_object, 2, row_indices);
		copy_cast(INTEGER(row_indices), spMat.row_indices, spMat.n_nonzero);

		SEXP values;
		PROTECT(values = Rf_allocVector(REALSXP, spMat.n_nonzero));
		SET_VECTOR_ELT(sexp_object, 3, values);
        arrayops::copy(REAL(values), spMat.values, spMat.n_nonzero);

		return sexp_object;
	}



	template<typename T>
	rObject(arma::field<T> const& field, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {

		PROTECT(exp = Rf_allocVector(VECSXP, field.n_elem)); // Creating a list with n_elem elements

		//Construct list
		unsigned int i;
		for (i = 0; i < field.n_elem; i++) {
			// attaching
			rObject tmp(field(i), true);
			//number_of_protects += tmp.n_protects();
			SET_VECTOR_ELT(exp, i, tmp);
		}

	}

	~rObject() {
		if (unprotect_on_destruction) {
			UNPROTECT(number_of_protects);
		}
	}

	operator SEXP() const {
		return getSEXP();
	}

	SEXP getSEXP() const {
		return exp;
	}

	arma::u32 n_protects() const {
		return number_of_protects;
	}

};

#endif /* ROBJECT_H_ */
