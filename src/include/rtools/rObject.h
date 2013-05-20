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

	R::SEXP exp;

	arma::u32 number_of_protects;

	bool const unprotect_on_destruction;

public:

	rObject(arma::u32 value, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {
		R::PROTECT(exp = R::allocVector(INTSXP, 1));
		R::INTEGER(exp)[0] = value;
	}

	rObject(double value, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {
		R::PROTECT(exp = R::allocVector(REALSXP, 1));
		R::REAL(exp)[0] = value;
	}

	rObject(arma::Mat<double> const& m, bool unprotect_on_destruction = true) :
			number_of_protects(2), unprotect_on_destruction(
					unprotect_on_destruction) {

		R::SEXP matrixDim;
		R::PROTECT(matrixDim = R::allocVector(INTSXP, 2));
		R::INTEGER(matrixDim)[0] = m.n_rows;
		R::INTEGER(matrixDim)[1] = m.n_cols;

		R::PROTECT(exp = R::allocVector(REALSXP, m.n_elem));

		//Copy data
		arrayops::copy(R::REAL(exp), m.mem, m.n_elem);

		setAttrib(exp, R::R_DimSymbol, matrixDim);
	}

	rObject(arma::Mat<arma::u32> const& m, bool unprotect_on_destruction = true) :
			number_of_protects(2), unprotect_on_destruction(
					unprotect_on_destruction) {

		R::SEXP matrixDim;
		R::PROTECT(matrixDim = R::allocVector(INTSXP, 2));
		R::INTEGER(matrixDim)[0] = m.n_rows;
		R::INTEGER(matrixDim)[1] = m.n_cols;

		R::PROTECT(exp = R::allocVector(INTSXP, m.n_rows * m.n_cols));

		//Copy data
		copy_cast(R::INTEGER(exp), m.mem, m.n_elem);

		setAttrib(exp, R::R_DimSymbol, matrixDim);
	}

	rObject(Col<double> const& v, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {

		R::PROTECT(exp = R::allocVector(REALSXP, v.n_elem));

		//Copy data
		arrayops::copy(R::REAL(exp), v.mem, v.n_elem);

	}

	rObject(Col<arma::u32> const& v, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {

		R::PROTECT(exp = R::allocVector(INTSXP, v.n_elem));

		//Copy data
		copy_cast(R::INTEGER(exp), v.mem, v.n_elem);
	}

	rObject(Indices i, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {

		R::PROTECT(exp = R::allocVector(INTSXP, i.size()));

		//Copy data
		arma::uvec tmp = i.getElements();
		copy_cast(R::INTEGER(exp), tmp.mem, tmp.n_elem);

	}

	rObject(arma::sp_mat const& m, bool unprotect_on_destruction = true) :
			number_of_protects(0), unprotect_on_destruction(
					unprotect_on_destruction) {
		exp = create_sparse_matrix_object(m, number_of_protects);
	}

	R::SEXP create_sparse_matrix_object(arma::sp_mat const& spMat,
			arma::u32 & number_of_protects) {

		number_of_protects += 5;

		R::SEXP sexp_object;
		R::PROTECT(sexp_object = R::allocVector(VECSXP, 4)); // Creating a list with 4 elements

		//TODO names on list

		R::SEXP dim;
		R::PROTECT(dim = R::allocVector(INTSXP, 2));
		R::SET_VECTOR_ELT(sexp_object, 0, dim);
		R::INTEGER(dim)[0] = spMat.n_rows;
		R::INTEGER(dim)[1] = spMat.n_cols;

		R::SEXP col_ptrs;
		R::PROTECT(col_ptrs = R::allocVector(INTSXP, spMat.n_cols+1));
		R::SET_VECTOR_ELT(sexp_object, 1, col_ptrs);
		copy_cast(INTEGER(col_ptrs), spMat.col_ptrs, spMat.n_cols+1);

		R::SEXP row_indices;
		R::PROTECT(row_indices = R::allocVector(INTSXP, spMat.n_nonzero));
		R::SET_VECTOR_ELT(sexp_object, 2, row_indices);
		copy_cast(INTEGER(row_indices), spMat.row_indices, spMat.n_nonzero);

		R::SEXP values;
		R::PROTECT(values = R::allocVector(REALSXP, spMat.n_nonzero));
		R::SET_VECTOR_ELT(sexp_object, 3, values);
        arrayops::copy(REAL(values), spMat.values, spMat.n_nonzero);

		return sexp_object;
	}



	template<typename T>
	rObject(arma::field<T> const& field, bool unprotect_on_destruction = true) :
			number_of_protects(1), unprotect_on_destruction(
					unprotect_on_destruction) {

		R::PROTECT(exp = R::allocVector(VECSXP, field.n_elem)); // Creating a list with n_elem elements

		//Construct list
		unsigned int i;
		for (i = 0; i < field.n_elem; i++) {
			// attaching
			rObject tmp(field(i), true);
			//number_of_protects += tmp.n_protects();
			R::SET_VECTOR_ELT(exp, i, tmp);
		}

	}

	~rObject() {
		if (unprotect_on_destruction) {
			R::UNPROTECT(number_of_protects);
		}
	}

	operator R::SEXP() const {
		return getSEXP();
	}

	R::SEXP getSEXP() const {
		return exp;
	}

	arma::u32 n_protects() const {
		return number_of_protects;
	}

};

#endif /* ROBJECT_H_ */
