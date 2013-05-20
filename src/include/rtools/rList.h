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


#ifndef RLIST_H_
#define RLIST_H_


class rList {

private:

	R::SEXP listEXP;
	R::SEXP names;

	unsigned int number_of_elements;
	unsigned int index;

	unsigned int const number_of_protects;

public:

	rList(unsigned int number_of_elements) :
		number_of_elements(number_of_elements), index(0), number_of_protects(2) {
		R::PROTECT(listEXP = R::allocVector(VECSXP, number_of_elements));
		R::PROTECT(names = R::allocVector(VECSXP, number_of_elements));
	}

	rList(R::SEXP list) :
		listEXP(list), names(getAttrib(listEXP, R::R_NamesSymbol)), number_of_elements(length(list)),
				index(number_of_elements), number_of_protects(0) {
	}

	~rList() {
		if (number_of_protects > 0) {
			R::UNPROTECT(number_of_protects);
		}
	}

	void attach(R::SEXP element, string const& name) {

		if (index >= number_of_elements) {
			throw runtime_error("Internal error - elements in r list exceed max number of elements.");
		}

		R::SET_VECTOR_ELT(listEXP, index, element);
		R::SET_VECTOR_ELT(names, index, R::mkChar(name.c_str()));
		++index;
	}

	R::SEXP get(unsigned int index) const {
		return VECTOR_ELT(listEXP, index);
	}

	int getIndex(string const& name) const {

		for (u32 index = 0; index < number_of_elements; ++index) {
			if (name.compare(R::CHAR(R::STRING_ELT(names, index))) == 0) {
				return index;
			}
		}

		return -1;

	}

	operator R::SEXP() const {
		return getSEXP();
	}

	R::SEXP getSEXP() const {
		setAttrib(listEXP, R::R_NamesSymbol, names);
		return listEXP;
	}
};


#endif /* RLIST_H_ */
