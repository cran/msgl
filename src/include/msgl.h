/* Routines for multinomial and logistic sparse group lasso regression.
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

#ifndef MSGL_H_
#define MSGL_H_

namespace msgl {

#include "msgl/msgl_matrix_data.h"

#include "msgl/msgl_multinomial_response.h"
#include "msgl/msgl_multinomial_predictor.h"
#include "msgl/msgl_mg_loss.h"

#ifdef MSGL_EXTENSIONS
#endif

}

#endif /* MSGL_H_ */
