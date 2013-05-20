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

#ifndef MSGL_R_INTERFACE_H_
#define MSGL_R_INTERFACE_H_

//Progress monitor
#include <progress.hpp>

#include <RcppCommon.h>
#include <Rconfig.h>
#include <RcppArmadilloConfig.h>

// Debugging
#ifdef SGL_DEBUG
// Do debugging
#ifdef ARMA_NO_DEBUG
#undef ARMA_NO_DEBUG
#endif
#ifdef NDEBUG
#undef NDEBUG
#endif
//#define SGL_DEBUG_SIMPLE
//#define SGL_DEBUG_COMPLEX
//#define SGL_DEBUG_INFO_ALL
#define PRINT_BACKTRANCE
//#define SGL_DEBUG_INFO_STEPSIZE
#else
// Do no debugging
#define ARMA_NO_DEBUG
#define NDEBUG
#endif

#ifdef PRINT_BACKTRANCE
#include <Backtrace.h>
#endif

#include <armadillo>
#include <Rcpp.h>

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
using boost::tuple;

#include "sgl/arma_additions.h"

#include <sgl.h>
#include <rtools.h>
#include <msgl.h>

namespace msgl {
#include "msgl/msgl_algorithm_config.h"
#include "msgl/msgl_mg.h"
#include "msgl/msgl_multinomial_weighted_sparse.h"
}

#ifdef MSGL_EXTENSIONS
#include "msgl/Extensions/msgl_debug_tools.h"
#endif


#endif /* MSGL_R_INTERFACE_H_ */
