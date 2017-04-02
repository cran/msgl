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

//Uncomment to turn on debuging
//#undef NDEBUG

//Configuration
//Debugging
#ifndef NDEBUG
#define SGL_DEBUG
#endif

//Runtime checking for numerical problems
#define SGL_RUNTIME_CHECKS

//Check dimension of input objects
#define SGL_DIM_CHECKS

//Converges checks
#define SGL_CONVERGENCE_CHECK

//Exception handling
#define SGL_CATCH_EXCEPTIONS

// print information abt convergence
//#define SGL_DEBUG_INFO_QUADRATIC

//Should the timers be activated (only needed for profiling the code)
//#define DO_TIMING

// Show function entering/leaving
//#define FUNC_ENTER

//Sgl optimizer
#include <sgl.h>
#include "pkg_c_config.h"

/**********************************
 *
 *  msgl dense module
 *
 *********************************/

// Module name
#define MODULE_NAME msgl_dense

//Objective
#include "multinomial_loss.h"

#define OBJECTIVE multinomial

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#include "multinomial_response.h"
#define PREDICTOR sgl::LinearPredictor < sgl::matrix , MultinomialResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_subsampling.h>

/*********************************
 *
 *  msgl sparse module
 *
 *********************************/
// Reset macros
#undef MODULE_NAME
#undef OBJECTIVE
#undef PREDICTOR

// Module name
#define MODULE_NAME msgl_sparse

//Objective
#include "multinomial_loss.h"
#define OBJECTIVE multinomial_spx

#include <sgl/RInterface/sgl_lambda_seq.h>
#include <sgl/RInterface/sgl_fit.h>

#include "multinomial_response.h"
#define PREDICTOR sgl::LinearPredictor < sgl::sparse_matrix , MultinomialResponse >

#include <sgl/RInterface/sgl_predict.h>
#include <sgl/RInterface/sgl_subsampling.h>

/* **********************************
 *
 *  Registration of methods
 *
 ***********************************/

#include <R_ext/Rdynload.h>

static const R_CallMethodDef sglCallMethods[] = {

  SGL_LAMBDA(msgl_dense), SGL_LAMBDA(msgl_sparse),
  SGL_FIT(msgl_dense), SGL_FIT(msgl_sparse),
  SGL_PREDICT(msgl_dense), SGL_PREDICT(msgl_sparse),
  SGL_SUBSAMPLING(msgl_dense), SGL_SUBSAMPLING(msgl_sparse),

  {"r_pkg_c_config", (DL_FUNC) &r_pkg_c_config, 0},

  {NULL, NULL, 0}
};

extern "C" {
  void R_init_msgl(DllInfo *info);
}

void R_init_msgl(DllInfo *info)
{
  // Register the .Call routines.
  R_registerRoutines(info, NULL, sglCallMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
