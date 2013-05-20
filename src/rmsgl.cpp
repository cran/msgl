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

//Configuration

//Debugging
#define ARMA_NO_DEBUG
#define SGL_DIM_CHECKS
//#define SGL_DEBUG_SIMPLE
//#define SGL_DEBUG_COMPLEX
//#define SGL_DEBUG_INFO_ALL
//#define DEBUG_BACKTRACE

//#define SGL_DEBUG_INFO_STEPSIZE
//#define SGL_DEBUG_INFO_QUADRATIC

//Runtime checking for numerical problems
#define SGL_RUNTIME_CHECKS

//Converges checks
#define SGL_CONVERGENCE_CHECK

//Exception handling
#define SGL_CATCH_EXCEPTIONS

//Should the timers be activated (only needed for profiling the code)
//#define SGL_TIMING

//Should openmp be used
#ifndef _OPENMP
//No openmp
//TODO a warning should be displayed when the package is loaded in R
//openmp (multithreading) not supported on this system - compiling without openmp support
#else
//Use openmp
#define SGL_USE_OPENMP
#endif

//Compile with extensions if present
//#define SGL_EXTENSIONS
//#define MSGL_EXTENSIONS

void report_error(const char *msg);

#include "include/msgl_R_interface.h"
#include <memory>

void report_error(const char *msg) {
	R::Rf_error(msg);
}
