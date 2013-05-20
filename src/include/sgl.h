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

#ifndef SGL_H_
#define SGL_H_

#ifdef SGL_USE_OPENMP
#include <omp.h>
#endif

#include <limits>
#include <time.h>
#include <stdexcept>

#include <boost/tuple/tuple.hpp>
using boost::tuple;

#include <armadillo>
#include "sgl/arma_additions.h"

#include "sgl/simple_timer.h"
#include "sgl/tools.h"

#include "Indices/Indices.h"
#include "Indices/GroupedIndices.h"

namespace sgl {
#include "sgl/numeric.h"
#include "sgl/config.h"
#include "sgl/AlgorithmConfigurationDefault.h"
#include "sgl/DimConfig.h"
#include "sgl/BlockVector.h"
#include "sgl/SglProblem.h"
#include "sgl/SglOptimizer.h"
#include "sgl/ObjectiveFunction.h"
#include "sgl/ObjectiveFunctionExpressionType.h"
#include "sgl/interface_basic.h"

#ifdef SGL_EXTENSIONS
//#warning compling with sgl extensions
#endif

}

#endif /* SGL_H_ */
