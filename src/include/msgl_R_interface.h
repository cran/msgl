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

#include <iostream>
#include <cstdio>

//TODO fix openmp problem with rstream
//TODO until fixed use write_msg in openmp loops

class rstream : public std::streambuf {
public:
protected:
  virtual std::streamsize xsputn(const char *s, std::streamsize n);
  virtual int overflow(int c = EOF);
  virtual int sync();
};

class rostream : public std::ostream
{
protected:
    rstream buf;

public:
    rostream() : std::ostream(&buf) {}
};

rostream rout;

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
