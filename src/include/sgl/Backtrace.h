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

#ifndef BACKTRACE_H_
#define BACKTRACE_H_

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <exception>
#include <cxxabi.h>
#include <string.h>

class Backtrace
{
private:
  void *array[10];
  size_t size;

  std::string
  demangle(const char* symbol) const
  {

    int status;
    char* demangled;

    //allocate mem

    char *tmp = (char *) malloc(strlen(symbol) * sizeof(char));

    if (tmp == NULL)
      {
        return symbol;
      }

    //first, try to demangle a c++ name

    if (1 == sscanf(symbol, "%*[^(]%*[^_]%[^)+]", tmp))
      {
        demangled = abi::__cxa_demangle(tmp, NULL, NULL, &status);
        if (status == 0)
          {
            std::string result(demangled);
            free(demangled);
            return result;
          }
      }

    return symbol;
  }

public:

  Backtrace()
  {
    size = backtrace(array, 10);
  }

  void
  print_trace() const
  {
    // print out all the frames to stderr

                char** symbols = backtrace_symbols(array, size);

                for (unsigned int i = 1; i < size; ++i) {
                        printf("%i) %s\n\n", i, demangle(symbols[i]).c_str());
                }

                free(symbols);
  }
};

#endif /* BACKTRACE_H_ */
