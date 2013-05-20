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

#ifndef TOOLS_H_
#define TOOLS_H_

//Return the index of the first element with maximal value
arma::u32
argmax(arma::vec x)
{

  arma::u32 j = 0;
  for (arma::u32 i = 1; i < x.n_elem; ++i)
    {
      if (x(i) > x(j))
        {
          j = i;
        }
    }

  return j;
}

std::string
create_error_msg(const char * msg, const char * file_name, int line_number)
{
  std::ostringstream error_msg;

  error_msg << msg << " (Assert failed in " << file_name << " at line " << line_number
      << " )";
  return error_msg.str();
}

#define ASSERT(condition, msg) if(!(condition)) throw std::runtime_error(create_error_msg(msg, __FILE__, __LINE__));

#endif /* TOOLS_H_ */
