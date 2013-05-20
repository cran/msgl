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


#ifndef RSTREAM_H_
#define RSTREAM_H_

std::streamsize 
rstream::xsputn(const char *s, std::streamsize n)
{
	R::Rprintf( "%.*s", n, s);

    return n;
}

int 
rstream::overflow(int c)
{
	 if (c != EOF) {
		R::Rprintf( "%.1s", &c);
	 }

	 return c;
}

int
rstream::sync(){
	R::R_FlushConsole();
    return 0;
}

#endif /* RSTREAM_H_ */
