/*
GeNuSys - computations with generalized number systems
Copyright (C) 2015-2017  Bence Németh
Copyright (C) 2017  Tamás Krutki

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GENUSYS_UTILS_H_
#define GENUSYS_UTILS_H_

#define GENUSYS_VERSION_MAJOR 1
#define GENUSYS_VERSION_MINOR 1

#ifndef NDEBUG
#define ASSERT_EXCEPTION(x,y) if(!(x)) throw y{#x}
#else
#define ASSERT_EXCEPTION(x,y)
#endif

#endif // GENUSYS_UTILS_H_
