/*
 * Copyright (C) 2015   Malte Brunn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "typedef.hpp"
#include "grid.hpp"
#include <cstdio>
//------------------------------------------------------------------------------
#ifndef __VTK_HPP
#define __VTK_HPP
//------------------------------------------------------------------------------
#if DIM != 2 && DIM != 3
#error VTK only usable with DIM = {2,3}
#endif // DIM
//------------------------------------------------------------------------------
/// A class to generate VTK files to visualize grids with Paraview etc.
class VTK {
public:
  /// Constructs a VTK instance
  VTK(const multi_real_t &h, const multi_index_t &size);
  VTK(const multi_real_t &h, const multi_index_t &size,
      const multi_real_t &offset);

  /// Initializes the file
  void Init(const char *path);
  /// Closes the file
  void Finish();

  /// Add a field of scalar values
  void AddScalar(const char *title, const Grid *grid);
  /// Add a field of 2D data
  void AddField(const char *title, const Grid *v1, const Grid *v2);
  /// Add a field of 3D data
  void AddField(const char *title, const Grid *v1, const Grid *v2,
                const Grid *v3);

private:
  const multi_real_t &_h;
  const multi_index_t &_size;
  multi_real_t _offset;
  FILE *_handle;

  static uint32_t _cnt;
};
//------------------------------------------------------------------------------
/*!     \class VTK
 *      This class creates a VTK conform file to visualize 1D, 2D and/or 3D data
 *      with Paraview etc.
 */
/*!     \fn     VTK::VTK (const multi_real_t& h, const multi_index_t& size)
 *      \param h                The mesh width of the data grids to visualize
 *      \param size             The overall size of the data domain
 *
 *      Constructs an instance of the VTK class. It initializes the underlining
 *      domains grid points depending on \p h and \p size.
 */
/*!     \fn     VTK::VTK (const multi_real_t& h, const multi_index_t& size,
 * const multi_real_t& offset)
 *      \param h                The mesh width of the data grids to visualize
 *      \param size             The overall size of the data domain
 *      \param offset   The position of the first grid point
 *
 *      Constructs an instance of the VTK class. It initializes the underlining
 *      domains grid points depending on \p h and \p size. All grid nodes are
 *      shifted by \p offset.
 */
/*!     \fn void VTK::Init (const char* path)
 *      \param path     The path and filename of the VTK files.
 *
 *  Initializes the file header and writes the domain points. Each time this
 *      method is called a counter is incremented that is attached to the
 * filename.
 *      The path given with path is attached by the number and the ending
 * ".vts".
 *      If the path is left empty, the files will be named like "field_xxx.vts".
 */
/*!     \fn void AddScalar (const char* title, const Grid* grid)
 *      \param title    The name of the field in the VTK file
 *      \param grid             The grid with the scalar values.
 *
 *      Writes a field of 1D values to the VTK file.
 */
/*!     \fn void AddField (const char* title, const Grid* v1, const Grid* v2)
 *      \param title    The name of the field in the VTK file
 *      \param v1               The first component of the 2D field
 *      \param v2               The second component of the 2D field
 *
 *      Writes a field of 2D values to the VTK file. The 2D field consists of
 *      a first and a second component.
 */
/*!     \fn void AddField (const char* title, const Grid* v1, const Grid* v2,
 * const Grid* v3)
 *      \param title    The name of the field in the VTK file
 *      \param v1               The first component of the 3D field
 *      \param v2               The second component of the 3D field
 *      \param v3               The third component of the 3D field
 *
 *      Writes a field of 3D values to the VTK file. The 3D field consists of
 *      a first, second and a third component.
 */
//------------------------------------------------------------------------------
#endif // __VTK_HPP
