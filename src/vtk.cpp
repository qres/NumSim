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

#include "vtk.hpp"
#include <cstring>
#include <cstdio>
//------------------------------------------------------------------------------
uint32_t VTK::_cnt = 0;
//------------------------------------------------------------------------------
VTK::VTK(const multi_real_t &h, const multi_index_t &size)
    : _h(h), _size(size) {
  _offset = multi_real_t(0.0);
  _handle = NULL;
}
//------------------------------------------------------------------------------
VTK::VTK(const multi_real_t &h, const multi_index_t &size,
         const multi_real_t &offset)
    : _h(h), _size(size), _offset(offset) {
  _handle = NULL;
}
//------------------------------------------------------------------------------
void VTK::Init(const char *path) {
  if (_handle)
    return;
  int flength = strlen(path) + 20;
  char *filename;
  filename = new char[flength];
  if (strlen(path))
    sprintf(filename, "%s_%i.vts", path, _cnt);
  else
    sprintf(filename, "%s_%i.vts", "field", _cnt);
  _handle = fopen(filename, "w");
  delete[] filename;

  fprintf(_handle, "<?xml version=\"1.0\"?>\n");
  fprintf(_handle, "<VTKFile type=\"StructuredGrid\">\n");
  fprintf(_handle, "<StructuredGrid WholeExtent=\"0 %i 0 %i 0 %i \">\n",
          _size[0] - 1, _size[1] - 1, (DIM == 3 ? _size[2] - 1 : 0));
  fprintf(_handle, "<Piece Extent=\"0 %i 0 %i 0 %i \">\n", _size[0] - 1,
          _size[1] - 1, (DIM == 3 ? _size[2] - 1 : 0));
  fprintf(_handle, "<Points>\n");
  fprintf(_handle, "<DataArray type=\"Float64\" format=\"ascii\" "
                   "NumberOfComponents=\"3\">\n");

  for (uint32_t z = 0; z < (DIM == 3 ? _size[2] : 1); ++z)
    for (uint32_t y = 0; y < _size[1]; ++y)
      for (uint32_t x = 0; x < _size[0]; ++x)
        fprintf(_handle, "%le %le %le\n", (double)x * _h[0] + _offset[0],
                (double)y * _h[1] + _offset[1],
                (DIM == 3 ? (double)z * _h[2] + _offset[2] : 0));

  fprintf(_handle, "</DataArray>\n");
  fprintf(_handle, "</Points>\n");
  fprintf(_handle, "<PointData>\n");
}
//------------------------------------------------------------------------------
void VTK::Finish() {
  if (!_handle)
    return;

  fprintf(_handle, "</PointData>\n");
  fprintf(_handle, "</Piece>\n");
  fprintf(_handle, "</StructuredGrid>\n");
  fprintf(_handle, "</VTKFile>\n");

  fclose(_handle);

  _handle = NULL;

  _cnt++;
}
//------------------------------------------------------------------------------
void VTK::AddScalar(const char *title, const Grid *grid) {
  if (!_handle)
    return;

  fprintf(_handle,
          "<DataArray Name=\"%s\" type=\"Float64\" format=\"ascii\">\n", title);

  multi_real_t pos;
  for (uint32_t z = 0; z < (DIM == 3 ? _size[2] : 1); ++z) {
#if DIM == 3
    pos[2] = (double)z * _h[2] + _offset[2];
#endif // DIM
    for (uint32_t y = 0; y < _size[1]; ++y) {
      pos[1] = (double)y * _h[1] + _offset[1];
      for (uint32_t x = 0; x < _size[0]; ++x) {
        pos[0] = (double)x * _h[0] + _offset[0];
        fprintf(_handle, "%le ", (double)grid->Interpolate(pos));
#if DIM == 3
        fprintf(_handle, "\n");
#endif // DIM
      }
#if DIM == 2
      fprintf(_handle, "\n");
#endif // DIM
    }
  }

  fprintf(_handle, "</DataArray>\n");
}
//------------------------------------------------------------------------------
void VTK::AddField(const char *title, const Grid *v1, const Grid *v2) {
  if (!_handle)
    return;

  fprintf(_handle, "<DataArray Name=\"%s\" type=\"Float64\" format=\"ascii\" "
                   "NumberOfComponents=\"3\">\n",
          title);

  multi_real_t pos;
#if DIM == 3
  pos[2] = 0;
#endif // DIM
  for (uint32_t y = 0; y < _size[1]; ++y) {
    pos[1] = (double)y * _h[1] + _offset[1];
    for (uint32_t x = 0; x < _size[0]; ++x) {
      pos[0] = (double)x * _h[0] + _offset[0];
      fprintf(_handle, "%le %le 0\n", (double)v1->Interpolate(pos),
              (double)v2->Interpolate(pos));
    }
  }

  fprintf(_handle, "</DataArray>\n");
}
//------------------------------------------------------------------------------
void VTK::AddField(const char *title, const Grid *v1, const Grid *v2,
                   const Grid *v3) {
  if (!_handle)
    return;

  fprintf(_handle, "<DataArray Name=\"%s\" type=\"Float64\" format=\"ascii\" "
                   "NumberOfComponents=\"3\">\n",
          title);

  multi_real_t pos;
#if DIM == 3
  pos[2] = 0;
#endif // DIM
  for (uint32_t y = 0; y < _size[1]; ++y) {
    pos[1] = (double)y * _h[1] + _offset[1];
    for (uint32_t x = 0; x < _size[0]; ++x) {
      pos[0] = (double)x * _h[0] + _offset[0];
      fprintf(_handle, "%le %le %le\n", (double)v1->Interpolate(pos),
              (double)v2->Interpolate(pos), (double)v3->Interpolate(pos));
    }
  }

  fprintf(_handle, "</DataArray>\n");
}
//------------------------------------------------------------------------------
