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
#include "flaggrid.hpp"
//------------------------------------------------------------------------------
#ifndef __GEOMETRY_HPP
#define __GEOMETRY_HPP
//------------------------------------------------------------------------------
class Geometry {
public:
  /// Constructs a default geometry:
  // driven cavity with 128 x 128 grid, no-slip boundary conditions
  // as shown below
  //
  //      u=1, v=0
  //    -------------
  //    |           |
  // u=0|           |u=0
  // v=0|           |v=0
  //    |           |
  //    |           |
  //    -------------
  //      u=0, v=0
  Geometry();
  Geometry(const Communicator *comm);

  ~Geometry();

  /// Loads a geometry from a file
  void Load(const char *file);

  /// reads the thread dim fromm communicator and sets the local grid sizes
  void split_for_comm();

  /// Returns the number of cells in each dimension
  const multi_index_t &Size() const;
  /// Returns the total number of cells in each dimension
  const multi_index_t &TotalSize() const;
  /// Returns the length of the domain
  const multi_real_t &Length() const;
  /// Returns the total length of the domain
  const multi_real_t &TotalLength() const;
  /// Returns the meshwidth
  const multi_real_t &Mesh() const;

  const multi_real_t &Velocity() const;
  const real_t &Pressure() const;

  const FlagGrid &Flags() const;
  FlagGrid &Flags();

  /// Updates the velocity field u boundaries
  void Update_U(Grid *u) const;
  /// Updates the velocity field v boundaries
  void Update_V(Grid *v) const;
  /// Updates the pressure field p boundaries
  void Update_P(Grid *p) const;

private:
  const Communicator *_comm;

  multi_index_t _size;
  multi_index_t _bsize; // block size
  multi_real_t _length;
  multi_real_t _blength; // block size
  multi_real_t _h;

  multi_real_t _velocity;
  real_t _pressure;

  FlagGrid *_flags;
};
//------------------------------------------------------------------------------
#endif // __GEOMETRY_HPP
