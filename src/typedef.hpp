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

#include <stdint.h>

//------------------------------------------------------------------------------

#ifndef __TYPEDEF_HPP
#define __TYPEDEF_HPP

//------------------------------------------------------------------------------

#define DIM 2
#define REAL_TYPE double
#define MPI_REAL_T MPI_DOUBLE
#define INDEX_TYPE uint32_t

//------------------------------------------------------------------------------

/// Typedef for reals
typedef REAL_TYPE real_t;

/// Typedef for integers
typedef INDEX_TYPE index_t;

//------------------------------------------------------------------------------

// we only want these in the nvcc compiler
#ifndef NVCC
#define __host__
#define __device__
#endif

/// Template for array/vector types
template <typename _type, uint32_t _dim>
struct array_t {
  // Constructors
  // Empty; initialize array_t to 0
  __host__ __device__
  array_t() {
    for (uint32_t i = 0; i < _dim; ++i)
      x[i] = 0;
  }

  // Value; initialize array_to to val
  __host__ __device__
  explicit array_t(const _type &v1) {
    for (uint32_t i = 0; i < _dim; ++i)
      x[i] = v1;
  }

  // Value 2: initialize first value to v1, all others to v2
  __host__ __device__
  array_t(const _type &v1, const _type &v2) {
    x[0] = v1;
    for (uint32_t i = 1; i < _dim; ++i)
      x[i] = v2;
  }

  // Copy-constructor: field of values
  __host__ __device__
  array_t(const _type(&cp)[_dim]) {
    for (uint32_t i = 0; i < _dim; ++i)
      x[i] = cp[i];
  }
  // Copy-constructor: array_t
  __host__ __device__
  array_t(const array_t<_type, _dim> &cp) {
    for (uint32_t i = 0; i < _dim; ++i)
      x[i] = cp.x[i];
  }

  // access operator
  __host__ __device__
  _type &operator[](uint32_t i) { return x[i]; }
  __host__ __device__
  const _type &operator[](uint32_t i) const { return x[i]; }

  // store values and size
  _type x[_dim];
  static const uint32_t dim = _dim;
};
/// Typedef for d-dimensional array of reals
typedef array_t<real_t, DIM> multi_real_t;
/// Typedef for d-dimensional array of integer
typedef array_t<index_t, DIM> multi_index_t;

//------------------------------------------------------------------------------

/// Forward declaration of classes used
class Iterator;
class Communicator;
class Geometry;
class Grid;
//class FlagGrid;
class Parameter;
class Solver;
class RedOrBlackSOR;
class MultiGrid;
class Cfg;
class Compute;
class ParticleLine;
class StreakLine;
class PathLine;

#endif // __TYPEDEF_HPP
