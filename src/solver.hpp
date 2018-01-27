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
//------------------------------------------------------------------------------
#include "typedef.hpp"
#include <iostream>
//------------------------------------------------------------------------------
#ifndef __SOLVER_HPP
#define __SOLVER_HPP
//------------------------------------------------------------------------------

#define assert(val) if (!val) { \
    std::cout << "ERR: ASSERT FAILED at " << __FILE__ << ':' << __LINE__ << std::endl; \
}

/** abstract base class for an iterative solver
*/
class Solver {
public:
  /// Constructor of the abstract Solver class
  Solver(const Geometry *geom);
  /// Destructor of the Solver Class
  virtual ~Solver();

  /// This function must be implemented in a child class
  // @param [in][out] grid current values
  // @param [in]      rhs  right hand side values
  // @returns accumulated residual
  virtual real_t Cycle(Grid *grid, const Grid *rhs) const = 0;

protected:
  const Geometry *_geom;

  /// Returns the residual at [it] for the pressure-Poisson equation
  real_t localRes(const Iterator &it, const Grid *grid, const Grid *rhs) const;
};

//------------------------------------------------------------------------------

/** concrete SOR solver
*/
class SOR : public Solver {
public:
  /// Constructs an actual SOR solver
  SOR(const Geometry *geom, const real_t &omega);
  /// Destructor
  ~SOR();

  /// Returns the total residual and executes a solver cycle
  // @param grid current pressure values
  // @param rhs right hand side
  real_t Cycle(Grid *grid, const Grid *rhs) const;

protected:
  real_t _omega;
};
//------------------------------------------------------------------------------

/** concrete Red or Balck SOR solver
 */
class RedOrBlackSOR : public SOR {
public:
  /// Constructs an actual SOR solver
  RedOrBlackSOR(const Geometry *geom, const real_t &omega);
  /// Destructor
  ~RedOrBlackSOR();

  real_t RedCycle(Grid *grid, const Grid *rhs) const;
  real_t BlackCycle(Grid *grid, const Grid *rhs) const;
};
//------------------------------------------------------------------------------

// to configure the inner most jacobi solver on the coarsest grid
class Cfg_jacobi {
public:
    Cfg_jacobi();
    Cfg_jacobi(unsigned int max_iters, double max_res);

    unsigned int max_iters;
    double       max_res;
};

class Cfg_multigrid {
public:
    Cfg_multigrid();

    unsigned int max_iters;
    double       max_res;
    unsigned int max_iters_jacobi_pre;
    double       max_res_jacobi_pre;
    unsigned int max_iters_jacobi_post;
    double       max_res_jacobi_post;
};

// to configure multigrid methods
class Cfg {
public:
    Cfg();
    ~Cfg();

    static Cfg v_cycle(unsigned int max_iters, double max_res, unsigned int levels, unsigned int max_iters_pre_post, Cfg_jacobi inner = Cfg_jacobi(5, 0));
    static Cfg two_grid(unsigned int max_iters, double max_res, unsigned int max_iters_pre_post, Cfg_jacobi inner = Cfg_jacobi());
    static Cfg jacobi(unsigned int max_iters, double max_res);

    unsigned int num_levels;
    Cfg_multigrid* levels;
    bool count_iters_like_single_grid;
};

#define USE_CUDA
//#define MIXED_PRECISION

#ifdef MIXED_PRECISION
#define solver_real_t float
#else
#define solver_real_t real_t
#endif

class MultiGrid : public Solver {
public:
    MultiGrid(const Geometry *geom, const Cfg* cfg);
    ~MultiGrid();
    // one MG cycle
    real_t Cycle(Grid *grid, const Grid *rhs) const;

private:
    template<typename F, typename T>
    unsigned int solve_mg_flat(T* _u0, const T* _b) const;

    const Cfg* _cfg;

    char** flags;

    // these are ignored if we don't use cuda
    mutable real_t* d_u0_f64;
    mutable real_t* d_b_f64;
    mutable solver_real_t* d_res_f32;
    mutable solver_real_t* d_err0_f32;

    // this contains all the data if we only use a single allcoation
    mutable solver_real_t* buffer;

    mutable solver_real_t** u0;
    mutable solver_real_t* scratch;
    mutable solver_real_t** b;
    mutable solver_real_t** r_0;
    mutable solver_real_t** e_0;
    mutable solver_real_t** res0;
    mutable solver_real_t** res1;
    mutable multi_index_t* N;
#ifdef MIXED_PRECISION
    mutable solver_real_t* res_f32;
    mutable solver_real_t* err_f32;
#endif
};

#endif // __SOLVER_HPP
