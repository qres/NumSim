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
struct Cfg_jacobi {
    unsigned int max_iters;
    double       max_res;

    Cfg_jacobi() :
        max_iters(100000000),
        max_res(1e-12) {
        //all done
    }

    Cfg_jacobi(unsigned int max_iters, double max_res) :
        max_iters(max_iters),
        max_res(max_res) {
        //all done
    }
};

struct Cfg_multigrid {
    unsigned int max_iters;
    double       max_res;
    unsigned int max_iters_jacobi_pre;
    double       max_res_jacobi_pre;
    unsigned int max_iters_jacobi_post;
    double       max_res_jacobi_post;

    Cfg_multigrid() :
        max_iters(100000000),
        max_res(1e-12),
        max_iters_jacobi_pre(10),
        max_res_jacobi_pre(1e-12),
        max_iters_jacobi_post(10),
        max_res_jacobi_post(1e-12) {
        //all done
    }
};

// to configure multigrid methods
struct Cfg {
    unsigned int num_levels;
    Cfg_multigrid* levels;
    bool count_iters_like_single_grid;

    Cfg() : num_levels(0), levels(0) {
        //all done
    }

    ~Cfg() {
        //delete[] levels; // TODO memory leak
        levels = (Cfg_multigrid*)0xDEADBEEF;
    }

    static Cfg v_cycle(unsigned int max_iters, double max_res, unsigned int levels, unsigned int max_iters_pre_post, Cfg_jacobi inner = Cfg_jacobi(5, 0)) {
        assert(levels > 0);
        Cfg cfg = Cfg();
        cfg.count_iters_like_single_grid = false;
        cfg.num_levels = levels;
        cfg.levels = new Cfg_multigrid[levels];
        cfg.levels[0].max_iters = max_iters;
        cfg.levels[0].max_res   = max_res;
        cfg.levels[0].max_iters_jacobi_pre  = max_iters_pre_post;
        cfg.levels[0].max_iters_jacobi_post = max_iters_pre_post;
        cfg.levels[0].max_res_jacobi_pre  = max_res;
        cfg.levels[0].max_res_jacobi_post = max_res;
        for (unsigned int i(1); i < levels-1; ++i) {
            cfg.levels[i].max_iters = 1;
            cfg.levels[i].max_res   = max_res;
            cfg.levels[i].max_iters_jacobi_pre  = max_iters_pre_post;
            cfg.levels[i].max_iters_jacobi_post = max_iters_pre_post;
            cfg.levels[i].max_res_jacobi_pre  = max_res;
            cfg.levels[i].max_res_jacobi_post = max_res;
        }
        if (levels > 1) {
            // if levels == 1 then levels[0] should contain max_iters, max_res
            // not inner.max_iters, inner.max_res
            // *_jacobi_{pre,post} is ignored any ways
            cfg.levels[levels-1].max_iters = inner.max_iters;
            cfg.levels[levels-1].max_res   = inner.max_res;
            cfg.levels[levels-1].max_iters_jacobi_pre  = 0; // ignored
            cfg.levels[levels-1].max_iters_jacobi_post = 0; // ignored
            cfg.levels[levels-1].max_res_jacobi_pre  = 0; // ignored
            cfg.levels[levels-1].max_res_jacobi_post = 0; // ignored
        }
        return cfg;
    }

    static Cfg two_grid(unsigned int max_iters, double max_res, unsigned int max_iters_pre_post, Cfg_jacobi inner = Cfg_jacobi()) {
        return v_cycle(max_iters, max_res, 2, max_iters_pre_post, inner);
    }

    static Cfg jacobi(unsigned int max_iters, double max_res) {
        // inner jacobi doesen't exist --------------,
        // and wont't be used                        v
        Cfg cfg = v_cycle(max_iters, max_res, 1, 0, Cfg_jacobi(0,0));
        cfg.count_iters_like_single_grid = true;
        return cfg;
    }
};

class MultiGrid : public Solver {
public:
    MultiGrid(const Geometry *geom, const Cfg* cfg);
    ~MultiGrid();
    // one MG cycle
    real_t Cycle(Grid *grid, const Grid *rhs) const;

private:
    const Cfg* _cfg;
};

#endif // __SOLVER_HPP
