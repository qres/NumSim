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
//------------------------------------------------------------------------------
#ifndef __COMPUTE_HPP
#define __COMPUTE_HPP


// use Red-Black-SOR
//#define USE_RB_SOR
// use Multi Grid
#define USE_MG


//------------------------------------------------------------------------------
class Compute {
public:
  /// Creates a compute instance with given geometry and parameter
  Compute(const Geometry *geom, const Parameter *param,
          const Communicator *comm = 0);
  /// Deletes all grids
  ~Compute();

  /// Execute one time step of the fluid simulation (with or without debug info)
  // @ param printInfo print information about current solver state (residual
  // etc.)
  void TimeStep(bool printInfo, uint32_t iter);

  /// Returns the simulated time in total
  const real_t &GetTime() const;

  /// Returns the pointer to U
  const Grid *GetU() const;
  /// Returns the pointer to V
  const Grid *GetV() const;
  /// Returns the pointer to P
  const Grid *GetP() const;
  /// Returns the pointer to RHS
  const Grid *GetRHS() const;

  /// Computes and returns the absolute velocity
  const Grid *GetVelocity();
  /// Computes and returns the vorticity
  const Grid *GetVorticity();
  /// Computes and returns the stream line values
  const Grid *GetStream();

  const StreakLine *GetStreakLine() const;
  const PathLine *GetPathLine() const;

private:
  // current timestep
  real_t _t;

  // donor-cell diffusion condition (p. 27)
  real_t _dtlimit;

  // limit for residual
  real_t _epslimit;

  // velocities
  Grid *_u;
  Grid *_v;

  // pressure
  Grid *_p;

  // prel. vel
  Grid *_F;
  Grid *_G;

  // right-hand side
  Grid *_rhs;

  // container for interpolating whichever values
  Grid *_velocity;
  Grid *_vorticity;
  Grid *_stream;

  StreakLine *_streak_line;
  PathLine *_path_line;

  #ifdef USE_MG
    MultiGrid *_solver;
    Cfg *_cfg;
  #elif defined USE_RB_SOR
    RedOrBlackSOR *_solver;
  #else
    Solver *_solver;
  #endif

  const Geometry *_geom;
  const Parameter *_param;
  const Communicator *_comm;

  /// Compute the new velocites u,v from F,G
  void NewVelocities(const real_t &dt);
  /// Compute the temporary velocites F,G
  void MomentumEqu(const real_t &dt);
  /// Compute the RHS of the poisson equation from F,G
  void RHS(const real_t &dt);
};
//------------------------------------------------------------------------------
#endif // __COMPUTE_HPP
