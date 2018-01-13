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
#include <vector>
//------------------------------------------------------------------------------
#ifndef __PARTICLE_HPP
#define __PARTICLE_HPP
//------------------------------------------------------------------------------
class Particle {
public:
  Particle (const multi_real_t& pos);
  ~Particle ();
  void TimeStep(const real_t& dt, const Grid* u, const Grid* v);
  const multi_real_t& Pos () const;
private:
  multi_real_t _pos;
};
//------------------------------------------------------------------------------
class ParticleLine {
public:
  ParticleLine();
  virtual void TimeStep(const real_t& dt, const Grid* u, const Grid* v) = 0;
  virtual void SaveVTK(const index_t& rank, const index_t& nump, const char* basename, const index_t& idx) const;
  const std::vector<Particle>& GetParticles() const;
protected:
  std::vector<Particle> _part;
};
//------------------------------------------------------------------------------
class PathLine : public ParticleLine {
public:
  PathLine (const multi_real_t& pos);
  virtual ~PathLine ();
  void TimeStep(const real_t& dt, const Grid* u, const Grid* v);
};
//------------------------------------------------------------------------------
class StreakLine : public ParticleLine {
public:
  StreakLine (const multi_real_t& pos);
  virtual ~StreakLine ();
  void TimeStep(const real_t& dt, const Grid* u, const Grid* v);
private:
  multi_real_t _org;
};
//------------------------------------------------------------------------------
#endif  // __PARTICLE_HPP
