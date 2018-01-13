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

#include "particle.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include <stdio.h>

Particle::Particle (const multi_real_t& pos) {
    this->_pos = pos;
}

Particle::~Particle () {

}

void Particle::TimeStep(const real_t& dt, const Grid* u, const Grid* v) {
    real_t uu = u->Interpolate(this->Pos());
    real_t vv = v->Interpolate(this->Pos());
    this->_pos[0] += dt * uu;
    this->_pos[1] += dt * vv;
}

const multi_real_t& Particle::Pos () const {
    return this->_pos;
}

ParticleLine::ParticleLine(): _part() {

}

const std::vector<Particle>& ParticleLine::GetParticles() const {
    return _part;
}


PathLine::PathLine(const multi_real_t& pos) {
    this->_part.push_back(Particle(pos));

}

PathLine::~PathLine() {

}

void PathLine::TimeStep(const real_t& dt, const Grid* u, const Grid* v) {
    Particle last = this->_part.back();
    last.TimeStep(dt, u, v);
    this->_part.push_back(last);
}

StreakLine::StreakLine(const multi_real_t& pos) {
    this->_part.push_back(Particle(pos));
    this->_org = pos;
}

StreakLine::~StreakLine() {

}

void StreakLine::TimeStep(const real_t& dt, const Grid* u, const Grid* v) {
    for (Particle& particle : this->_part) {
        particle.TimeStep(dt, u, v);
    }
    this->_part.push_back(Particle(this->_org));
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void ParticleLine::SaveVTK (const index_t& rank, const index_t& nump, const char* basename, const index_t& idx) const {
  char fname[1000];
  const char *filebase = basename;
  for (int i = 0; basename[i] != 0; ++i) {
    if (basename[i] == '/') filebase = &basename[i+1];
    if (basename[i] == '\\') filebase = &basename[i+1];
  }
  if (rank == 0) {
    sprintf(fname, "%s_%05i.pvtp", basename, idx);
    FILE *handle = fopen(fname, "w");
    fprintf(handle, "<?xml version=\"1.0\"?>\n"
                    "<VTKFile type=\"PPolyData\" version=\"0.1\" "
                    "byte_order=\"LittleEndian\">\n\t"
                    "<PPolyData GhostLevel=\"0\">\n\t\t<PPoints>\n\t\t\t"
                    "<PDataArray type=\"Float32\" Name=\"Position\" "
                    "NumberOfComponents=\"3\" format=\"ascii\"/>\n"
                    "\t\t</PPoints>\n\t\t"
                    "<PPointData Scalars=\"mpirank\">\n\t\t\t"
                    "<PDataArray type=\"Int32\" Name=\"mpirank\" "
                    "format=\"ascii\"/>\n\t\t\t"
                    "</PPointData>\n");
    for (unsigned int p = 0; p < nump; ++p)
      fprintf(handle, "\t\t<Piece Source=\"%s_%05i_%04i.vtp\"/>\n", filebase, idx, rank);
    fprintf(handle, "\t</PPolyData>\n</VTKFile>\n");
    fclose(handle);
  }
  sprintf(fname, "%s_%05i_%04i.vtp", basename, idx, rank);
  FILE *handle = fopen(fname, "w");
  fprintf(handle, "<?xml version=\"1.0\"?>\n"
                  "<VTKFile type=\"PolyData\" version=\"0.1\" "
                  "byte_order=\"LittleEndian\">\n\t");
  fprintf(handle, "<PolyData>\n\t\t<Piece NumberOfPoints=\"%i\" NumberOfVerts=\"0\" ", _part.size());
  fprintf(handle, "NumberOfLines=\"0\" NumberOfStrips=\"0\" "
                  "NumberOfPolys=\"0\">\n\t\t\t<Points>\n\t\t\t\t"
                  "<DataArray type=\"Float32\" Name=\"Position\" "
                  "NumberOfComponents=\"3\" format=\"ascii\">\n");
  for (unsigned int p = 0; p < _part.size(); ++p) {
    fprintf(handle, "\t\t\t\t\t%le %le %le\n", _part[p].Pos()[0], _part[p].Pos()[1], 0.0);
  }
  fprintf(handle, "\t\t\t\t</DataArray>\n\t\t\t</Points>\n\t\t\t"
                  "<PointData Scalars=\"mpirank\">\n\t\t\t\t"
                  "<DataArray type=\"Int32\" Name=\"mpirank\" "
                  "format=\"ascii\">\n\t\t\t\t\t");
  for (unsigned int p = 0; p < _part.size(); ++p) {
    fprintf(handle, "%i ", rank);
  }
  fprintf(handle, "\n\t\t\t\t</DataArray>\n\t\t\t</PointData>\n"
                  "\t\t</Piece>\n\t</PolyData>\n</VTKFile>\n");
  fclose(handle);
}
//------------------------------------------------------------------------------
