/*
 *  Copyright (C) 2015   Malte Brunn
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//------------------------------------------------------------------------------
#include "typedef.hpp"
#include "communicator.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
//#include "visu.hpp"
//#include "vtk.hpp"
#include "particle.hpp"

#include <iostream>
#include <sys/stat.h>

#include <fenv.h>

#undef USE_DEBUG_VISU
#define NO_VTK

int main(int argc, char **argv) {
  // we want floating point exeptions
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  // Create parameter and geometry instances with default values
  Communicator comm(&argc, &argv);
  comm.set_boundary_comm(CommBoundary::Sweep);
  Parameter param;
  Geometry geom(&comm);
  if (argc >= 2) {
      param.Load(argv[1], comm.getRank() == 0);
  }
  if (argc >= 3) {
      geom.Load(argv[2]);
  }
  // resize the therads according to the grid
  comm.opt_geom(&geom);
  // Create the fluid solver
  Compute comp(&geom, &param, &comm);

  if (comm.getRank() == 0) {
    // check if folder "VTK" exists
    struct stat info;

    if (stat("VTK", &info) != 0) {
      system("mkdir VTK");
    }
  }

#ifdef USE_DEBUG_VISU
  // Create and initialize the visualization
  Renderer visu(geom.Length(), geom.Mesh());
  visu.Init(1200 / comm.ThreadDim()[0], 1200 / comm.ThreadDim()[1] * geom.Length()[1] / geom.Length()[0],
            comm.getRank() + 1);
#endif // USE_DEBUG_VISU

#ifndef NO_VTK
  // Create a VTK generator
  // use offset as the domain shift
  multi_real_t offset;
  offset[0] = comm.ThreadIdx()[0] * (geom.Mesh()[0] * (double)(geom.Size()[0]));
  offset[1] = comm.ThreadIdx()[1] * (geom.Mesh()[1] * (double)(geom.Size()[1]));
  VTK vtk(geom.Mesh(), geom.Length(), geom.TotalLength(), offset, comm.getRank(),
          comm.getSize(), comm.ThreadDim());
#endif

#ifdef USE_DEBUG_VISU
  const Grid *visugrid;

  visugrid = comp.GetVelocity();
#endif // USE_DEBUG_VISU

  // Run the time steps until the end is reached
  uint32_t iter = 1;
  while (comp.GetTime() < param.Tend()) {
#ifdef USE_DEBUG_VISU
    // Render and check if window is closed
    int vis_mode = visu.Render(visugrid, comp.GetStreakLine());

    // bcast vis_mode to other ranks
    vis_mode = comm.bcast(vis_mode, 0);

    switch (vis_mode) {
    case -1:
      return -1;
    case 0:
      visugrid = comp.GetVelocity();
      break;
    case 1:
      visugrid = comp.GetU();
      break;
    case 2:
      visugrid = comp.GetV();
      break;
    case 3:
      visugrid = comp.GetP();
      break;
    case 4:
      visugrid = comp.GetVorticity();
      break;
    case 5:
      visugrid = comp.GetStream();
      break;
    default:
      break;
    };
#endif // DEBUG_VISU

#ifndef NO_VTK
    // Create VTK Files in the folder VTK
    // Note that when using VTK module as it is you first have to write cell
    // information, then call SwitchToPointData(), and then write point data.
    vtk.Init("VTK/field");
    vtk.AddRank();
    vtk.AddCellField("Cell Velocity", comp.GetU(), comp.GetV());
    vtk.SwitchToPointData();
    vtk.AddPointField("Velocity", comp.GetU(), comp.GetV());
    vtk.AddPointScalar("Pressure", comp.GetP());
    vtk.AddPointScalar("Vorticity", comp.GetVorticity());
    vtk.AddPointScalar("StreamLines", comp.GetStream());
    vtk.Finish();
    comp.GetPathLine()->SaveVTK(comm.getRank(), comm.getSize(), "VTK/path_line", iter);
    comp.GetStreakLine()->SaveVTK(comm.getRank(), comm.getSize(), "VTK/streak_line", iter);
#endif

    // Run a few steps
    for (uint32_t i = 0; i < 9; ++i) {
      comp.TimeStep(false, iter);
      iter++;
    }
    bool printOnlyOnMaster = !comm.getRank();
    comp.TimeStep(printOnlyOnMaster, iter);
    iter++;
  }
  return 0;
}
