#include "compute.hpp"
#include "grid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
#include "solver.hpp"
#include "typedef.hpp"

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <algorithm>


/// Creates a compute instance with given geometry and parameter
Compute::Compute(const Geometry *geom, const Parameter *param) {
    this->_t = 0;
    this->_dtlimit = param->Dt();
    this->_epslimit = param->Eps();

    /*

        Upper left:                     Upper Right

        +----v----+--                   --+----v----+
        |         |                       |         |
        |    p    u                       |    p    u
        |         |                       |         |
        +---------+--   --boundary--    --+---------+
        |         |                       |         |
                boundary              boundary
                  |                       |

        Origins for grids:               Lower Right:

        |         | inner points          |         |
        +----v----+--    --boundary--   --+----v----+
        |         |                       |         |
        |    p    u                       |    p    u
        |         |                       |         |
        +---------+--                   --+---------+
        \___global position (0,0)
    */

    multi_real_t h = geom->Mesh();

    // velocities
    this->_u = new Grid(geom, multi_real_t(h[0],       h[1] / 2.0));
    this->_v = new Grid(geom, multi_real_t(h[0] / 2.0, h[1]      ));
    // pressure
    this->_p = new Grid(geom, multi_real_t(h[0] / 2.0, h[1] / 2.0));

    // prel. vel
    this->_F = new Grid(geom); // TODO same offset as u and v?
    this->_G = new Grid(geom);

    // right-hand side
    this->_rhs = new Grid(geom); // TODO same offset as p?

    this->_tmp = new Grid(geom);

    _solver = new SOR(geom, param->Omega());

    this->_geom = geom;
    this->_param = param;

    // TODO init values to zero?
    this->_u->Initialize(0);
    this->_v->Initialize(0);
    this->_p->Initialize(0);
    this->_rhs->Initialize(0);
    this->_F->Initialize(0);
    this->_G->Initialize(0);
}

/// Deletes all grids
Compute::~Compute() {
    delete this->_u;
    delete this->_v;
    delete this->_p;
    delete this->_F;
    delete this->_G;
    delete this->_rhs;
    delete this->_tmp;
    delete this->_solver;
}

/// Execute one time step of the fluid simulation (with or without debug info)
// @ param printInfo print information about current solver state (residual
// etc.)
void Compute::TimeStep(bool printInfo) {
    // update boundaries -- in no specific order
    this->_geom->Update_U(this->_u);
    this->_geom->Update_V(this->_v);
    this->_geom->Update_P(this->_p);

    // timestep
    multi_real_t h = this->_geom->Mesh();
    real_t u_max = this->_u->AbsMax();
    real_t v_max = this->_v->AbsMax();
    real_t re = this->_param->Re();
    real_t dt = std::min(std::min(
        // stability conditions for convection operator:
        h[0] / u_max,
        h[1] / v_max
    ),  // stability condition for diffusion operator:
        re / 2 * (h[0]*h[0]*h[1]*h[1]) / (h[0]*h[0] + h[1]*h[1])
    ) * this->_param->Tau(); // scaling in range (0, 1)

    // if an explicit timestep is given, use this timestep instead of the computed timestep
    if (this->_param->Dt() > 0.0) {
        dt = this->_param->Dt();
    }

    // compute F and G
    this->MomentumEqu(dt);
    // solve poisson equation for p
    this->RHS(dt);
    real_t res = 12345678.9;
    index_t iter = 0;
    while (res > this->_epslimit && iter < this->_param->IterMax()) {
        res = this->_solver->Cycle(this->_p, this->_rhs);
        iter++;
    }
    // compute new velocities from F and G and the p derivatives
    this->NewVelocities(dt);

    // increment simulated time
    this->_t += dt;

    if (printInfo) {
        std::cout << "[t=" << this->_t << "] Solved poisson eq in " << iter << " iterations and final residual of " << res << std::endl;
    }
    
}

/// Returns the simulated time in total
const real_t &Compute::GetTime() const {
    return this->_t;
}


/// Returns the pointer to U
const Grid *Compute::GetU() const {
    return this->_u;
}
/// Returns the pointer to V
const Grid *Compute::GetV() const {
    return this->_v;
}
/// Returns the pointer to P
const Grid *Compute::GetP() const {
    return this->_p;
}
/// Returns the pointer to RHS
const Grid *Compute::GetRHS() const {
    return this->_rhs;
}

/// Computes and returns the absolute velocity
const Grid *Compute::GetVelocity() {
    Iterator it = Iterator(this->_geom);

    for (it.First(); it.Valid(); it.Next()) {
        real_t u = (this->_u->Cell(it) + this->_u->Cell(it.Left())) / 2.0;
        real_t v = (this->_v->Cell(it) + this->_v->Cell(it.Down())) / 2.0;
		_tmp->Cell(it) = sqrt(u*u + v*v);
	}

    return _tmp;
}
/// Computes and returns the vorticity
const Grid *Compute::GetVorticity() {
    throw std::runtime_error("Unimplemented");
}
/// Computes and returns the stream line values
const Grid *Compute::GetStream() {
    throw std::runtime_error("Unimplemented");
}


/// Compute the temporary velocites F,G
void Compute::MomentumEqu(const real_t &dt) {
    // calculate F
    InteriorIterator it = InteriorIterator(this->_geom);

    multi_real_t h = this->_geom->Mesh();
    real_t Re_inv = 1.0 / this->_param->Re();
    real_t alpha = this->_param->Alpha();
    real_t gx = 0.0;
    real_t gy = 0.0;

    for (it.First(); it.Valid(); it.Next()) {
        real_t laplace_u = _u->dxx(it) + _u->dyy(it);
        real_t duudx = _u->DC_duu_dx(it, alpha);
        real_t duvdy = _u->DC_duv_dy(it, alpha, _v);

        real_t A_ij = Re_inv * laplace_u - duudx - duvdy + gx;
        this->_F->Cell(it) = this->_u->Cell(it) + dt * A_ij;

        real_t laplace_v = _v->dxx(it) + _v->dyy(it);
        real_t dvvdy = _v->DC_dvv_dy(it, alpha);
        real_t duvdx = _v->DC_duv_dx(it, alpha, _u);

        real_t B_ij = Re_inv * laplace_v - duvdx - dvvdy + gy;
        this->_G->Cell(it) = this->_v->Cell(it) + dt * B_ij;
    }
}

/// Compute the new velocites u,v
void Compute::NewVelocities(const real_t &dt) {
    InteriorIterator it = InteriorIterator(this->_geom);
    
    multi_real_t h = this->_geom->Mesh();
    
    for (it.First(); it.Valid(); it.Next()) {
        // u_m+1 = F_n - dt dpdx_n+1
        real_t dpdx = this->_p->dx_l(it);
        real_t F    = this->_F->Cell(it);
        this->_u->Cell(it) = F - dt * dpdx;
        
        real_t dpdy = this->_p->dy_l(it);
        real_t G    = this->_G->Cell(it);
        this->_v->Cell(it) = G - dt * dpdy;
    }
}

/// Compute the RHS of the poisson equation
void Compute::RHS(const real_t &dt) {
    InteriorIterator it = InteriorIterator(this->_geom);

    multi_real_t h = this->_geom->Mesh();

    for (it.First(); it.Valid(); it.Next()) {
        real_t F_i   = this->_F->Cell(it);
        real_t F_im1 = this->_F->Cell(it.Left());
        real_t G_j   = this->_G->Cell(it);
        real_t G_jm1 = this->_G->Cell(it.Down());
        
        this->_rhs->Cell(it) = 1.0/dt * (
            (F_i - F_im1) / h[0] + 
            (G_j - G_jm1) / h[1]
        );
    }
}
