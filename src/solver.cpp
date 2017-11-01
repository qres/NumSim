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

Solver::Solver(const Geometry *geom) {
    this->_geom = geom;
}

Solver::~Solver() {

}

real_t Solver::localRes(const Iterator &it, const Grid *p, const Grid *rhs) const {
    /*
          1
        1 -4 1  p  = rhs
          1
    */
    real_t ddp_dxx = p->dxx(it);
    real_t ddp_dyy = p->dyy(it);
    return rhs->Cell(it) - (ddp_dxx + ddp_dyy);
}

SOR::SOR(const Geometry *geom, const real_t &omega) : Solver(geom) {
    this->_omega = omega;
}

SOR::~SOR() {

}


real_t SOR::Cycle(Grid *p, const Grid *rhs) const {
    InteriorIterator it(this->_geom);
    real_t omega = this->_omega;
    multi_real_t h = this->_geom->Mesh();

    for(it.First(); it.Valid(); it.Next()) {
        real_t p_ij = p->Cell(it);
        real_t p_im1j = p->Cell(it.Left());
        real_t p_ip1j = p->Cell(it.Right());
        real_t p_ijm1 = p->Cell(it.Down());
        real_t p_ijp1 = p->Cell(it.Top());
        real_t rhs_ij = rhs->Cell(it);

        real_t dxsq = h[0] * h[0];
        real_t dysq = h[1] * h[1];

        p->Cell(it) = 
            (1 - omega) * p_ij + 
            omega * 0.5 * (dxsq*dysq) / (dxsq+dysq) * (
                (p_im1j + p_ip1j) / dxsq +
                (p_ijm1 + p_ijp1) / dysq +
                -rhs_ij
        );
    }

    // compute residual
    real_t res = 0;
    for(it.First(); it.Valid(); it.Next()) {
        res = std::max(res, std::abs(localRes(it, p, rhs)));
    }

    return res;
}
