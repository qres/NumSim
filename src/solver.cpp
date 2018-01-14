#include "compute.hpp"
#include "grid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
#include "solver.hpp"
#include "typedef.hpp"
#include "flaggrid.hpp"

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <algorithm>
#include <math.h>

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

template<bool red, bool black>
void solve_SOR(const Geometry* geom, real_t omega, Grid *p, const Grid *rhs) {
    InteriorIterator it(geom);
    multi_real_t h = geom->Mesh();
    if (omega <= 0.0) {
        if (h[0] == h[1]) {
            // choose optimal omega for poisson with hx = hy
            omega =  2.0 / (1.0 + std::sin(3.14159265359 * h[0]));
        } else {
            std::cerr << "ERR: optimal omega parameter can only be computed for h.x = h.y" << std::endl;
        }
    }

    for(it.First(); it.Valid(); it.Next()) {
        if (p->getGeometry()->Flags().Cell(it) == Flags::Fluid) {
            multi_index_t ix_ij = it.Pos();
            bool even = (ix_ij[0] + ix_ij[1]) % 2 == 0;

            if ((even && red) || (!even && black)) {
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
        }
    }
}

real_t SOR::Cycle(Grid *p, const Grid *rhs) const {
    solve_SOR<true,true>(this->_geom, this->_omega, p, rhs);

    InteriorIterator it(this->_geom);
    // compute residual
    real_t res = 0;
    index_t count = 0;
    for(it.First(); it.Valid(); it.Next()) {
        if (p->getGeometry()->Flags().Cell(it) == Flags::Fluid) {
            real_t loc = localRes(it, p, rhs);
            res += loc*loc;
            count++;
        }
    }
    res /= count;
    res = sqrt(res);

    return res;
}


RedOrBlackSOR::RedOrBlackSOR(const Geometry *geom, const real_t &omega) : SOR(geom, omega) {

}

RedOrBlackSOR::~RedOrBlackSOR() {

}

real_t RedOrBlackSOR::RedCycle(Grid *p, const Grid *rhs) const {
    solve_SOR<false,true>(this->_geom, this->_omega, p, rhs);

    InteriorIterator it(this->_geom);
    // compute residual
    real_t res = 0;
    index_t count = 0;
    for(it.First(); it.Valid(); it.Next()) {
        if (p->getGeometry()->Flags().Cell(it) == Flags::Fluid) {
            real_t loc = localRes(it, p, rhs);
            res += loc*loc;
            count++;
        }
    }
    res /= count;
    res = sqrt(res);

    return res;
}

real_t RedOrBlackSOR::BlackCycle(Grid *p, const Grid *rhs) const {
    solve_SOR<true,false>(this->_geom, this->_omega, p, rhs);

    InteriorIterator it(this->_geom);
    // compute residual
    real_t res = 0;
    index_t count = 0;
    for(it.First(); it.Valid(); it.Next()) {
        if (p->getGeometry()->Flags().Cell(it) == Flags::Fluid) {
            real_t loc = localRes(it, p, rhs);
            res += loc*loc;
            count++;
        }
    }
    res /= count;
    res = sqrt(res);

    return res;

}

#include "mg_impl.cpp"


Cfg_jacobi::Cfg_jacobi() :
    max_iters(100000000),
    max_res(1e-12) {
    //all done
}

Cfg_jacobi::Cfg_jacobi(unsigned int max_iters, double max_res) :
    max_iters(max_iters),
    max_res(max_res) {
    //all done
}


Cfg_multigrid::Cfg_multigrid() :
    max_iters(100000000),
    max_res(1e-12),
    max_iters_jacobi_pre(10),
    max_res_jacobi_pre(1e-12),
    max_iters_jacobi_post(10),
    max_res_jacobi_post(1e-12) {
    //all done
}

Cfg::Cfg() : num_levels(0), levels(0) {
    //all done
}

Cfg::~Cfg() {
    //delete[] levels; // TODO memory leak
    levels = (Cfg_multigrid*)0xDEADBEEF;
}

Cfg Cfg::v_cycle(unsigned int max_iters, double max_res, unsigned int levels, unsigned int max_iters_pre_post, Cfg_jacobi inner) {
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

Cfg Cfg::two_grid(unsigned int max_iters, double max_res, unsigned int max_iters_pre_post, Cfg_jacobi inner) {
    return v_cycle(max_iters, max_res, 2, max_iters_pre_post, inner);
}

Cfg Cfg::jacobi(unsigned int max_iters, double max_res) {
    // inner jacobi doesen't exist --------------,
    // and wont't be used                        v
    Cfg cfg = v_cycle(max_iters, max_res, 1, 0, Cfg_jacobi(0,0));
    cfg.count_iters_like_single_grid = true;
    return cfg;
}

MultiGrid::MultiGrid(const Geometry* geom, const Cfg* cfg) : Solver(geom) {
    this->_cfg = cfg;
}

MultiGrid::~MultiGrid() {

}

real_t MultiGrid::Cycle(Grid *p, const Grid *rhs) const {
    solve_v_flat<Fn_laplace<real_t>>(*this->_cfg, this->_geom->Size()[0], p->Data(), rhs->Data()); // TODO: different resolution in each dimenson
    InteriorIterator it(this->_geom);
    // compute residual
    real_t res = 0;
    index_t count = 0;
    for(it.First(); it.Valid(); it.Next()) {
        if (p->getGeometry()->Flags().Cell(it) == Flags::Fluid) {
            real_t loc = localRes(it, p, rhs);
            res += loc*loc;
            count++;
        }
    }
    res /= count;
    res = sqrt(res);

    return res;
}
