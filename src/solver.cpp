#include "compute.hpp"
#include "grid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
#include "solver.hpp"
#include "typedef.hpp"
#include "flaggrid.hpp"
#include "mg_impl.hpp"

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
    typedef Fn_laplace<solver_real_t> F;
    typedef Fn_laplace<char> FLAGS;

    multi_index_t N_fine = this->_geom->Size();
    multi_index_t N_i;

    this->_cfg = cfg;

    u0  = new solver_real_t*[_cfg->num_levels];
    scratch  = F::malloc_typed(F::size_N(N_fine)); //can be used as ping pong buffer by all levels
    b = new solver_real_t*[_cfg->num_levels];
    r_0 = new solver_real_t*[_cfg->num_levels];
    e_0 = new solver_real_t*[_cfg->num_levels];
    res0 = new solver_real_t*[_cfg->num_levels];
    res1 = new solver_real_t*[_cfg->num_levels];
    N = new multi_index_t[_cfg->num_levels];

    // preallocate memory
#define single_buffer
#ifdef single_buffer
    // use one large allocation for all memory needed
    unsigned int total_size = 0;
    N_i = N_fine;
    for (unsigned int i = 0; i < _cfg->num_levels; ++i, N_i = F::coarsen(N_i)) {
        total_size += F::size_N(N_i);
    }
    // 6 buffer hirachies (minus _u0 and _b if they can be used directly)
    buffer = F::malloc_typed(total_size * 6 - (F::mem_device_host_equal()?2*F::size_N(N_fine):0));

    solver_real_t* ptr = buffer;
    N_i = N_fine;
    for (unsigned int i = 0; i < _cfg->num_levels; ++i, N_i = F::coarsen(N_i)) {
        if (i==0 && F::mem_device_host_equal()) {
            // set pointer later
            u0[i] = 0;
            b[i] = 0;
        } else {
            u0[i] = ptr;
            ptr += F::size_N(N_i);
            b[i] = ptr;
            ptr += F::size_N(N_i);
        }

        r_0[i] = ptr;
        ptr += F::size_N(N_i);
        e_0[i] = ptr;
        ptr += F::size_N(N_i);
        res0[i] = ptr;
        ptr += F::size_N(N_i);
        res1[i] = ptr;
        ptr += F::size_N(N_i);

        N[i] = N_i;
    }
#else
    // use different allcations
    N_i = N_fine;
    for (unsigned int i = 0; i < _cfg->num_levels; ++i, N_i = F::coarsen(N_i)) {
        u0[i]  = (i==0 && F::mem_device_host_equal() ? 0 : F::malloc_typed(F::size_N(N_i)));
        b[i] = (i==0 && F::mem_device_host_equal() ? 0 : F::malloc_typed(F::size_N(N_i)));
        r_0[i] = F::malloc_typed(F::size_N(N_i));
        e_0[i] = F::malloc_typed(F::size_N(N_i));
        res0[i] = F::malloc_typed(F::size_N(N_i));
        res1[i] = F::malloc_typed(F::size_N(N_i));
        N[i] = N_i;
    }
#endif

#ifdef MIXED_PRECISION
    unsigned int size_N = Fn_laplace<real_t>::size_N(this->_geom->Size());

    res_f32 = new solver_real_t[size_N];
    err_f32 = new solver_real_t[size_N];
    std::fill_n(res_f32, size_N, 0.0);
#endif

    // prepare flags for the hierarchy
    unsigned int num_levels = cfg->num_levels;
    this->flags = new char*[num_levels];
    N_i = N_fine;
    multi_index_t N_im1 = N_fine;
    for (unsigned int i = 0; i < num_levels; ++i, N_i = F::coarsen(N_i), N_im1 = F::coarsen(N_im1)) {
        if (i==1) N_im1 = N_fine;
        if (i != 0) {
            this->flags[i] = FLAGS::malloc_typed(F::size_N(N_i));
            // run restriction
            restrict_flags_2D(N_im1, flags[i-1], flags[i]);
        } else {
            this->flags[0] = const_cast<char*>(this->_geom->Flags().Data()); // we won't touch this data
        }
    }
}

MultiGrid::~MultiGrid() {
    typedef Fn_laplace<solver_real_t> F;

#ifdef single_buffer
    F::free_typed(buffer);
#else
    N_i = N_fine;
    for (unsigned int i = 0; i < _cfg->num_levels; ++i, N_i = N_i/2) {
        if (i!=0 || !F::mem_device_host_equal()) {
            F::free_typed(u0[i]);
            F::free_typed(b[i]);
        }
        F::free_typed(r_0[i]);
        F::free_typed(e_0[i]);
        F::free_typed(res0[i]);
        F::free_typed(res1[i]);
    }
#endif

#ifdef MIXED_PRECISION
    delete[] res_f32;
    delete[] err_f32;
#endif


    F::free_typed(scratch);
    delete[] N;

    delete[] u0 ;
    delete[] b;
    delete[] r_0;
    delete[] e_0;
    delete[] res0;
    delete[] res1;

    // prepare flags for the hierarchy
    unsigned int num_levels = this->_cfg->num_levels;
    for (unsigned int i = 1; i < num_levels; ++i) {
        // flags[0] points at the initial flag grid and must not be deleted
        delete[] this->flags[i];
    }
    delete this->flags;
}

real_t MultiGrid::Cycle(Grid *p, const Grid *rhs) const {
    InteriorIterator it(this->_geom);

#ifdef MIXED_PRECISION
    // compute f32 residual
    for(it.First(); it.Valid(); it.Next()) {
        if (p->getGeometry()->Flags().Cell(it) == Flags::Fluid) {
            res_f32[it.Value()] = (solver_real_t)localRes(it, p, rhs);
        }
    }

    unsigned int size_N = Fn_laplace<real_t>::size_N(this->_geom->Size());
    std::fill_n(err_f32, size_N, 0.0);
    // solve for f32 error
    solve_mg_flat<Fn_laplace<solver_real_t>>(err_f32, res_f32);

    for(it.First(); it.Valid(); it.Next()) {
        if (p->getGeometry()->Flags().Cell(it) == Flags::Fluid) {
            p->Cell(it) += (real_t)err_f32[it.Value()];
        }
    }
#else
    solve_mg_flat<Fn_laplace<real_t>>(p->Data(), rhs->Data());
#endif

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
