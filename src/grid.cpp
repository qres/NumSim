#include "grid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "typedef.hpp"

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <algorithm>

using std::max;
using std::min;
using std::abs;

/// \brief Allocates the grid defined by \a geom with no offset
/// the first boundary cell is therefore aligned to (0,0)
Grid::Grid(const Geometry *geom) : Grid(geom, multi_real_t(0.0, 0.0)) {

}

/// \brief Allocates the grid defined by \a geom
/// \param offset position offset of the first boundary grid point (this is needed for staggered grids)
Grid::Grid(const Geometry *geom, const multi_real_t &offset) {
    this->_geom = geom;
    this->_offset = offset;

    multi_index_t N = geom->Size();
    // add boundary
    size_t size = (N[0] + 2) * (N[1] + 2);

    this->_data = new real_t[size];
}

/// Deletes the grid
Grid::~Grid() {
    delete[] this->_data;
}

/// Initializes the grid with a value
void Grid::Initialize(const real_t &value) {
    multi_index_t N = this->_geom->Size();
    // add boundary
    size_t size = (N[0] + 2) * (N[1] + 2); // TODO duplicate
    std::fill_n(this->_data, size, value);
}


/// Write access to the grid cell at position [it]
real_t &Grid::Cell(const Iterator &it) {
    return this->_data[it.Value()];
}

/// Read access to the grid cell at position [it]
const real_t &Grid::Cell(const Iterator &it) const {
    return this->_data[it.Value()];
}


/// Interpolate the value at a arbitrary position
real_t Grid::Interpolate(const multi_real_t &pos) const {
    real_t x = std::max(pos[0] - this->_offset[0], 0.0);
    real_t y = std::max(pos[1] - this->_offset[1], 0.0);
    multi_real_t h = this->_geom->Mesh();

    multi_index_t size = this->_geom->Size();
    index_t stride_y = size[0] + 2;

    // pos is in range [i,i+1)x[j,j+1)
    index_t i = (index_t)(x / h[0]);
    index_t j = (index_t)(y / h[1]);

    // distance of x to x1,x2,y1,y2
    real_t dx1 = x - i * h[0];
    real_t dy1 = y - j * h[1];
    real_t dx2 = (i+1) * h[0] - x;
    real_t dy2 = (j+1) * h[1] - y;

    real_t u      = this->_data[i      +  j       * stride_y];
    real_t u_ip1  = this->_data[i + 1  +  j       * stride_y];
    real_t u_jp1  = this->_data[i      +  (j + 1) * stride_y];
    real_t u_ijp1 = this->_data[i + 1  +  (j + 1) * stride_y];

    // bilinear interpolation
    return 1.0/(h[0]*h[1]) * (
        u      * dx2 * dy2 +
        u_ip1  * dx1 * dy2 +
        u_jp1  * dx2 * dy1 +
        u_ijp1 * dx1 * dy1
    );
}


/// Computes the left-sided difference quatient in x-dim at [it]
real_t Grid::dx_l(const Iterator &it) const {
    real_t f_i   = this->Cell(it);
    real_t f_im1 = this->Cell(it.Left());
    real_t h = this->_geom->Mesh()[0];
    return (f_i - f_im1) / h;
}

/// Computes the right-sided difference quatient in x-dim at [it]
real_t Grid::dx_r(const Iterator &it) const {
    real_t f_i   = this->Cell(it);
    real_t f_ip1 = this->Cell(it.Right());
    real_t h = this->_geom->Mesh()[0];
    return (f_ip1 - f_i) / h;
}

/// Computes the left-sided difference quatient in y-dim at [it]
real_t Grid::dy_l(const Iterator &it) const {
    real_t f_j   = this->Cell(it);
    real_t f_jm1 = this->Cell(it.Down());
    real_t h = this->_geom->Mesh()[1];
    return (f_j - f_jm1) / h;
}

/// Computes the right-sided difference quatient in x-dim at [it]
real_t Grid::dy_r(const Iterator &it) const {
    real_t f_j   = this->Cell(it);
    real_t f_jp1 = this->Cell(it.Top());
    real_t h = this->_geom->Mesh()[1];
    return (f_jp1 - f_j) / h;
}

/// Computes the central difference quatient of 2nd order in x-dim at [it]
real_t Grid::dxx(const Iterator &it) const {
    real_t f_i   = this->Cell(it);
    real_t f_im1 = this->Cell(it.Left());
    real_t f_ip1 = this->Cell(it.Right());
    real_t h = this->_geom->Mesh()[0];
    return (f_im1 - 2 * f_i + f_ip1) / (h*h);
}

/// Computes the central difference quatient of 2nd order in y-dim at [it]
real_t Grid::dyy(const Iterator &it) const {
    real_t f_j   = this->Cell(it);
    real_t f_jm1 = this->Cell(it.Down());
    real_t f_jp1 = this->Cell(it.Top());
    real_t h = this->_geom->Mesh()[1];
    return (f_jm1 - 2 * f_j + f_jp1) / (h*h);
}


/// Computes d(u^2)/dx with the donor cell method
real_t Grid::DC_duu_dx(const Iterator &it, const real_t &alpha) const {
    real_t h = this->_geom->Mesh()[0];
    real_t u_i   = this->Cell(it);
    real_t u_im1 = this->Cell(it.Left());
    real_t u_ip1 = this->Cell(it.Right());
    real_t u_ip12 = (u_i + u_ip1)/2.0;
    real_t u_im12 = (u_i + u_im1)/2.0;

    real_t central = 1.0/h * (u_ip12*u_ip12 - u_im12*u_im12);

    real_t donor_cell = 1.0/h * (
        ((u_ip12 > 0) ? (u_ip12*u_i)   : (u_ip12*u_ip1)) -
        ((u_im12 > 0) ? (u_im12*u_im1) : (u_im12*u_i))
    );

    return (1 - alpha) * central + alpha * donor_cell;
}

/// Computes d(uv)/dy with the donor cell method
real_t Grid::DC_duv_dy(const Iterator &it, const real_t &alpha, const Grid *v_grid) const {
    real_t h = this->_geom->Mesh()[1];
    real_t u      = this->Cell(it);
    real_t u_jm1  = this->Cell(it.Down());
    real_t u_jp1  = this->Cell(it.Top());
    real_t u_jp12 = (u + u_jp1)/2.0;
    real_t u_jm12 = (u + u_jm1)/2.0;
    real_t v         = v_grid->Cell(it);
    real_t v_ip1     = v_grid->Cell(it.Right());
    real_t v_jm1     = v_grid->Cell(it.Down());
    real_t v_ip1jm1  = v_grid->Cell(it.Right().Down());
    real_t v_ip12    = (v + v_ip1)/2.0;
    real_t v_ip12jm1 = (v_jm1 + v_ip1jm1)/2.0;

    real_t central = 1.0 / h * (v_ip12*u_jp12 - v_ip12jm1*u_jm12);

    real_t a;
    real_t b;
    if (v_ip12 > 0) {
        a = v_ip12*u;
    } else {
        a = v_ip12*u_jp1;
    }
    if (v_ip12jm1 > 0) {
        b = v_ip12jm1*u_jm1;
    } else {
        b = v_ip12jm1*u;
    }
    real_t donor_cell = 1.0 / h * (a - b);

    return (1 - alpha) * central + alpha * donor_cell;
}

/// Computes d(uv)/dx with the donor cell method
real_t Grid::DC_duv_dx(const Iterator &it, const real_t &alpha, const Grid *u_grid) const {
    real_t h = this->_geom->Mesh()[0];
    real_t u         = u_grid->Cell(it);
    real_t u_jp1     = u_grid->Cell(it.Top());
    real_t u_im1     = u_grid->Cell(it.Left());
    real_t u_im1jp1  = u_grid->Cell(it.Left().Top());
    real_t u_im1jp12 = (u_im1 + u_im1jp1)/2.0;
    real_t u_jp12    = (u + u_jp1)/2.0;
    real_t v      = this->Cell(it);
    real_t v_ip1  = this->Cell(it.Right());
    real_t v_im1  = this->Cell(it.Left());
    real_t v_im12 = (v + v_im1)/2.0;
    real_t v_ip12 = (v + v_ip1)/2.0;

    real_t central = 1.0 / h * (v_ip12*u_jp12 - v_im12*u_im1jp12);

    real_t a;
    real_t b;
    if (u_jp12 > 0) {
        a = u_jp12 * v;
    } else {
        a = u_jp12 * v_ip1;
    }
    if (u_im1jp12 > 0) {
        b = u_im1jp12 * v_im1;
    } else {
        b = u_im1jp12 * v;
    }
    real_t donor_cell = 1.0 / h * (a - b);

    return (1 - alpha) * central + alpha * donor_cell;
}

/// Computes d(v^2)/dy with the donor cell method
real_t Grid::DC_dvv_dy(const Iterator &it, const real_t &alpha) const {
    real_t h = this->_geom->Mesh()[1];
    real_t v_j   = this->Cell(it);
    real_t v_jm1 = this->Cell(it.Down());
    real_t v_jp1 = this->Cell(it.Top());
    real_t v_jp12 = (v_j + v_jp1)/2.0;
    real_t v_jm12 = (v_j + v_jm1)/2.0;

    real_t central = 1.0/h * (v_jp12*v_jp12 - v_jm12*v_jm12);

    real_t donor_cell = 1.0/h * (
        ((v_jp12 > 0) ? (v_jp12*v_j)   : (v_jp12*v_jp1)) -
        ((v_jm12 > 0) ? (v_jm12*v_jm1) : (v_jm12*v_j))
    );

    return (1 - alpha) * central + alpha * donor_cell;
}


/// Returns the maximal value of the grid
real_t Grid::Max() const {
    multi_index_t N = this->_geom->Size();
    // add boundary
    size_t size = (N[0] + 2) * (N[1] + 2); // TODO duplicate
    return *std::max_element(this->_data, this->_data + size);
}

/// Returns the minimal value of the grid
real_t Grid::Min() const {
    multi_index_t N = this->_geom->Size();
    // add boundary
    size_t size = (N[0] + 2) * (N[1] + 2); // TODO duplicate
    return *std::min_element(this->_data, this->_data + size);
}

/// Returns the absolute maximal value
real_t Grid::AbsMax() const {
    multi_index_t N = this->_geom->Size();
    // add boundary
    size_t size = (N[0] + 2) * (N[1] + 2); // TODO duplicate
    return abs(*std::max_element(this->_data, this->_data + size, [=](real_t i, real_t j){return abs(i) < abs(j);}));
}



real_t *Grid::Data() {
    return this->_data;
}

const Geometry *Grid::getGeometry() const {
    return this->_geom;
}
