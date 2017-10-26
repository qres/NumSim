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
    real_t x = pos[0] - this->_offset[0];
    real_t y = pos[1] - this->_offset[1];
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
    return (f_i - f_ip1) / h;
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
    return (f_j - f_jp1) / h;
}

/// Computes the central difference quatient of 2nd order in x-dim at [it]
real_t Grid::dxx(const Iterator &it) const {
    real_t f_i   = this->Cell(it);
    real_t f_im1 = this->Cell(it.Left());
    real_t f_ip1 = this->Cell(it.Right());
    real_t h = this->_geom->Mesh()[0];
    return (f_im1 - 2 * f_i + f_im1) / (h*h);
}

/// Computes the central difference quatient of 2nd order in y-dim at [it]
real_t Grid::dyy(const Iterator &it) const {
    real_t f_j   = this->Cell(it);
    real_t f_jm1 = this->Cell(it.Down());
    real_t f_jp1 = this->Cell(it.Top());
    real_t h = this->_geom->Mesh()[1];
    return (f_jm1 - 2 * f_j + f_jm1) / (h*h);
}


/// Computes u*du/dx with the donor cell method
real_t Grid::DC_udu_x(const Iterator &it, const real_t &alpha) const {
    // TODO DC
    real_t u    = this->Cell(it);
    real_t dudx = this->dx_l(it);
    return u*dudx; 
}

/// Computes v*du/dy with the donor cell method
real_t Grid::DC_vdu_y(const Iterator &it, const real_t &alpha, const Grid *v_grid) const {
    // TODO DC
    real_t v    = v_grid->Cell(it);
    real_t dudy = this->dy_l(it);
    return v*dudy;
}

/// Computes u*dv/dx with the donor cell method
real_t Grid::DC_udv_x(const Iterator &it, const real_t &alpha, const Grid *u_grid) const {
    // TODO DC
    real_t u    = u_grid->Cell(it);
    real_t dvdx = this->dx_l(it);
    return u*dvdx;
}

/// Computes v*dv/dy with the donor cell method
real_t Grid::DC_vdv_y(const Iterator &it, const real_t &alpha) const {
    // TODO DC
    real_t v    = this->Cell(it);
    real_t dvdx = this->dy_l(it);
    return v*dvdx; 
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
    return *std::max_element(this->_data, this->_data + size, [=](real_t i, real_t j){return max(abs(i), abs(j));});
}



real_t *Grid::Data() {
    return this->_data;
}

