#include "grid.hpp"
#include "geometry.hpp"
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
    
}

/// Read access to the grid cell at position [it]
const real_t &Grid::Cell(const Iterator &it) const {
    
}


/// Interpolate the value at a arbitrary position
real_t Grid::Interpolate(const multi_real_t &pos) const {
    
}


/// Computes the left-sided difference quatient in x-dim at [it]
real_t Grid::dx_l(const Iterator &it) const {
    
}

/// Computes the right-sided difference quatient in x-dim at [it]
real_t Grid::dx_r(const Iterator &it) const {
    
}

/// Computes the left-sided difference quatient in y-dim at [it]
real_t Grid::dy_l(const Iterator &it) const {
    
}

/// Computes the right-sided difference quatient in x-dim at [it]
real_t Grid::dy_r(const Iterator &it) const {
    
}

/// Computes the central difference quatient of 2nd order in x-dim at [it]
real_t Grid::dxx(const Iterator &it) const {
    
}

/// Computes the central difference quatient of 2nd order in y-dim at [it]
real_t Grid::dyy(const Iterator &it) const {
    
}


/// Computes u*du/dx with the donor cell method
real_t Grid::DC_udu_x(const Iterator &it, const real_t &alpha) const {
    
}

/// Computes v*du/dy with the donor cell method
real_t Grid::DC_vdu_y(const Iterator &it, const real_t &alpha, const Grid *v) const {
    
}

/// Computes u*dv/dx with the donor cell method
real_t Grid::DC_udv_x(const Iterator &it, const real_t &alpha, const Grid *u) const {
    
}

/// Computes v*dv/dy with the donor cell method
real_t Grid::DC_vdv_y(const Iterator &it, const real_t &alpha) const {
    
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

