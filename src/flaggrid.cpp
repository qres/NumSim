
#include "flaggrid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include <iostream>

FlagGrid::FlagGrid(const Geometry *geom) {
    this->_geom = geom;

    if (geom) {
        multi_index_t N = geom->Size();
        // add boundary
        size_t size = (N[0] + 2) * (N[1] + 2);

        this->_data = new char[size];
        this->Initialize(' ');
    } else {
        this->_data = 0;
    }
}

FlagGrid::~FlagGrid() {
    delete[] this->_data;
}

/// Initializes the grid with a value
void FlagGrid::Initialize(const char &value) {
    multi_index_t N = this->_geom->Size();
    // add boundary
    size_t size = (N[0] + 2) * (N[1] + 2); // TODO duplicate
    std::fill_n(this->_data, size, value);
}

/// Write access to the grid cell at position [it]
char &FlagGrid::Cell(const Iterator &it) {
    return this->_data[it.Value()];
}

/// Read access to the grid cell at position [it]
const char &FlagGrid::Cell(const Iterator &it) const {
    return this->_data[it.Value()];
}

/// fluid     = .....
/// hor t     = 1...1
/// vert r    = 1  1.
/// vert l    = 1 1..
/// hor b     = 11...
/// corner tl = 1.1.1
/// corner tr = 1..11
/// corner bl = 111..
/// corner br = 11.1.
/// interior  = 1....
char FlagGrid::BoundaryOrientation(const Iterator &it) const {
    auto self   = this->Cell(it);
    if (self == Flags::Fluid) return 0;
    auto left   = this->Cell(it.Left());
    auto right  = this->Cell(it.Right());
    auto bottom = this->Cell(it.Down());
    auto top    = this->Cell(it.Top());

    // we alredy know that self is not a fluid
    bool vert_r  = right  == Flags::Fluid;
    bool vert_l  = left   == Flags::Fluid;
    bool horiz_t = top    == Flags::Fluid;
    bool horiz_b = bottom == Flags::Fluid;

    return 1<<4 | horiz_b << 3 | vert_l << 2 | vert_r << 1 | horiz_t << 0;
}

void FlagGrid::set_driven_cavity() {
    BoundaryIterator bit (this->_geom);
    bit.SetBoundary(Boundary::Bottom);
    for (bit.First(); bit.Valid(); bit.Next()) {
        this->Cell(bit) = Flags::Noslip;
    }
    bit.SetBoundary(Boundary::Left);
    for (bit.First(); bit.Valid(); bit.Next()) {
        this->Cell(bit) = Flags::Noslip;
    }
    bit.SetBoundary(Boundary::Right);
    for (bit.First(); bit.Valid(); bit.Next()) {
        this->Cell(bit) = Flags::Noslip;
    }
    bit.SetBoundary(Boundary::Top);
    for (bit.First(); bit.Valid(); bit.Next()) {
        this->Cell(bit) = Flags::Inflow;
    }
}

char* FlagGrid::Data() {
    return this->_data;
}
