#include "geometry.hpp"
#include "grid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
#include "typedef.hpp"

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <algorithm>

/// \brief sets geometry for Driven Cavity
Geometry::Geometry() : 
    _size(128, 128),
    _length(1.0, 1.0),
    _h(1.0/(128 + 1), 1.0/(128 + 1)),
    // velocity at upper boundary
    _velocity(1.0, 0.0),
    _pressure(0.1)  {
    
}

/// \brief load geometry settings form file
void Geometry::Load(const char *file) {
    // TODO
}


const multi_index_t &Geometry::Size() const {
    return this->_size;
}

const multi_real_t &Geometry::Length() const {
    return this->_length;
}

const multi_real_t &Geometry::Mesh() const {
    return this->_h;
}


void Geometry::Update_U(Grid *u) const {
    BoundaryIterator it (this);

    it.SetBoundary(Boundary::Bottom);
    for(it.First(); it.Valid(); it.Next()) {
        // no slip => define boundary s.t. the interpolated value is zero
        u->Cell(it) = - u->Cell(it.Top());
    }

    it.SetBoundary(Boundary::Left);
    for(it.First(); it.Valid(); it.Next()) {
        // no slip
        u->Cell(it) = 0.0;
    }

    it.SetBoundary(Boundary::Right);
    for(it.First(); it.Valid(); it.Next()) {
        // no slip, staggered grid => set left (= boundary) point to zero
        u->Cell(it) = 0.0; // TODO what happens with this point?
        u->Cell(it.Left()) = 0.0;
    }

    it.SetBoundary(Boundary::Top);
    for(it.First(); it.Valid(); it.Next()) {
        // velocity given
        u->Cell(it) = 2* this->_velocity[0] - u->Cell(it.Down());
    }
}

void Geometry::Update_V(Grid *v) const {
    BoundaryIterator it (this);

    it.SetBoundary(Boundary::Bottom);
    for(it.First(); it.Valid(); it.Next()) {
        // no slip
        v->Cell(it) = 0.0;
    }

    it.SetBoundary(Boundary::Left);
    for(it.First(); it.Valid(); it.Next()) {
        // no slip => define boundary s.t. the interpolated value is zero
        v->Cell(it) = - v->Cell(it.Right());
    }

    it.SetBoundary(Boundary::Right);
    for(it.First(); it.Valid(); it.Next()) {
        // no slip => define boundary s.t. the interpolated value is zero
        v->Cell(it) = - v->Cell(it.Left());
    }

    it.SetBoundary(Boundary::Top);
    for(it.First(); it.Valid(); it.Next()) {
        // velocity given
        v->Cell(it) = this->_velocity[1]; // TODO what happens with this point?
        v->Cell(it.Down()) = this->_velocity[1];
    }
}

void Geometry::Update_P(Grid *p) const {
    BoundaryIterator it (this);

    it.SetBoundary(Boundary::Bottom);
    for(it.First(); it.Valid(); it.Next()) {
        // after choosong F_0,j := u_0,j we get homogeneous Neumann boundary conditions
        p->Cell(it) = p->Cell(it.Top());
    }

    it.SetBoundary(Boundary::Left);
    for(it.First(); it.Valid(); it.Next()) {
        // after choosong F_0,j := u_0,j we get homogeneous Neumann boundary conditions
        p->Cell(it) = p->Cell(it.Right());
    }

    it.SetBoundary(Boundary::Right);
    for(it.First(); it.Valid(); it.Next()) {
        // after choosong F_0,j := u_0,j we get homogeneous Neumann boundary conditions
        p->Cell(it) = p->Cell(it.Left());
    }

    it.SetBoundary(Boundary::Top);
    for(it.First(); it.Valid(); it.Next()) {
        // after choosong F_0,j := u_0,j we get homogeneous Neumann boundary conditions
        p->Cell(it) = p->Cell(it.Down());

    }
}
