#include "geometry.hpp"
#include "grid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
#include "communicator.hpp"
#include "typedef.hpp"

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <algorithm>

#define N (12)

/// \brief sets geometry for Driven Cavity
Geometry::Geometry() : Geometry(0) {

}

Geometry::Geometry(const Communicator *comm) :
    _size(N, N),
    _length(1.0, 1.0),
    _h(1.0/(N + 1), 1.0/(N + 1)),
    // velocity at upper boundary
    _velocity(1.0, 0.0),
    _pressure(0.1),
    _comm(comm) {

    this->split_for_comm();
}

void Geometry::split_for_comm() {

    // n cells
    // |_ n/p _| per p
    // n % p left \in 0...p
    // if .%. != 0: distribute the rest over the processes
    index_t block_size_x = _size[0] / _comm->ThreadDim()[0];
    if (_size[0] % _comm->ThreadDim()[0] != 0) {
        block_size_x += 1;
    }
    index_t block_size_y = _size[1] / _comm->ThreadDim()[1];
    if (_size[1] % _comm->ThreadDim()[1] != 0) {
        block_size_y += 1;
    }


    if (this->_comm->isRight() && _size[0] % block_size_x != 0) {
        _bsize[0] = _size[0] % block_size_x;
        _blength[0] = _length[0] / _comm->ThreadDim()[0] * _bsize[0] / block_size_x;
    } else {
        _bsize[0] = _size[0] / _comm->ThreadDim()[0];
        _blength[0] = _length[0] / _comm->ThreadDim()[0];
    }
    if (this->_comm->isTop() && _size[1] % block_size_y != 0) {
        _bsize[1] = _size[1] % block_size_y;
        _blength[1] = _length[1] / _comm->ThreadDim()[1] * _bsize[1] / block_size_y;
    } else {
        _bsize[1] = _size[1] / _comm->ThreadDim()[1];
        _blength[1] = _length[1] / _comm->ThreadDim()[1];
    }
}

/// \brief load geometry settings form file
void Geometry::Load(const char *file) {
    if (this->_comm->getRank() == 0) std::cout << "Loading geometry from " << file << std::endl;
    std::ifstream fin (file);
    std::string param;
    while (fin >> param) {
        if (param == "size") {
            fin >> this->_size[0] >> this->_size[1];
            this->_h[0] = this->_length[0] / (this->_size[0] + 1);
            this->_h[1] = this->_length[1] / (this->_size[1] + 1);
        } else if (param == "length") {
            fin >> this->_length[0] >> this->_length[1];
            this->_h[0] = this->_length[0] / (this->_size[0] + 1);
            this->_h[1] = this->_length[1] / (this->_size[1] + 1);
        } else if (param == "velocity") {
            fin >> this->_velocity[0] >> this->_velocity[1];
        } else if (param == "pressure") {
            fin >> this->_pressure;
        } else if (param == "geometry") {
            std::cout << "'geometry' parameter not yet supported. Skipping..." << std::endl;
            break;
        } else {
            throw std::runtime_error("unsupported parameter in geometry file");
        }
    }
    if (this->_comm->getRank() == 0) {
        std::cout << "  size:     " << this->_size[0] << " " << this->_size[1] << std::endl;
        std::cout << "  length:   " << this->_length[0] << " " << this->_length[1] << std::endl;
        std::cout << "  -> h:     " << this->_h[0] << " " << this->_h[1] << std::endl;
        std::cout << "  velocity: " << this->_velocity[0] << " " << this->_velocity[1] <<std::endl;
        std::cout << "  pressure: " << this->_pressure << std::endl;
        std::cout << "  geometry: " << "<skipped>" << std::endl;
    }

    this->split_for_comm();
}


const multi_index_t &Geometry::TotalSize() const {
    return this->_size;
}

const multi_real_t &Geometry::TotalLength() const {
    return this->_length;
}

const multi_index_t &Geometry::Size() const {
    return this->_bsize;
}

const multi_real_t &Geometry::Length() const {
    return this->_blength;
}

const multi_real_t &Geometry::Mesh() const {
    return this->_h;
}


void Geometry::Update_U(Grid *u) const {
    BoundaryIterator it (this);

    if (this->_comm->isBottom()) {
        it.SetBoundary(Boundary::Bottom);
        for(it.First(); it.Valid(); it.Next()) {
            // no slip => define boundary s.t. the interpolated value is zero
            u->Cell(it) = - u->Cell(it.Top());
        }
    }

    if (this->_comm->isLeft()) {
        it.SetBoundary(Boundary::Left);
        for(it.First(); it.Valid(); it.Next()) {
            // no slip
            u->Cell(it) = 0.0;
        }
    }

    if (this->_comm->isRight()) {
        it.SetBoundary(Boundary::Right);
        for(it.First(); it.Valid(); it.Next()) {
            // no slip, staggered grid => set left (= boundary) point to zero
            u->Cell(it) = 0.0;
            u->Cell(it.Left()) = 0.0;
        }
    }

    if (this->_comm->isTop()) {
        it.SetBoundary(Boundary::Top);
        for(it.First(); it.Valid(); it.Next()) {
            // velocity given
            u->Cell(it) = 2* this->_velocity[0] - u->Cell(it.Down());
        }
    }
}

void Geometry::Update_V(Grid *v) const {
    BoundaryIterator it (this);

    if (this->_comm->isBottom()) {
        it.SetBoundary(Boundary::Bottom);
        for(it.First(); it.Valid(); it.Next()) {
            // no slip
            v->Cell(it) = 0.0;
        }
    }

    if (this->_comm->isLeft()) {
        it.SetBoundary(Boundary::Left);
        for(it.First(); it.Valid(); it.Next()) {
            // no slip => define boundary s.t. the interpolated value is zero
            v->Cell(it) = - v->Cell(it.Right());
        }
    }

    if (this->_comm->isRight()) {
        it.SetBoundary(Boundary::Right);
        for(it.First(); it.Valid(); it.Next()) {
            // no slip => define boundary s.t. the interpolated value is zero
            v->Cell(it) = - v->Cell(it.Left());
        }
    }

    if (this->_comm->isTop()) {
        it.SetBoundary(Boundary::Top);
        for(it.First(); it.Valid(); it.Next()) {
            // velocity given
            v->Cell(it) = this->_velocity[1];
            v->Cell(it.Down()) = this->_velocity[1];
        }
    }
}

void Geometry::Update_P(Grid *p) const {
    BoundaryIterator it (this);

    if (this->_comm->isBottom()) {
        it.SetBoundary(Boundary::Bottom);
        for(it.First(); it.Valid(); it.Next()) {
            // after choosong F_0,j := u_0,j we get homogeneous Neumann boundary conditions
            p->Cell(it) = p->Cell(it.Top());
        }
    }

    if (this->_comm->isLeft()) {
        it.SetBoundary(Boundary::Left);
        for(it.First(); it.Valid(); it.Next()) {
            // after choosong F_0,j := u_0,j we get homogeneous Neumann boundary conditions
            p->Cell(it) = p->Cell(it.Right());
        }
    }

    if (this->_comm->isRight()) {
        it.SetBoundary(Boundary::Right);
        for(it.First(); it.Valid(); it.Next()) {
            // after choosong F_0,j := u_0,j we get homogeneous Neumann boundary conditions
            p->Cell(it) = p->Cell(it.Left());
        }
    }

    if (this->_comm->isTop()) {
        it.SetBoundary(Boundary::Top);
        for(it.First(); it.Valid(); it.Next()) {
            // after choosong F_0,j := u_0,j we get homogeneous Neumann boundary conditions
            p->Cell(it) = p->Cell(it.Down());
        }
    }
}
