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
    _comm(comm),
    _flags(0) {

    this->split_for_comm();

    this->_flags = new FlagGrid(this);
}

Geometry::~Geometry() {
    delete this->_flags;
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
            fin >> param;
            if (param == "free") {
                this->split_for_comm();
                delete this->_flags;
                this->_flags = new FlagGrid(this);

                char x;
                fin.read(&x, 1); // Newline
                index_t stride_y = this->_size[0] + 2;
                for (int j = this->_size[1] + 2 - 1; j >= 0; j--) {
                    fin.read(&this->Flags().Data()[j * stride_y], stride_y);
                    fin.read(&x, 1); // Newline
                }
                // print
                for (int j = this->_size[1] + 2 - 1; j >= 0; j--) {
                    for (int i = 0; i <=  this->_size[0] + 2 - 1; i++) {
                        std::cout << Flags().Data()[j*stride_y + i];
                    }
                    std::cout << std::endl;
                }
                break;
            } else {
                std::cout << "parameter for 'geometry' differs from 'free'. Canceling..." << std::endl;
                break;
            }
        } else {
            std::cout << '|' << param << '|' << std::endl;
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

const FlagGrid &Geometry::Flags() const {
    return *this->_flags;
}

FlagGrid &Geometry::Flags() {
    return *this->_flags;
}

void Geometry::Update_U(Grid *u) const {
    /*
          |         |         |         |
        --+----v----+----v----+----v----+--
          |         |         |         |
          u    p    u    p    u    p    u
          |         |         |         |
        --+----v----+====V====+----v----+--
          |         ||       ||         |
          u    p    u    P    U    p    u
          |         ||       ||         |
        --+----v----+====v====+----v----+--
          |         |         |         |
          u    p    u    p    u    p    u
          |         |         |         |
        --+----v----+----v----+----v----+--
          |         |         |         |
    */
    Iterator it (this);
    for(it.First(); it.Valid(); it.Next()) {
        auto b_type = this->Flags().Cell(it);
        if (b_type == Flags::Fluid) continue;

        auto b_ori  = this->Flags().BoundaryOrientation(it);
        switch (b_type) {
        case Flags::Fluid:
            // nothing to do
            continue; break; // :P
        case Flags::Noslip:
            // set u to zero
            if (b_ori & BoundaryOrientation::Vert_l) {
                // global left boundary => left.flag = self.flag => ori != vert_l -> OK
                // left = fluid => we have to set the left boundary
                // left = boundary => it doesent matter what we do...
                u->Cell(it.Left()) = 0.0;
            } else if (b_ori & BoundaryOrientation::Vert_r) {
                u->Cell(it) = 0.0;
            } else if (b_ori & BoundaryOrientation::Horiz_t) {
                u->Cell(it) = -u->Cell(it.Top());
            } else if (b_ori & BoundaryOrientation::Horiz_b) {
                u->Cell(it) = -u->Cell(it.Down());
            }
            break;
        case Flags::Inflow:
        case Flags::InflowHorizontal:
        case Flags::InflowVertical: // TODO parabolo
        case Flags::SlipHorizontal: // TODO parabola?
            if (b_ori & BoundaryOrientation::Vert_l) {
                u->Cell(it.Left()) = this->_velocity[0];
            } else if (b_ori & BoundaryOrientation::Vert_r) {
                u->Cell(it) = this->_velocity[0];
            } else if (b_ori & BoundaryOrientation::Horiz_t) {
                u->Cell(it) = 2 * this->_velocity[0] - u->Cell(it.Top());
            } else if (b_ori & BoundaryOrientation::Horiz_b) {
                u->Cell(it) = 2 * this->_velocity[0] - u->Cell(it.Down());
            }
            break;
        case Flags::Outflow:
        case Flags::SlipVertical:
            // d/dn (u v)^T = 0
            if (b_ori & BoundaryOrientation::Vert_l) {
                u->Cell(it.Left()) = u->Cell(it.Left().Left());
                u->Cell(it)        = u->Cell(it.Left().Left()); // nice picture
            } else if (b_ori & BoundaryOrientation::Vert_r) {
                u->Cell(it) = u->Cell(it.Right());
            } else if (b_ori & BoundaryOrientation::Horiz_t) {
                u->Cell(it) = 2 * u->Cell(it.Top()) - u->Cell(it.Top());
            } else if (b_ori & BoundaryOrientation::Horiz_b) {
                u->Cell(it) = 2 * u->Cell(it.Down()) - u->Cell(it.Down());
            }
            break;
        }
    }
}

void Geometry::Update_V(Grid *v) const {
    Iterator it (this);
    for(it.First(); it.Valid(); it.Next()) {
        auto b_type = this->Flags().Cell(it);
        if (b_type == Flags::Fluid) continue;

        auto b_ori  = this->Flags().BoundaryOrientation(it);
        switch (b_type) {
        case Flags::Fluid:
            // nothing to do
            continue; break; // :P
        case Flags::Noslip:
            // set v to zero
            if (b_ori & BoundaryOrientation::Horiz_b) {
                v->Cell(it.Down()) = 0.0;
            } else if (b_ori & BoundaryOrientation::Horiz_t) {
                v->Cell(it) = 0.0;
            } else if (b_ori & BoundaryOrientation::Vert_l) {
                v->Cell(it) = -v->Cell(it.Left());
            } else if (b_ori & BoundaryOrientation::Vert_r) {
                v->Cell(it) = -v->Cell(it.Right());
            }
            break;
        case Flags::Inflow:
        case Flags::InflowHorizontal: // TODO parabola
        case Flags::InflowVertical:
        case Flags::SlipVertical: // TODO parabola?
            if (b_ori & BoundaryOrientation::Horiz_b) {
                v->Cell(it.Down()) = this->_velocity[1];
            } else if (b_ori & BoundaryOrientation::Horiz_t) {
                v->Cell(it) = this->_velocity[1];
            } else if (b_ori & BoundaryOrientation::Vert_l) {
                v->Cell(it) = 2 * this->_velocity[1] - v->Cell(it.Left());
            } else if (b_ori & BoundaryOrientation::Vert_r) {
                v->Cell(it) = 2 * this->_velocity[1] - v->Cell(it.Right());
            }
            break;
        case Flags::Outflow:
        case Flags::SlipHorizontal:
            // d/dn (u v)^T = 0
            if (b_ori & BoundaryOrientation::Horiz_b) {
                v->Cell(it.Down()) = v->Cell(it.Down().Down());
                v->Cell(it)        = v->Cell(it.Down().Down()); // nice picture
            } else if (b_ori & BoundaryOrientation::Horiz_t) {
                v->Cell(it) = v->Cell(it.Top());
            } else if (b_ori & BoundaryOrientation::Vert_l) {
                v->Cell(it) = 2 * v->Cell(it.Left()) - v->Cell(it.Left());
            } else if (b_ori & BoundaryOrientation::Vert_r) {
                v->Cell(it) = 2 * v->Cell(it.Right()) - v->Cell(it.Right());
            }
            break;
        }
    }
}

void Geometry::Update_P(Grid *p) const {
    Iterator it (this);
    for(it.First(); it.Valid(); it.Next()) {
        auto b_type = this->Flags().Cell(it);
        if (b_type == Flags::Fluid) continue;

        auto b_ori  = this->Flags().BoundaryOrientation(it);
        switch (b_type) {
        case Flags::Fluid:
            // nothing to do
            continue; break; // :P
        case Flags::Noslip:
        case Flags::Inflow:
        case Flags::InflowHorizontal:
        case Flags::InflowVertical:
            // just choose one side
            if (b_ori & BoundaryOrientation::Horiz_b) {
                p->Cell(it) = p->Cell(it.Down());
            } else if (b_ori & BoundaryOrientation::Horiz_t) {
                p->Cell(it) = p->Cell(it.Top());
            } else if (b_ori & BoundaryOrientation::Vert_l) {
                p->Cell(it) = p->Cell(it.Left());
            } else if (b_ori & BoundaryOrientation::Vert_r) {
                p->Cell(it) = p->Cell(it.Right());
            }
            break;
        case Flags::Outflow:
            if (b_ori & BoundaryOrientation::Horiz_b) {
                p->Cell(it) = -p->Cell(it.Down());
            } else if (b_ori & BoundaryOrientation::Horiz_t) {
                p->Cell(it) = -p->Cell(it.Top());
            } else if (b_ori & BoundaryOrientation::Vert_l) {
                p->Cell(it) = -p->Cell(it.Left());
            } else if (b_ori & BoundaryOrientation::Vert_r) {
                p->Cell(it) = -p->Cell(it.Right());
            }
            break;
        case Flags::SlipHorizontal:
        case Flags::SlipVertical:
            // set pressure
            if (b_ori & BoundaryOrientation::Horiz_b) {
                p->Cell(it) = 2 * this->_pressure - p->Cell(it.Down());
            } else if (b_ori & BoundaryOrientation::Horiz_t) {
                p->Cell(it) = 2 * this->_pressure - p->Cell(it.Top());
            } else if (b_ori & BoundaryOrientation::Vert_l) {
                p->Cell(it) = 2 * this->_pressure - p->Cell(it.Left());
            } else if (b_ori & BoundaryOrientation::Vert_r) {
                p->Cell(it) = 2 * this->_pressure - p->Cell(it.Right());
            }
        }
    }
}
