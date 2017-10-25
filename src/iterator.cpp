#include "iterator.hpp"
#include "geometry.hpp"
#include "typedef.hpp"

#include <iostream>
#include <stdexcept>

/// Constructs a new Iterator depending on a geometry
Iterator::Iterator(const Geometry *geom) {
    this->_geom = geom;
    this->First();
}

/// Constructs a new Iterator on a geometry with a defined starting value
Iterator::Iterator(const Geometry *geom, const index_t &value) {
    this->_geom = geom;
    this->_value = value;
    this->_valid = true;
}

///     Returns the current position value
const index_t &Iterator::Value() const {
    return this->_value;
}
/// Cast operator to convert Iterators to integers
Iterator::operator const index_t &() const {
    return this->Value();
}
/// Returns the position coordinates
multi_index_t Iterator::Pos() const {
    multi_index_t size = this->_geom->Size();
    size_t ix = this->Value() % (size[0] + 2); //x-minor
    size_t jy = this->Value() / (size[0] + 2); //y-major
    
    //multi_real_t offset = this->_geom->_offset; // TODO if we cannot access this here... whats the offset for?
    multi_real_t h = this->_geom->Mesh();
    return multi_index_t(/*offset[0] + */ix * h[0], /*offset[1] + */jy * h[1]);
}

/// Sets the iterator to the first element
void Iterator::First() {
    this->_value = 0;
    this->_valid = true;
}
/// Goes to the next element of the iterator, disables it if position is end
void Iterator::Next() {
    multi_index_t size = this->_geom->Size();
    if (this->_value < (size[0] + 2)*(size[1] + 2) - 1) {
        this->_value += 1;
    } else {
        this->_valid = false;
    }
}

/// Checks if the iterator still has a valid value
bool Iterator::Valid() const {
    return this->_valid;
}

/// Returns an Iterator that is located left from this one.
// if we are at the left boundary, the cell sees itself
Iterator Iterator::Left() const {
    multi_index_t size = this->_geom->Size();
    size_t ix = this->Value() % (size[0] + 2); //x-minor
    if (ix != 0) {
        return Iterator(this->_geom, this->_value - 1);
    } else {
        // we are on a boundary
        return *this; // TODO valid flag?
    }
}

/// Returns an Iterator that is located right from this one
// If we are at the right boundary, the cell sees itself
Iterator Iterator::Right() const {
    multi_index_t size = this->_geom->Size();
    size_t ix = this->Value() % (size[0] + 2); //x-minor
    if (ix != this->_geom->Size()[0] + 1) { // + boundary
        return Iterator(this->_geom, this->_value + 1);
    } else {
        // we are on a boundary
        return *this; // TODO valid flag?
    }
}


/// Returns an Iterator that is located below this one
// If we are at the lower domain boundary, the cell sees itself
Iterator Iterator::Down() const {
    multi_index_t size = this->_geom->Size();
    index_t stride_y = size[0] + 2;
    size_t jy = this->Value() / stride_y; //y-major
    if (jy != 0) {
        return Iterator(this->_geom, this->_value - stride_y);
    } else {
        // we are on a boundary
        return *this; // TODO valid flag?
    }
}

/// Returns an Iterator that is located above this one
// If we are at the upper domain boundary, the cell sees itself
Iterator Iterator::Top() const {
    multi_index_t size = this->_geom->Size();
    index_t stride_y = size[0] + 2;
    size_t jy = this->Value() / stride_y; //y-major
    if (jy != this->_geom->Size()[1] + 1) {
        return Iterator(this->_geom, this->_value + stride_y);
    } else {
        // we are on a boundary
        return *this; // TODO valid flag?
    }
}

InteriorIterator::InteriorIterator(const Geometry *geom) : Iterator(geom) {
    
}

/// Sets the iterator to the first element
void InteriorIterator::First() {
    multi_index_t size = this->_geom->Size();
    index_t stride_y = size[0] + 2;
    // go to (1,1) as (0,0) is boundary
    this->_value = 1 + stride_y * 1;
    this->_valid = true;
}
/// Goes to the next element of the iterator, disables it if position is end
void InteriorIterator::Next() {
    multi_index_t size = this->_geom->Size();
    index_t stride_y = size[0] + 2;
    this->_value += 1;
    size_t ix = this->Value() % stride_y; //x-minor
    size_t jy = this->Value() / stride_y; //y-major
    if (ix == size[0] + 1) {
        // skip boundary
        this->_value += 2;
        if (jy == size[1]) {
            // we are at the last cell
            this->_valid = false;
        }
    }
}

/// Constructs a new BoundaryIterator for the bottom
BoundaryIterator::BoundaryIterator(const Geometry *geom) : Iterator(geom) {
    this->_boundary = Boundary::Bottom;
}

/// Sets the boundary to iterate
/// \param boundary 0: bottom, 1: left, 2: right, 3: top
void BoundaryIterator::SetBoundary(const index_t &boundary) {
    if (boundary >= Boundary::_MaxBoundary) {
        throw std::runtime_error("Invalid boundary");
    }

    this->_boundary = boundary;
}

/// Sets the iterator to the first element
void BoundaryIterator::First() {
    multi_index_t size = this->_geom->Size();
    index_t stride_y = size[0] + 2;

    switch (this->_boundary) {
    case Boundary::Bottom:
        this->_value = 1;
        break;
    case Boundary::Left:
        this->_value = stride_y;
        break;
    case Boundary::Right:
        this->_value = 2*stride_y - 1;
        break;
    case Boundary::Top:
        this->_value = 1 + (size[1] + 1)*stride_y;
        break;
    default:
        throw std::runtime_error("Unreachable");
    }
    this->_valid = true;
}
/// Goes to the next element of the iterator, disables it if position is end
void BoundaryIterator::Next() {
    multi_index_t size = this->_geom->Size();
    index_t stride_y = size[0] + 2;

    switch (this->_boundary) {
    case Boundary::Bottom:
        if (this->_value < size[0]) {
            this->_value += 1;
        } else {
            this->_valid = false;
        }
        break;
    case Boundary::Left:
        if (this->_value < size[1] * stride_y) {
            this->_value += stride_y;
        } else {
            this->_valid = false;
        }
        break;
    case Boundary::Right:
        if (this->_value < size[0] + 1 + size[1] * stride_y) {
            this->_value += stride_y;
        } else {
            this->_valid = false;
        }
        break;
    case Boundary::Top:
        if (this->_value < size[0] + (size[1] + 1) * stride_y) {
            this->_value += 1;
        } else {
            this->_valid = false;
        }
        break;
    default:
        throw std::runtime_error("Unreachable");
    }
}

