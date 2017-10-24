#include "parameter.hpp"

#include <iostream>
#include <stdexcept>
#include <fstream>

/// \brief sets parameters for Driven Cavity 
Parameter::Parameter() {
    this->_re      = 1000;
    this->_omega   = 1.7;
    this->_tau     = 0.5;
    this->_eps     = 0.001;
    this->_alpha   = 0.9;
    this->_itermax = 100;
    this->_tend    = 16.4;
}

/// \brief load parameters from file
/// whitespace seperated list of:
/// re omega alpha dt tend eps tau itermax
void Parameter::Load(const char *file) {
    std::ifstream fin(file);

    if (fin.fail()) {
        throw std::runtime_error("Parameter::Load failed");
    } else {
        // TODO
        fin >> _re;
        fin >> _omega;
        fin >> _alpha;
        fin >> _dt;
        fin >> _tend;
        fin >> _eps;
        fin >> _tau;
        fin >> _itermax;
    }
}

const real_t &Parameter::Re() const {
    return this->_re;
}

const real_t &Parameter::Omega() const {
    return this->_omega;
}

const real_t &Parameter::Alpha() const {
    return this->_alpha;
}

const real_t &Parameter::Dt() const {
    return this->_dt;
}

const real_t &Parameter::Tend() const {
    return this->_tend;
}

const index_t &Parameter::IterMax() const {
    return this->_itermax;
}

const real_t &Parameter::Eps() const {
    return this->_eps;
}

const real_t &Parameter::Tau() const {
    return this->_tau;
}
