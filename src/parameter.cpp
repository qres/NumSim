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
    this->_dt      = 0;
}

/// \brief load parameters from file
/// whitespace seperated list of:
/// re omega alpha dt tend eps tau itermax
void Parameter::Load(const char *file) {
    std::cout << "Loading parameters from " << file << std::endl;
    std::ifstream fin (file);
    std::string eq;
    std::string param;
    while (fin >> param >> eq) {
        if (param == "re") {
            fin >> this->_re;
        } else if (param == "omega") {
            fin >> this->_omega;
        } else if (param == "tau") {
            fin >> this->_tau;
        } else if (param == "eps") {
            fin >> this->_eps;
        } else if (param == "alpha") {
            fin >> this->_alpha;
        } else if (param == "iter") {
            fin >> this->_itermax;
        } else if (param == "tend") {
            fin >> this->_tend;
        } else if (param == "dt") {
            fin >> this->_dt;
        } else {
            std::cout << param << std::endl;
            throw std::runtime_error("unsupported parameter in parameter file");
        }
    }
    std::cout << "  re:    " << this->_re << std::endl;
    std::cout << "  omega: " << this->_omega << std::endl;
    std::cout << "  tau:   " << this->_tau << std::endl;
    std::cout << "  eps:   " << this->_eps << std::endl;
    std::cout << "  alpha: " << this->_alpha << std::endl;
    std::cout << "  iter:  " << this->_itermax << std::endl;
    std::cout << "  tend:  " << this->_tend << std::endl;
    std::cout << "  dt:    " << this->_dt << std::endl;
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
