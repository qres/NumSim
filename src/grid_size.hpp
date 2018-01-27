#ifndef GRID_SIZE_HPP
#define GRID_SIZE_HPP

#include "typedef.hpp"

struct Grid2D {
    static unsigned int dim() {
        return 2;
    }

    static unsigned int size_N(multi_index_t N) {
        return (N[0]+2)*(N[1]+2);
    }

    static unsigned int size_n(multi_index_t N) {
        return (N[0]/2+2)*(N[1]/2+2);
    }

    static multi_index_t coarsen(multi_index_t N) {
        return multi_index_t(N[0]/2, N[1]/2);
    }
};

#endif
