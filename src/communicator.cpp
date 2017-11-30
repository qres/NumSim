#include <mpi/mpi.h> # keep this at the top. Otherwise we get problems with MPI_REAL

#include "communicator.hpp"
#include "geometry.hpp"
#include <iostream>
#include "typedef.hpp"

Communicator::Communicator(int *argc, char ***argv) {
    MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &this->_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &this->_size);

    // compute based on _size
    this->_tdim = multi_index_t(1,1); // TODO
    this->_tidx = multi_index_t(_rank % this->_tdim[0], _rank / this->_tdim[0]);
    this->_evenodd = (_tidx[0] + _tidx[1]) % 2 == 0;
}

Communicator::~Communicator() {
    MPI_Finalize();
}

/** Returns the position of the current process with respect to the
 *  fields lower left corner
 */
const multi_index_t& Communicator::ThreadIdx() const {
    return this->_tidx;
}

/** Returns the way the domain is partitioned among all processes
 */
const multi_index_t& Communicator::ThreadDim() const {
    return this->_tdim;
}

/** Returns whether this process is a red or a black field
 */
const bool& Communicator::EvenOdd() const {
    return this->_evenodd;
}

/** Gets the sum of all values and distributes the result among all
 *  processes
 *
 * \param [in] val The data over which the sum is to be calculated
 */
real_t Communicator::gatherSum(const real_t &val) const {
    real_t reduced;
    MPI_Allreduce(&val, &reduced, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
    return reduced;
}

/** Finds the minimum of the values and distributes the result among
 *  all processes
 *
 * \param [in] val The data over which to find the minimum
 */
real_t Communicator::gatherMin(const real_t &val) const {
    real_t reduced;
    MPI_Allreduce(&val, &reduced, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD);
    return reduced;
}

/** Finds the maximum of the values and distributes the result among
 *  all processes
 *
 * \param [in] val The data over which to find the maximum
 */
real_t Communicator::gatherMax(const real_t &val) const {
    real_t reduced;
    MPI_Allreduce(&val, &reduced, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD);
    return reduced;
}

/** Synchronizes ghost layer
 *
 * \param [in] grid  The values to sync
 */
void Communicator::copyBoundary(Grid *grid) const {
    // exchange every second vertical boundary
    const multi_index_t grid_size = grid->getGeometry()->Size();
    const index_t Nx = grid_size[0];
    const index_t Ny = grid_size[1];
    MPI_Datatype grid_row;
    MPI_Datatype grid_col;
    const multi_index_t grid_stride(1, Nx+2);
    // MPI_Type_vector(number of blocks, elems per block, stride for blocks, type in block, new_type)
    //MPI_Type_vector(1, grid_size[0], 0, MPI_REAL, &grid_row); alternative as the block is contineous
    MPI_Type_vector(grid_size[0] + 2, 1, grid_stride[0], MPI_REAL, &grid_row);
    MPI_Type_vector(grid_size[1] + 2, 1, grid_stride[1], MPI_REAL, &grid_col);
    MPI_Type_commit(&grid_row);
    MPI_Type_commit(&grid_col);

    real_t* data = grid->Data();

    /*
        Version 1:
            copy all data to the right
            copy all data to the left
            copy all data up
            copy all data down
        Version 2:
            swap the even-uneven boundary horizontally
            swap the even-uneven boundary vertically
            swap the uneven-even boundary horizontally
            swap the uneven-even boundary vertically
    */
    const bool BOUNDARY_SWEEP = true;
    const bool BOUNDARY_SWAP  = !BOUNDARY_SWEEP;

    if (BOUNDARY_SWEEP) {
        if (ThreadDim()[0] > 1) {
            // copy to the right, copy to the left
            if (isLeft()) {
                // LEFT
                MPI_Send(&data[Nx],   1, grid_col, getRank() + 1, 1, MPI_COMM_WORLD);                    // send N -> right
                MPI_Recv(&data[Nx+1], 1, grid_col, getRank() + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recive right -> N+1
            } else if (isRight()) {
                // RIGHT
                MPI_Recv(&data[0], 1, grid_col, getRank() - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recive left -> 0
                MPI_Send(&data[1], 1, grid_col, getRank() - 1, 1, MPI_COMM_WORLD);                    // send 1 -> left
            } else {
                // INNER
                MPI_Sendrecv(
                    &data[Nx], 1, grid_col, getRank() + 1, 1, // send N -> right
                    &data[0],  1, grid_col, getRank() - 1, 1, // recive left -> 0
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                );
                MPI_Sendrecv(
                    &data[1],    1, grid_col, getRank() - 1, 1, // send 1 -> left
                    &data[Nx+1], 1, grid_col, getRank() + 1, 1, // recive right -> N+1
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                );
            }
        }
        // copy to the top, copy down
        if (ThreadDim()[1] > 1) {
            const index_t stride_y = grid_size[0];
            const index_t stride_ry = ThreadDim()[0];
            const index_t row0   = 0;
            const index_t row1   = stride_y;
            const index_t rowN   = Ny*stride_y;
            const index_t rowNp1 = (Ny+1)*stride_y;
            if (isTop()) {
                // TOP
                MPI_Send(&data[rowN],   1, grid_row, getRank() + stride_ry, 1, MPI_COMM_WORLD);                    // send N -> top
                MPI_Recv(&data[rowNp1], 1, grid_row, getRank() + stride_ry, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recive top -> N+1
            } else if (isBottom()) {
                // BOTTOM
                MPI_Recv(&data[row0], 1, grid_row, getRank() - stride_ry, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recive down -> 0
                MPI_Send(&data[row1], 1, grid_row, getRank() - stride_ry, 1, MPI_COMM_WORLD);                    // send 1 -> down
            } else {
                // INNER
                MPI_Sendrecv(
                    &data[rowN], 1, grid_row, getRank() + stride_ry, 1, // send N -> top
                    &data[row0], 1, grid_row, getRank() - stride_ry, 1, // recive down -> 0
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                );
                MPI_Sendrecv(
                    &data[row1],   1, grid_row, getRank() - stride_ry, 1, // send 1 -> down
                    &data[rowN+1], 1, grid_row, getRank() + stride_ry, 1, // recive top -> N+1
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                );
            }
        }
    } else if (BOUNDARY_SWAP) {
        // swap right side of Even and top of Even
        /*if (EvenOdd()) {
            if (!isRight()) {
                send N-1 -> right, recive right -> N
            }
            if (!isTop()) {
                send N-1 -> top, recive top -> N
            }
        } else {
            if (!isLeft()) {
                send 1 -> left, recive left -> 0
            }
            if (!isBottom()) {
                send 1 -> down, recive down -> 0
            }
        }
        // swap left side of Even and down of Even
        if (EvenOdd()) {
            if (!isLeft()) {
                send 1 -> left, recive left -> 0
            }
            if (!isBottom()) {
                send 1 -> down, recive down -> 0
            }
        } else (!EvenOdd() && !isRight() {
            if (!isRight()) {
                send N-1 -> right, recive right -> N
            }
            if (!isTop()) {
                send N-1 -> top, recive top -> N
            }
        }*/
    } else {
        std::cout << "ERR: no boundary excange defined" << std::endl;
    }

    MPI_Type_free(&grid_row);
    MPI_Type_free(&grid_col);
}


/** Decide whether our left boundary is a domain boundary
*/
const bool Communicator::isLeft() const {
    return this->_tidx[0] == 0;
}

/** Decide whether our right boundary is a domain boundary
*/
const bool Communicator::isRight() const {
    return this->_tidx[0] == this->_tdim[0] - 1;
}

/** Decide whether our top boundary is a domain boundary
*/
const bool Communicator::isTop() const {
    return this->_tidx[1] == 0;
}

/** Decide whether our bottom boundary is a domain boundary
*/
const bool Communicator::isBottom() const {
    return this->_tidx[1] == this->_tdim[1] - 1;
}

/** Get MPI rank of current process
*/
const int &Communicator::getRank() const {
    return this->_rank;
}

/** Get number of MPI processes
*/
const int &Communicator::getSize() const {
    return this->_size;
}
