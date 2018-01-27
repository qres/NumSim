#include "grid_size.hpp"

#include <iostream>

void checkForErrors(const cudaError_t status, const char *label, const int line, const char *file)
{
    if (status != cudaSuccess) {
        std::cerr << "CUDA ERROR (" << label << ") ";
        std::cerr << "at " << file << ":" << line << std::endl;
        std::cerr << cudaGetErrorString(status) << ". Exiting..." << std::endl;
        exit(1);
    }
}

//#define CUDA_CHECK_KERNEL_ERRORS

#define cuchck(func) checkForErrors(func, #func, __LINE__, __FILE__)
#ifdef CUDA_CHECK_KERNEL_ERRORS
    #warning !!! !!! !!! CUDA SYNC IS ACTIVE !!! !!! !!!
    #define cuchck_last() cudaDeviceSynchronize(); cuchck(cudaGetLastError())
#else
    #define cuchck_last()
#endif

float normoo_cu(const float* vec, unsigned int count) {
    // TODO
    return abs(1);
}

double normoo_cu(const double* vec, unsigned int count) {
    // TODO
    return abs(1);
}

float norm2_cu(const float* vec, unsigned int count) {
    return 1; // TODO
}

double norm2_cu(const double* vec, unsigned int count) {
    return 1; // TODO
}

template<typename T>
T norm_cu(const T* vec, unsigned int count) {
    return normoo_cu(vec, count);
}

namespace cuda {

    template<typename T>
    __global__
    void add(T* dst, const T* lhs, const T* rhs, const unsigned int count) {
        const int i (threadIdx.x + blockDim.x * blockIdx.x);
        if (i < count) {
            dst[i] = lhs[i] + rhs[i];
        }
    }

    template<typename T>
    __global__
    void sub(T* dst, const T* lhs, const T* rhs, const unsigned int count) {
        const int i (threadIdx.x + blockDim.x * blockIdx.x);
        if (i < count) {
            dst[i] = lhs[i] - rhs[i];
        }
    }
}

namespace cuda {
    template<typename T>
    __global__
    void zerof64(T* dst, const unsigned int count) {
        const int i (threadIdx.x + blockDim.x * blockIdx.x);
        if (i < count) {
            dst[i] = 0.0;
        }
    }

    template<typename T>
    __global__
    void addAssign(T* dst, const T* vec, const unsigned int count) {
        const int i (threadIdx.x + blockDim.x * blockIdx.x);
        if (i < count) {
            dst[i] += vec[i];
        }
    }

    //http://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/reduction/doc/reduction.pdf
    void reduce_sum() {
        std::cout << "ERR reduce sum" << std::endl;
    }
}


namespace cuda {
    /** x,y are block indices, i,j are thread indices in block x,y
     *  storage is row major: 00 01 02 03 ... 10 11 12 13 ... ...
     *
     *   o--> y,j
     *   |
     *   V x,i
     *
     *      +-------+-------+  \     ~.
     *      | 00 01 | 02 03 |   | Bx  |
     *      | 10 11 | 12 13 |   |     |
     *      +-------+-------+  /      | Gx
     *      | 20 21 | 22 23 |         |
     *      | 30 31 | 32 33 |         |
     *      +-------+-------+        ~'
     *
     *      \__By___/
     *
     *      ^*-----Gy------*^
     *
     *
     */
    template<typename T>
    __global__
    void residuum_neg_laplace(const T* u, const T* b, T* res, multi_index_t N, multi_real_t length) {
        const int i (threadIdx.x + blockIdx.x * blockDim.x);
        const int j (threadIdx.y + blockIdx.y * blockDim.y);
        const int ix (i + j * (N[0] + 2));

        if ((i <= N[0]+1 && j <= N[1]+1) && (i == 0 || j == 0 || i == N[0]+1 || j == N[1]+1)) {
            res[ix] = 0;
        } else if (i <= N[0] && j <= N[1]) {
            const int im1_j = ix - 1;
            const int ip1_j = ix + 1;
            const int i_jm1 = ix - (N[0] + 2);
            const int i_jp1 = ix + (N[0] + 2);
            res[ix] = b[ix] - (N[0]+1)*(N[1]+1)/length[0]/length[1]*(4*u[ix] - 1*u[im1_j] - 1*u[ip1_j] - 1*u[i_jm1] - 1*u[i_jp1]);
        }
    }

    template<typename T>
    __global__
    void jacobi_step_neg_laplace(double omega, const T* u0, T* u1, const T* b, multi_index_t N, multi_real_t length) {
        const int i (threadIdx.x + blockIdx.x * blockDim.x);
        const int j (threadIdx.y + blockIdx.y * blockDim.y);
        const int ix (i + j * (N[0] + 2));

        if ((i <= N[0]+1 && j <= N[1]+1) && (i == 0 || j == 0 || i == N[0]+1 || j == N[1]+1)) {
            // boundary is known
            u1[ix] = u0[ix];
        } else if (i <= N[0] && j <= N[1]) {
            const int im1_j = ix - 1;
            const int ip1_j = ix + 1;
            const int i_jm1 = ix - (N[0] + 2);
            const int i_jp1 = ix + (N[0] + 2);
            const float u_ = 0.25 * (b[ix] / ((N[0]+1)*(N[1]+1)) * length[0] * length[1] + u0[im1_j] + u0[ip1_j] + u0[i_jm1] + u0[i_jp1]);
            u1[ix] = (1-omega) * u0[ix] + omega * u_;
        }
    }

    // treads: N*N, one thread per per inner grid point
    template<typename T>
    __global__
    void jacobi_step_neg_laplace_inner(double omega, const T* u0, T* u1, const T* b, multi_index_t N, multi_real_t length) {
        const int i (1 + threadIdx.x + blockIdx.x * blockDim.x);
        const int j (1 + threadIdx.y + blockIdx.y * blockDim.y);
        const int ix (i + j * (N[0] + 2));

        if (i <= N[0] && j <= N[1]) {
            const int im1_j = ix - 1;
            const int ip1_j = ix + 1;
            const int i_jm1 = ix - (N[0] + 2);
            const int i_jp1 = ix + (N[0] + 2);
            const float u_ = 0.25 * (b[ix] / ((N[0]+1)*(N[1]+1)) * length[0] * length[1] + u0[im1_j] + u0[ip1_j] + u0[i_jm1] + u0[i_jp1]);
            u1[ix] = (1-omega) * u0[ix] + omega * u_;
        }
    }

    template<typename T>
    __global__
    void restrict_fw_2D(const T* v_N, T* v_n, multi_index_t N, multi_index_t n) {
        // big indices for 'N'-Matrices, small indices for 'n'-Matrices
        const int i (threadIdx.x + blockIdx.x * blockDim.x);
        const int j (threadIdx.y + blockIdx.y * blockDim.y);
        const int ix (i + j * (n[0] + 2));
        const int I (2*i);
        const int J (2*j);
        const int IX (I + J * (N[0] + 2));

        /*
            O . O . O
            . . . . .
            O . O . O
            . . . . .
            O . O . O
        */

        if ((i <= n[0]+1 && j <= n[1]+1) && (i == 0 || j == 0 || i == n[0]+1 || j == n[1]+1)) {
            // boundary is known
            v_n[ix] = v_N[IX];
        } else if (i <= n[0] && j <= n[1]) {
            const int dI = 1;
            const int dJ = N[0]+2;
            v_n[ix] = 1.0/16.0 * (
                1*v_N[IX-dI-dJ] + 2*v_N[IX-dJ] + 1*v_N[IX+dI-dJ] +
                2*v_N[IX-dI   ] + 4*v_N[IX   ] + 2*v_N[IX+dI   ] +
                1*v_N[IX-dI+dJ] + 2*v_N[IX+dJ] + 1*v_N[IX+dI+dJ]
            );
        }
    }

    template<typename T>
    __global__
    void interpolate_2D(const T* v_n, T* v_N, multi_index_t n, multi_index_t N) {
        // big indices for 'N'-Matrices, small indices for 'n'-Matrices
        const int i (threadIdx.x + blockIdx.x * blockDim.x);
        const int j (threadIdx.y + blockIdx.y * blockDim.y);
        const int ix (i + j * (N[0] + 2));
        const int I (2*i);
        const int J (2*j);
        const int IX (I + J * (N[0] + 2));

        /*
            O . O . O
            . . . . .
            O . O . O
            . . . . .
            O . O . O
        */

        if ((i <= n[0]+1 && j <= n[1]+1) && (i == 0 || j == 0)) {
            // boundary is known
            // no update assuming that the boundary of v_N is already set
        } else if (i <= n[0] && j <= n[1]) {
            /* working on patches like this when all '.' are inner points
                . .
                . O
            */
            const int di = 1;
            const int dj = n[0] + 2;
            const int dI = 1;
            const int dJ = N[0] + 2;
            v_N[IX] = v_n[ix];
            v_N[IX - dI] = 0.5 * (v_n[ix - di] + v_n[ix]);
            v_N[IX - dJ] = 0.5 * (v_n[ix - dj] + v_n[ix]);
            v_N[IX - dI - dJ] = 0.25 * (v_n[ix - di] + v_n[ix - dj] + v_n[ix - di - dj] + v_n[ix]);
        } else if (i <= n[0]+1 && j <= n[1]+1) {
            /* working on patches like this when some '.' are bondary points

                . .
                . O
                ----

                . . |
                . O |

                . . |
                . O |
                ----+
            */
            const int di = 1;
            const int dj = n[0] + 2;
            const int dI = 1;
            const int dJ = N[0] + 2;
            // don't do v_N[IX] = v_n[ix]; as the point is on the border for sure -> no update
            if (J != N[1]+1) v_N[IX - dI] = 0.5 * (v_n[ix - di] + v_n[ix]);
            if (I != N[0]+1) v_N[IX - dJ] = 0.5 * (v_n[ix - dj] + v_n[ix]);
            v_N[IX - dI - dJ] = 0.25 * (v_n[ix - di] + v_n[ix - dj] + v_n[ix - di - dj] + v_n[ix]);
        }


    }

    // treads N*N: one thrad per inner grid point
    template<typename T>
    __global__
    void addAssign2Dinner(T* dst, const T* vec, multi_index_t N) {
        // index 0 is on the border
        const int i (1 + threadIdx.x + blockIdx.x * blockDim.x);
        const int j (1 + threadIdx.y + blockIdx.y * blockDim.y);
        const int ix (i + j * (N[0] + 2));
        if (i < N[0]+1 && j < N[1]+1) {
            dst[ix] += vec[ix];
        }
    }

}

template<typename T>
struct Fn_CUDA_mem {
    static bool mem_device_host_equal() {
        return false;
    }

    static T* malloc_typed(unsigned int count) {
        T* ptr = 0;
        cuchck(cudaMalloc(&ptr, sizeof(T) * count));
        return ptr;
    }

    static void free_typed(T* ptr) {
        cuchck(cudaFree(ptr));
    }

    static void memcpy_typed_HostToDevice(T* dst, const T* src, unsigned int count) {
        cuchck(cudaMemcpy(dst, src, sizeof(T) * count, cudaMemcpyHostToDevice));
    }

    static void memcpy_typed_DeviceToHost(T* dst, const T* src, unsigned int count) {
        cuchck(cudaMemcpy(dst, src, sizeof(T) * count, cudaMemcpyDeviceToHost));
    }

    static void memcpy_typed_DeviceToDevice(T* dst, const T* src, unsigned int count) {
        cuchck(cudaMemcpy(dst, src, sizeof(T) * count, cudaMemcpyDeviceToDevice));
    }

    static void zero(T* dst, const unsigned int count) {
        dim3 block(1024, 1, 1);
        dim3 grid((unsigned)ceil(count/(double)block.x), 1, 1);
        cuda::zerof64<<<grid, block>>>(dst, count);
        cuchck_last();
    }
};

template<typename T>
struct Fn_laplace_cuda : Fn_CUDA_mem<T>, Grid2D {
    static void restrict(multi_index_t N, const T* v_N, T* v_n) {
        const multi_index_t n = coarsen(N);
        dim3 block(1, 512, 1);
        dim3 grid((unsigned)ceil((N[0]+2)/(double)block.x), (unsigned)ceil((N[1]+2)/(double)block.y), 1);
        cuda::restrict_fw_2D<<<grid, block>>>(v_N, v_n, N, n);
        cuchck_last();
    }

    static void interpolate(multi_index_t n, const T* v_n, T* v_N, const char* mask) {
        const multi_index_t N (n[0]*2 + 1, n[1]*2 + 1);
        dim3 block(1, 512, 1);
        dim3 grid((unsigned)ceil((N[0]+2)/(double)block.x), (unsigned)ceil((N[1]+2)/(double)block.y), 1);
        cuda::interpolate_2D<<<grid, block>>>(v_n, v_N, n, N);
        cuchck_last();
    }

    static void residuum(multi_index_t N, multi_real_t length, const T* u, const T* b, const char* mask, T* res) {
        dim3 block(2, 64, 1);
        dim3 grid((unsigned)ceil((N[0]+2)/(double)block.x), (unsigned)ceil((N[1]+2)/(double)block.y), 1);
        cuda::residuum_neg_laplace<<<grid, block>>>(u, b, res, N, length);
        cuchck_last();
    }

    static T norm(multi_index_t N, const T* vec0) {
        return ::norm_cu(vec0, size_N(N));
    }

    static T norm_residuum(multi_index_t N, multi_real_t length, const T* u, const T* b, const char* mask, T* scratch) {
      residuum(N, length, u, b, mask, scratch);
      return norm_cu(scratch, size_N(N));
    }

    // TODO scratch buffer might be obsolate
    static unsigned int _jacobi(unsigned int max_iters, double max_r, double omega, T* u0, T* u1, const T* b, multi_index_t N, multi_real_t length, const char* mask, T* scratch) {
        // as we copy (only if we have to) the result back to u0, which has its
        // initialized, we can also only interate on the inner points, if
        // we initialize the u1 borders first

        dim3 block(1, 512, 1);
        dim3 grid_all  ((unsigned)ceil((N[0]+2)/(double)block.x), (unsigned)ceil((N[1]+2)/(double)block.y), 1);
        dim3 grid_inner((unsigned)ceil( N[0]   /(double)block.x), (unsigned)ceil( N[1]   /(double)block.y), 1);

        T * const dst = u0;
        unsigned int iters = 1;
        T r = norm_residuum(N, length, u0, b, mask, u1);
        while (iters <= max_iters && r >= max_r) {
            T r_old = r;

            // as mcpy will cause a sync, we spend one more block in the first iteration
            // to initialize the borders in u1
            if (iters == 1) {
                cuda::jacobi_step_neg_laplace<<<grid_all, block>>>(omega, u0, u1, b, N, length);
            } else {
                cuda::jacobi_step_neg_laplace_inner<<<grid_inner, block>>>(omega, u0, u1, b, N, length);
            }
            cuchck_last();


            // swap input and output
            std::swap(u0, u1);

            r = norm_residuum(N, length, u0, b, mask, scratch);
            ++iters;

            if (r == r_old) {
                if (max_r != 0) break; // max_res == 0 means no residuum checks
            }
        }
        if (dst != u0) {
            cudaMemcpy(dst, u0, size_N(N) * sizeof(T), cudaMemcpyDeviceToDevice);
        }
        return iters - 1;
    }

    static unsigned int solve(unsigned int max_iters, double max_r, multi_index_t N, multi_real_t length, T* u0, T* u1, const T* b, const char* mask, T* scratch) {
        return _jacobi(max_iters, max_r, 1.0, u0, u1, b, N, length, mask, scratch);
    }

    static unsigned int smooth(unsigned int max_iters, double max_r, multi_index_t N, multi_real_t length, T* u0, T* u1, const T* b, const char* mask, T* scratch) {
        return _jacobi(max_iters, max_r, 4.0/5.0, u0, u1, b, N, length, mask, scratch);
    }

    static double norm_sub(multi_index_t N, const T* vec0, const T* vec1, T* scratch) {
        dim3 block(1, 512, 1);
        dim3 grid((unsigned)ceil((N[0]+2)/(double)block.x), (unsigned)ceil((N[1]+2)/(double)block.y), 1);
        cuda::sub<<<grid, block>>>(scratch, vec0, vec1, size_N(N));
        cuchck_last();
        return norm_cu(scratch, size_N(N));
    }

    static void add_correction(multi_index_t N, T* u, const T* e, const char* mask) {
        dim3 block(1, 512, 1);
        dim3 grid((unsigned)ceil(N[0]/(double)block.x), (unsigned)ceil(N[1]/(double)block.y), 1);
        cuda::addAssign2Dinner<<<grid, block>>>(u, e, N); //update only inner points
        cuchck_last();
    }
};

// explicit instantiations
template struct Fn_CUDA_mem<float>;
template struct Fn_CUDA_mem<double>;
template struct Fn_laplace_cuda<float>;
template struct Fn_laplace_cuda<double>;
