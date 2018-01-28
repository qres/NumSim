#ifndef MG_IMP_CUH
#define MG_IMP_CUH

#include "grid_size.hpp"
#include "cuda.h"
#include "cuda_runtime.h"

void print_cuda_stats() {
    cudaDeviceProp prop;
    cudaError_t status = cudaGetDeviceProperties (&prop, 0); //device 0

    std::cout << "Device name           : " << prop.name << std::endl;
    std::cout << "Compute capability    : " << prop.major << "." << prop.minor << std::endl;
    std::cout << "Warp size             : " << prop.warpSize << std::endl;
    std::cout << "Max threads per core  : " << prop.maxThreadsPerMultiProcessor << std::endl;
    std::cout << "Max threads per block : " << prop.maxThreadsPerBlock << std::endl;
    std::cout << "Max block dimensions  : " << prop.maxThreadsDim[0] << "x" << prop.maxThreadsDim[1] << "x"<< prop.maxThreadsDim[2] << std::endl;
    std::cout << "Max grid dimensions   : " << prop.maxGridSize[0] << "x" << prop.maxGridSize[1] << "x"<< prop.maxGridSize[2] << std::endl;
    std::cout << "Overlap copy & execute: " << prop.deviceOverlap << std::endl;
    std::cout << std::endl;
}

void init_cuda() {
    auto result = cuInit(0);
    if (result != cudaSuccess) {
        std::cout << "FAILED to init cuda (" << result << ")" << std::endl;
        exit(1);
    }
    print_cuda_stats();
}


template<typename T>
struct Fn_CUDA_mem {
    static bool mem_device_host_equal();
    static T* malloc_typed(unsigned int count);
    static void free_typed(T* ptr);
    static void memcpy_typed_HostToDevice(T* dst, const T* src, unsigned int count);
    static void memcpy_typed_DeviceToHost(T* dst, const T* src, unsigned int count);
    static void memcpy_typed_DeviceToDevice(T* dst, const T* src, unsigned int count);
    static void zero(T* dst, const unsigned int count);
};

template<typename T>
struct Fn_laplace_cuda : Fn_CUDA_mem<T>, Grid2D {
    static void restrict(multi_index_t N, const T* v_N, T* v_n);
    static void interpolate(multi_index_t n, const T* v_n, T* v_N, const char* mask);
    static void residuum(multi_index_t N, multi_real_t length, const T* u, const T* b, const char* mask, T* res);
    static double norm(multi_index_t N, const T* vec0);
    static T norm_residuum(multi_index_t N, multi_real_t length, const T* u, const T* b, const char* mask, T* scratch);

    // TODO scratch buffer might be obsolate
    static unsigned int _jacobi(unsigned int max_iters, double max_r, double omega, T* u0, T* u1, const T* b, multi_index_t N, multi_real_t length, const char* mask, T* scratch);
    static unsigned int solve(unsigned int max_iters, double max_r, multi_index_t N, multi_real_t length, T* u0, T* u1, const T* b, const char* mask, T* scratch);
    static unsigned int smooth(unsigned int max_iters, double max_r, multi_index_t N, multi_real_t length, T* u0, T* u1, const T* b, const char* mask, T* scratch);

    static double norm_sub(multi_index_t N, const T* vec0, const T* vec1, T* scratch);
    static void add_correction(multi_index_t N, T* u, const T* e, const char* mask);
};

#endif
