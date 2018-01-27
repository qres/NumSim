#include "grid_size.hpp"

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
