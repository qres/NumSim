
#include <cstring>
#include <algorithm>
#include <math.h>
#include <limits>

// access matrix indices
#define row_major (1)
#if row_major
#define get_matrix(a,i,j,N) ((a)[(i) * ((N)+2) + (j)])
#else
#define get_matrix(a,i,j,N) ((a)[(i) + (j) * ((N)+2)])
#endif
#define get_matrix_NxN(a,i,j) (get_matrix(a,i,j,N))

using std::max;
using std::abs;

bool verbose = false;
bool verbose_levels = false;

struct Grid2D {
    static unsigned int dim() {
        return 2;
    }

    static unsigned int size_N(unsigned int N) {
        return (N+2)*(N+2);
    }

    static unsigned int size_n(unsigned int N) {
        unsigned int n = N/2;
        return (n+2)*(n+2);
    }
};

template<typename T>
struct Fn_CPU_mem {
    static bool mem_device_host_equal() {
        return true;
    }

    /// allocates `new T[count]`
    static T* malloc_typed(unsigned int count) {
        return new T[count];
    }

    static void free_typed(T* ptr) {
        delete[] ptr;
    }

    static void memcpy_typed_HostToDevice(T* dst, const T* src, unsigned int count) {
        memcpy(dst, src, sizeof(T) * count);
    }

    static void memcpy_typed_DeviceToHost(T* dst, const T* src, unsigned int count) {
        memcpy(dst, src, sizeof(T) * count);
    }

    static void memcpy_typed_DeviceToDevice(T* dst, const T* src, unsigned int count) {
        memcpy(dst, src, sizeof(T) * count);
    }

    /// sets `dst[0..count] = (T)0`
    static void zero(double* dst, const unsigned int count) {
        std::fill_n(dst, count, (T)0 );
    }
};

template<typename T>
void print_vec(unsigned int count, const T* v, bool trailing_comma, std::ostream& out = std::cout) {
    out << "[";
    if (count > 0) out << v[0] << ",";
    for (unsigned int i(1); i<count-1; ++i) {
        out << v[i] << ", ";
    }
    if (count-1 > 0) out << v[count-1];
    out << "]";
    if (trailing_comma) out << ",";
    out << std::endl;
}

template<typename T>
T norm2_residuum_laplace(unsigned int N, const T* u, const T* b) {
    #define m get_matrix_NxN
    T sum = 0;
    // only inner values as boundary is given
    for (unsigned int i(1); i<=N; ++i) {
        for (unsigned int j(1); j<=N; ++j) {
            T entry = m(b,i,j) + (N+1)*(N+1)*(4*m(u,i,j) - 1*m(u,i-1,j) - 1*m(u,i+1,j) - 1*m(u,i,j-1) - 1*m(u,i,j+1));
            sum += entry*entry;
        }
    }
    return sqrt(sum);
    #undef m
}

template<typename T>
T norm2w_residuum_laplace(unsigned int N, const T* u, const T* b) {
    #define m get_matrix_NxN
    T sum = 0;
    // only inner values as boundary is given
    for (unsigned int i(1); i<=N; ++i) {
        for (unsigned int j(1); j<=N; ++j) {
            T entry = m(b,i,j) + (N+1)*(N+1)*(4*m(u,i,j) - 1*m(u,i-1,j) - 1*m(u,i+1,j) - 1*m(u,i,j-1) - 1*m(u,i,j+1));
            sum += entry*entry;
        }
    }
    return sqrt(sum/N);
    #undef m
}

template<typename T>
T normoo_residuum_laplace(unsigned int N, const T* u, const T* b) {
    #define m get_matrix_NxN
    T norm = 0;
    // only inner values as boundary is given
    for (unsigned int i(1); i<=N; ++i) {
        for (unsigned int j(1); j<=N; ++j) {
            T entry = m(b,i,j) + (N+1)*(N+1)*(4*m(u,i,j) - 1*m(u,i-1,j) - 1*m(u,i+1,j) - 1*m(u,i,j-1) - 1*m(u,i,j+1));
            norm = max(norm,abs(entry));
        }
    }
    return norm;
    #undef m
}

template<typename T>
T norm_residuum_laplace(unsigned int N, const T* u, const T* b) {
    return normoo_residuum_laplace(N,u,b);
}

template<typename T>
T norm2_sub(unsigned int N, const T* u_h, const T* u) {
    T sum = 0;
    // only inner values as boundary is given
    for (unsigned int i(0); i<=N+1; ++i) {
        T entry = u_h[i] - u[i];
        sum += entry*entry;
    }
    return sqrt(sum);
}

template<typename T>
T norm2w_sub(const T* u_h, const T* u, unsigned int count) {
    T sum = 0;
    // only inner values as boundary is given
    for (unsigned int i(0); i < count; ++i) {
        T entry = u_h[i] - u[i];
        sum += entry*entry;
    }
    return sqrt(sum/count);
}

template<typename T>
T normoo_sub(const T* u_h, const T* u, unsigned int count) {
    T norm = 0;
    // only inner values as boundary is given
    for (unsigned int i(0); i < count; ++i) {
        T entry = u_h[i] - u[i];
        norm = max(norm,abs(entry));
    }
    return norm;
}

template<typename T>
T normoo(const T* u_h, unsigned int count) {
    T norm = 0;
    // only inner values as boundary is given
    for (unsigned int i(0); i < count; ++i) {
        T entry = u_h[i];
        norm = max(norm,abs(entry));
    }
    return norm;
}

template<typename T>
T norm(const T* u_h, unsigned int count) {
    return normoo(u_h, count);
}

template<typename T>
T norm_sub(const T* u_h, const T* u, unsigned int count) {
    return normoo_sub(u_h, u, count);
}

template<typename T>
void residuum_laplace(unsigned int N, const T* u, const T* b, T* res) {
    #define m get_matrix_NxN
    for (unsigned int i(0); i<=N+1; ++i) {
        m(res,i,0) = 0;//m(b,i,0) - m(u,i,0);
        m(res,0,i) = 0;//m(b,0,i) - m(u,0,i);
        m(res,i,N+1) = 0;//m(b,i,N+1) - m(u,i,N+1);
        m(res,N+1,i) = 0;//m(b,N+1,i) - m(u,N+1,i);
    }
    for (unsigned int i(1); i<=N; ++i) {
        for (unsigned int j(1); j<=N; ++j) {
            m(res,i,j) = m(b,i,j) + (N+1)*(N+1)*(4*m(u,i,j) - 1*m(u,i-1,j) - 1*m(u,i+1,j) - 1*m(u,i,j-1) - 1*m(u,i,j+1));
        }
    }
    #undef m
}

//#define JACOBI_DONT_CHECK_RESIDUUM
#ifdef JACOBI_DONT_CHECK_RESIDUUM
    #warning !!! !!! !!! JACOBI DOES NOT CHECK RESIDUUMS !!! !!! !!!
#endif

/*
 * jacobi interation for Ax = b with the Matrix corresponding to -d_xx-d_yy
 * -1 * [[-4  1  0       ... 1     ...0]
 *       [ 1 -4  1  0    ...    1  ...0]
 *       [ 0  1 -4  1  0 ...          0]
 *       [...            ...        ...]
 *       [ 1  ...        ...     ...  1]
 *       [...            ...        ...]
 *       [ 0             ... 1 -4  1  0]
 *       [ 0...     1    ... 0  1 -4  1]
 *       [ 0...        1 ...    0  1 -4]] / h^2
 *
 * v0: initial guess
 * v1: 2nd buffer
 * N: number of inner values [boundary, val, val, boundary]=>N=2
 * b: rhs
 * Returns: #iters, v0 as approximation, v1 with the previous approximation
 */
template<typename T>
unsigned int jacobi_laplace(unsigned int max_iters, double max_r, double omega, unsigned int N, T* u0, T* u1, const T* b) {
    T* const dst = u0;
    unsigned int iters = 1;
    T r = norm_residuum_laplace(N, u0, b);
    while (iters <= max_iters && r >= max_r) {
        T r_old = r;

        // jacobi update
        // x' = D_inv (b - (A - D) x)
        // here: D_inv = 1/(4*N*N), (A-D)x = (-v_-1 + -v_+1)*N*N
        // => x' =  1 / (4*N*N) (b - (-v_-1 + -v_+1)*N*N)
        // => x' = 1/4 (b/(N*N) + v_-1 + v_+1)

        #define m get_matrix_NxN

        // boundary is known
        for (unsigned int i(0); i<=N+1; ++i) {
            m(u1,i,0) = m(u0,i,0);
            m(u1,0,i) = m(u0,0,i);
            m(u1,i,N+1) = m(u0,i,N+1);
            m(u1,N+1,i) = m(u0,N+1,i);
        }
        for (unsigned int i(1); i<=N; ++i) {
            for (unsigned int j(1); j<=N; ++j) {
                m(u1,i,j) = 0.25 * (-m(b,i,j) / ((N+1)*(N+1)) + m(u0,i-1,j) + m(u0,i+1,j) + m(u0,i,j-1) + m(u0,i,j+1));
                m(u1,i,j) = (1.0-omega) * m(u0,i,j) + omega * m(u1,i,j);
            }
        }

        #undef m

        // swap input and output
        T* swap = u1;
        u1 = u0;
        u0 = swap;

#ifndef JACOBI_DONT_CHECK_RESIDUUM
        r = norm_residuum_laplace(N, u0, b);
        if (r == r_old) {
          if (verbose) std::cerr << "breaking jacobi interation: no change. N=" << N << " iters=" << iters <<  " residuum=" << r << std::endl;
          if (max_r != 0) break;
        }
#endif

        ++iters;
    }
    if (verbose && iters == max_iters) {
            std::cerr << "breaking jacobi interation: too many iterations. N=" << N << " iters=" << iters <<  " residuum=" << r << std::endl;
    }
    if (dst != u0) {
        memcpy(dst, u0, (N+2)*(N+2) * sizeof(T));
    }
    return iters - 1;
}

/*
 * restriction unsing the 1/16*[1 2 1
 *                              2 4 2
 *                              1 2 1] stencil
 * N: number of inner values of the fine grid [boundary, val, val boundary] => N = 2
 */
template<typename T>
unsigned int restrict_fw_2D(unsigned int N, const T* v_N, T* v_n) {
    /*     01     ....    NN+1
       |   11111111111111111   17 15
       Y   1 1 1 1 1 1 1 1 1   9  7
       Y   1   1   1   1   1   5  3
       Y   1       1       1   3  1
       V   1               1   2  0
    */
    #define _v_n(i,j) (get_matrix(v_n,i,j,n))
    #define _v_N(i,j) (get_matrix(v_N,i,j,N))
    const unsigned int n = N/2;
    // keep boundary condition
    for (unsigned int i(0); i<=N+1; i+=2) {
        _v_n(i/2, 0) = _v_N(i,0);
        _v_n(0, i/2) = _v_N(0,i);
        _v_n(i/2, n+1) = _v_N(i,N+1);
        _v_n(n+1, i/2) = _v_N(N+1,i);
    }
    _v_n(n+1, n+1) = _v_N(N+1,N+1);

    // stencil
    for (unsigned int i(2); i<=N; i+=2) {
        for (unsigned int j(2); j<=N; j+=2) {
            _v_n(i/2,j/2) = 1.0/16.0 * (
                1*_v_N(i-1, j-1) + 2*_v_N(i, j-1) + 1*_v_N(i+1, j-1) +
                2*_v_N(i-1, j  ) + 4*_v_N(i, j  ) + 2*_v_N(i+1, j  ) +
                1*_v_N(i-1, j+1) + 2*_v_N(i, j+1) + 1*_v_N(i+1, j+1)
            );
        }
    }
    return n;
    #undef _v_n
    #undef _v_N
}

/*
 * linear interpolation
 * n : number of inner values of the coarse grid [boundary, val, val boundary] => N = 2
 */
template<typename T>
unsigned int interplolate_2D(unsigned int n, const T* v_n, T* v_N) {
    const unsigned int N = n*2 + 1;
    #define _v_n(i,j) (get_matrix(v_n,i,j,n))
    #define _v_N(i,j) (get_matrix(v_N,i,j,N))
    // keep boundary condition
    // interpolate
    /*
    x 0 x 0 x 0 x
    0 0 0 0 0 0 0
    x 0 x 0 x 0 x
    0 0 0 0 0 0 0
    x 0 x 0 x 0 x
    */
    for (unsigned int i(1); i<=n+1; ++i) {
        for (unsigned int j(1); j<=n+1; ++j) {
            if (2*i != N+1 && 2*j != N+1) _v_N(2*i,     2*j    ) = _v_n(i,j);
            if (2*j != N+1)               _v_N(2*i - 1, 2*j    ) = 0.5 * (_v_n(i-1,j) + _v_n(i,j));
            if (2*i != N+1)               _v_N(2*i,     2*j - 1) = 0.5 * (_v_n(i,j-1) + _v_n(i,j));
            _v_N(2*i - 1, 2*j - 1) = 0.25 * (_v_n(i-1,j) + _v_n(i,j-1) + _v_n(i-1,j-1) + _v_n(i,j));
        }
    }

    return N;
    #undef _v_n
    #undef _v_N
}

template<typename T>
struct Fn_laplace : Fn_CPU_mem<T>, Grid2D {

    static void restrict(unsigned int N, const T* v_N, T* v_n) {
        restrict_fw_2D(N, v_N, v_n);
    }

    static void interpolate(unsigned int n, const T* v_n, T* v_N) {
        interplolate_2D(n, v_n, v_N);
    }

    static T norm_residuum(unsigned int N, const T* u, const T* b, T* scratch = 0) {
        return  norm_residuum_laplace(N, u, b);
    }

    static double norm(unsigned int N, const T* vec0) {
        return ::norm(vec0, size_N(N));
    }

    static double norm_sub(unsigned int N, const T* vec0, const T* vec1, T* scratch = 0) {
        return ::norm_sub(vec0, vec1, size_N(N));
    }

    static void residuum(unsigned int N, const T* u, const T* b, T* res) {
        residuum_laplace(N, u, b, res);
    }

    static unsigned int solve(unsigned int max_iters, double max_r, unsigned int N, T* u0, T* u1, const T* b, T* scratch = 0) {
        return jacobi_laplace(max_iters, max_r, 1.0, N, u0, u1, b);
    }

    static unsigned int smooth(unsigned int max_iters, double max_r, unsigned int N, T* u0, T* u1, const T* b, T* scratch = 0) {
        return jacobi_laplace(max_iters, max_r, 4.0/5.0, N, u0, u1, b);
    }

    static void add_correction(unsigned int N, T* u, const T* e) {
        #define m get_matrix_NxN
        for (unsigned int i(1); i <= N; ++i) {
            for (unsigned int j(1); j <= N; ++j) {
                m(u,i,j) += m(e,i,j);
            }
        }
        #undef m
    }
};

template<typename F, typename T>
unsigned int solve_mg_flat(const Cfg& cfg, unsigned int N_coarse, T* _u0, const T* _b) {
    T** u0  = new T*[cfg.num_levels];
    T* scratch  = F::malloc_typed(F::size_N(N_coarse)); //can be used as ping pong buffer by all levels
    T** b = new T*[cfg.num_levels];
    T** r_0 = new T*[cfg.num_levels];
    T** e_0 = new T*[cfg.num_levels];
    T** res0 = new T*[cfg.num_levels];
    T** res1 = new T*[cfg.num_levels];
    unsigned int* N = new unsigned int[cfg.num_levels];

    // preallocate memory
#define single_buffer
#ifdef single_buffer
    // use one large allocation for all memory needed
    unsigned int total_size = 0;
    for (unsigned int i = 0, N_i = N_coarse; i < cfg.num_levels; ++i, N_i = N_i/2) {
        total_size += F::size_N(N_i);
    }
    // 6 buffer hirachies (minus _u0 and _b if they can be used directly)
    T* buffer = F::malloc_typed(total_size * 6 - (F::mem_device_host_equal()?2*F::size_N(N_coarse):0));

    T* ptr = buffer;
    for (unsigned int i = 0, N_i = N_coarse; i < cfg.num_levels; ++i, N_i = N_i/2) {
        if (i==0 && F::mem_device_host_equal()) {
            u0[i] = _u0;
            b[i] = const_cast<T*>(_b); // we will never overwrite this value
        } else {
            u0[i] = ptr;
            ptr += F::size_N(N_i);
            b[i] = ptr;
            ptr += F::size_N(N_i);
        }

        r_0[i] = ptr;
        ptr += F::size_N(N_i);
        e_0[i] = ptr;
        ptr += F::size_N(N_i);
        res0[i] = ptr;
        ptr += F::size_N(N_i);
        res1[i] = ptr;
        ptr += F::size_N(N_i);

        N[i] = N_i;
    }
    if (!F::mem_device_host_equal()) {
        F::memcpy_typed_HostToDevice(u0[0], _u0, F::size_N(N_coarse));
        F::memcpy_typed_HostToDevice(b[0],  _b,  F::size_N(N_coarse));
    }
#else
    // use different allcations
    for (unsigned int i = 0, N_i = N_coarse; i < cfg.num_levels; ++i, N_i = N_i/2) {
        u0[i]  = (i==0 && F::mem_device_host_equal() ? _u0 : F::malloc_typed(F::size_N(N_i)));
        b[i] = (i==0 && F::mem_device_host_equal() ? _b : F::malloc_typed(F::size_N(N_i)));
        r_0[i] = F::malloc_typed(F::size_N(N_i));
        e_0[i] = F::malloc_typed(F::size_N(N_i));
        res0[i] = F::malloc_typed(F::size_N(N_i));
        res1[i] = F::malloc_typed(F::size_N(N_i));
        N[i] = N_i;
    }
#endif

    T* r = new T[cfg.num_levels];
    T* r_old = new T[cfg.num_levels];

    unsigned int* iters = new unsigned int[cfg.num_levels];
    memset(iters, 0, cfg.num_levels * sizeof(unsigned int));

    unsigned int current_level = 0;

    r_old[current_level] = std::numeric_limits<T>::infinity();
    F::residuum(N[current_level], u0[current_level], b[current_level], res1[current_level]);
    r[current_level] = F::norm(N[current_level], res1[current_level]);

    const int SWEEP_UP = 1;
    const int SWEEP_DOWN = 2;
    const int SOLVE = 4;
    int next = SWEEP_DOWN;

    while(true) {

        if (current_level == cfg.num_levels - 1) {
            if (iters[current_level] == 0) {
                next = SOLVE;
            } else {
                if (current_level == 0) {
                    break;
                } else {
                    next = SWEEP_UP;
                }
            }
        } else {
            if (iters[current_level] < cfg.levels[current_level].max_iters && r[current_level] >= cfg.levels[current_level].max_res) {
                next = SWEEP_DOWN;
            } else {
                if (current_level == 0) {
                    break;
                } else {
                    next = SWEEP_UP;
                }
            }
        }

        if (next == SWEEP_DOWN) {
            // down-sweep
            iters[current_level]++;
            if(verbose) {for (unsigned int i = 0; i < current_level; ++i) std::cout << " "; std::cout << std::scientific << r[current_level] << std::endl;}
            if(verbose_levels) {for (unsigned int i = 0; i < current_level; ++i) std::cout << " "; std::cout << "." << std::endl;}
            F::smooth(
                cfg.levels[current_level].max_iters_jacobi_pre,
                cfg.levels[current_level].max_res_jacobi_pre,
                N[current_level],
                u0[current_level], scratch, b[current_level],
                r_0[current_level]); //scratch: will be overwritten soon

            //residuum = b-Lu
            F::residuum(N[current_level], u0[current_level], b[current_level], r_0[current_level]);

            // solve Le = r
            F::restrict(N[current_level], r_0[current_level], b[current_level + 1]);
            F::zero(u0[current_level+1], F::size_n(N[current_level]));
            if(verbose_levels) {for (unsigned int i = 0; i < current_level; ++i) std::cout << " "; std::cout << "\\" << std::endl;}

            //init new level
            iters[current_level+1] = 0;
            F::residuum(N[current_level+1], u0[current_level+1], b[current_level+1], res1[current_level+1]);
            r[current_level+1] = F::norm(N[current_level+1], res1[current_level+1]);
            r_old[current_level+1] = std::numeric_limits<T>::infinity();

            current_level++;
            continue;

        } else if (next == SWEEP_UP) {
            // up-sweep
            current_level--;
            F::interpolate(N[current_level+1], u0[current_level+1], e_0[current_level]);
            F::add_correction(N[current_level], u0[current_level], e_0[current_level]);
            if(verbose_levels) {for (unsigned int i = 0; i < current_level; ++i) std::cout << " "; std::cout << "/" << std::endl;}

            if(verbose_levels) {for (unsigned int i = 0; i < current_level; ++i) std::cout << " "; std::cout << "." << std::endl;}
            F::smooth(
                cfg.levels[current_level].max_iters_jacobi_post,
                cfg.levels[current_level].max_res_jacobi_post,
                N[current_level],
                u0[current_level], scratch, b[current_level],
                res1[current_level]); //scratch: will be overwritten soon

            std::swap(res0[current_level], res1[current_level]);
            F::residuum(N[current_level], u0[current_level], b[current_level], res1[current_level]);

            r_old[current_level] = r[current_level];
            r[current_level] = F::norm(N[current_level], res1[current_level]);
            if(verbose) {for (unsigned int i = 0; i < current_level; ++i) std::cout << " "; std::cout << std::scientific << r[current_level] << std::endl;}

            T res_diff = F::norm_sub(N[current_level], res0[current_level], res1[current_level], scratch);
            if (res_diff < cfg.levels[0].max_res/pow(10,current_level)  || res_diff/ r_old[current_level] < cfg.levels[0].max_res/pow(10,current_level)) {
                if (verbose) {
                    for (unsigned int i = 0; i < current_level; ++i) std::cout << " ";
                    std::cout << (res_diff < 1e-12 ? "" : "res_diff | ") << (res_diff/r_old[current_level] < 1e-12 ? "" : "res_diff_rel") << std::endl;
                }
                r[current_level] = 0; //break current level interation
            }

            continue;
        } else {
            //solve
            iters[current_level]++;
            if(verbose_levels) {for (unsigned int i = 0; i < current_level; ++i) std::cout << " "; std::cout << "." << std::endl;}
            unsigned int it = F::solve(
                cfg.levels[current_level].max_iters,
                cfg.levels[current_level].max_res,
                N[current_level],
                u0[current_level], scratch, b[current_level],
                e_0[current_level]); //scratch: will be overwritten soon
            if (cfg.count_iters_like_single_grid) {
                iters[current_level] = it;
            }
        }
    }

    F::free_typed(scratch);
    delete[] N;

#ifdef single_buffer
    if (!F::mem_device_host_equal()) {
        F::memcpy_typed_DeviceToHost(_u0, u0[0], F::size_N(N_coarse));
    }
    F::free_typed(buffer);
#else
    for (unsigned int i = 0, N_i = N_coarse; i < cfg.num_levels; ++i, N_i = N_i/2) {
        if (i!=0 || !F::mem_device_host_equal()) {
            F::free_typed(u0[i]);
            F::free_typed(b[i]);
        }
        F::free_typed(r_0[i]);
        F::free_typed(e_0[i]);
        F::free_typed(res0[i]);
        F::free_typed(res1[i]);
    }
#endif

    delete[] u0 ;
    delete[] b;
    delete[] r_0;
    delete[] e_0;
    delete[] res0;
    delete[] res1;

    delete[] r;
    delete[] r_old;

    int iters_0 = iters[0];
    delete[] iters;

    return iters_0;
}