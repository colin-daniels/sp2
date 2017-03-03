#ifndef SP2_BLAS_HPP
#define SP2_BLAS_HPP

#include <vector>
#include <cmath>

namespace sp2 {

////////////////////////////////////////////////////////////////////////////////
// Miscellaneous BLAS-like routines                                           //
////////////////////////////////////////////////////////////////////////////////

/// x <- alpha * x
template<class T>
inline void vscal(T alpha, std::vector<T> &x);

/// y <- alpha * x + y
template<class T>
inline void vaxpy(T alpha, const std::vector<T> &x, std::vector<T> &y);

/// x . y
template<class T>
inline T vdot(const std::vector<T> &x, const std::vector<T> &y);

/// x <- x / ||x||
template<class T>
inline void vnormalize(std::vector<T> &x);

/// maximum norm, aka max{abs(x)}
template<class T>
inline T max_norm(const std::vector<T> &x);

/// euclidean distance
template<class T>
inline T vdist(const std::vector<T> &x, const std::vector<T> &y);

////////////////////////////////////////////////////////////////////////////////
// Loop/Template BLAS operations                                              //
////////////////////////////////////////////////////////////////////////////////

// level 1
template<unsigned int N>
static inline void loop_dscal(const double alpha, double *x)
    __attribute__((always_inline));

template<unsigned int N>
static inline void loop_daxpy(const double alpha, const
double *__restrict__ x, double *__restrict__ y)
    __attribute__((always_inline));

template<unsigned int N>
static inline void loop_dcopy(const double *__restrict__ x,
                              double *__restrict__ y)
    __attribute__((always_inline));

template<unsigned int N>
static inline double loop_ddot(const double *__restrict__ x,
                               const double *__restrict__ y)
    __attribute__((always_inline));

// level 2
template<unsigned int M, unsigned int N>
static inline void loop_dger(const double *__restrict__ x,
                             const double *__restrict__ y,
                             double *__restrict__ result)
    __attribute__((always_inline));

template<unsigned int M, unsigned int N>
static inline void loop_dgemv(const double *__restrict__ A,
                              const double *__restrict__ x,
                              double *__restrict__ result)
    __attribute__((always_inline));


// partially specialized
template<unsigned int N>
static inline void loop_dgemv(const double *__restrict__ A,
                              const unsigned int M,
                              const double *__restrict__ x,
                              double *__restrict__ result)
    __attribute__((always_inline));

// level 3

// specialized dgemm for A B' = C
//   C/result(MxN) A(MxK) B(NxK), should only use where B is small
template<unsigned int N, unsigned int K>
static inline void loop_dgemm(const double *__restrict__ A,
                              const unsigned int M,
                              const double *__restrict__ B,
                              double *__restrict__ result)
    __attribute__((always_inline));

////////////////////////////////////////////////////////////////////////////////
// Miscellaneous implementations                                              //
////////////////////////////////////////////////////////////////////////////////

/// x <- alpha * x
template<class T>
inline void vscal(T alpha, std::vector<T> &x)
{
    for (auto &v : x)
        v *= alpha;
}

/// y <- alpha * x + y
template<class T>
inline void vaxpy(T alpha, const std::vector<T> &x, std::vector<T> &y)
{
    for (std::size_t i = 0; i < x.size(); ++i)
        y[i] += alpha * x[i];
}

/// x . y
template<class T>
inline T vdot(const std::vector<T> &x, const std::vector<T> &y)
{
    auto result = T{};
    for (std::size_t i = 0; i < x.size(); ++i)
        result += x[i] * y[i];
    return result;
}

/// x <- x / ||x||
template<class T>
inline void vnormalize(std::vector<T> &x)
{
    vscal(T{1} / static_cast<T>(std::sqrt(vdot(x, x))), x);
}

/// maximum norm, aka max{abs(x)}
template<class T>
inline T max_norm(const std::vector<T> &x)
{
    if (x.empty())
        return T{};

    auto max_elem = std::abs(x[0]);
    for (const auto &elem : x)
        max_elem = std::max(max_elem, std::abs(elem));
    return max_elem;
}

/// euclidean distance
template<class T>
inline T vdist(const std::vector<T> &x, const std::vector<T> &y)
{
    auto sum = T{};
    for (std::size_t i = 0; i < x.size(); ++i)
        sum += (x[i] - y[i]) * (x[i] - y[i]);
    return static_cast<T>(std::sqrt(sum));
}

////////////////////////////////////////////////////////////////////////////////
// Loop/Template BLAS operation implementations                               //
////////////////////////////////////////////////////////////////////////////////

// level 1
template<unsigned int N>
static inline void loop_dscal(const double alpha, double *x)
{
    for (auto i = 0u; i < N; ++i)
        x[i] *= alpha;
}

template<unsigned int N>
static inline void loop_daxpy(const double alpha,
                              const double *__restrict__ x,
                              double *__restrict__ y)
{
    for (auto i = 0u; i < N; ++i)
        y[i] += alpha * x[i];
}

template<unsigned int N>
static inline void
loop_dcopy(const double *__restrict__ x, double *__restrict__ y)
{
    for (auto i = 0u; i < N; ++i)
        y[i] = x[i];
}

template<unsigned int N>
static inline double loop_ddot(const double *__restrict__ x,
                               const double *__restrict__ y)
{
    double result = 0;
    for (auto i = 0u; i < N; ++i)
        result += x[i] * y[i];
    return result;
}

// level 2
template<unsigned int M, unsigned int N>
static inline void loop_dger(const double *__restrict__ x,
                             const double *__restrict__ y,
                             double *__restrict__ result)
{
    for (unsigned int i = 0; i < M; ++i)
        for (unsigned int j = 0; j < N; ++j)
            result[i * N + j] = x[i] * y[j];
}

template<unsigned int M, unsigned int N>
static inline void loop_dgemv(const double *__restrict__ A,
                              const double *__restrict__ x,
                              double *__restrict__ result)
{
    for (unsigned int i = 0; i < M; ++i)
    {
        result[i] = 0;
        for (unsigned int j = 0; j < N; ++j)
            result[i] += A[i * N + j] * x[j];
    }
}

template<unsigned int N>
static inline void
loop_dgemv(const double *__restrict__ A, const unsigned int M,
           const double *__restrict__ x, double *__restrict__ result)
{
    for (unsigned int i = 0; i < M; ++i)
    {
        result[i] = 0;
        for (unsigned int j = 0; j < N; ++j)
            result[i] += A[i * N + j] * x[j];
    }
}

// level 3
// specialized dgemm for A B' = C
//   C/result(MxN) A(MxK) B(NxK)
template<unsigned int N /* rows of B */, unsigned int K /* cols of A and B */>
static inline void
loop_dgemm(const double *__restrict__ A, const unsigned int M, // rows of a
           const double *__restrict__ B, double *__restrict__ result)
{
    for (unsigned int i = 0; i < M; ++i)
    {
        for (unsigned int j = 0; j < N; ++j)
        {
            result[i * N + j] = 0;
            for (unsigned int k = 0; k < K; ++k)
                result[i * N + j] += A[i * K + k] * B[j * K + k];
        }
    }
}

} // namespace sp2

#endif // SP2_BLAS_HPP
