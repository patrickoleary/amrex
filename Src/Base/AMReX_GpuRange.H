#ifndef AMREX_GPU_RANGE_H_
#define AMREX_GPU_RANGE_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuError.H>
#include <AMReX_Box.H>
#include <AMReX_TypeTraits.H>

namespace amrex {

template <typename T, typename = amrex::EnableIf_t<std::is_integral<T>::value> >
bool isEmpty (T n) noexcept { return n <= 0; }

AMREX_FORCE_INLINE bool isEmpty (Box const& b) noexcept { return b.isEmpty(); }

namespace Gpu {

namespace range_detail {

//! integer version
template <typename T, typename = amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_GPU_HOST_DEVICE
long size (T const& b) noexcept { return static_cast<long>(b); }

template <typename T, typename = amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_GPU_HOST_DEVICE
long at (T const& b, long offset) noexcept { return offset; }

//! Box version
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE long size (Box const& b) noexcept
{
#if AMREX_DEVICE_COMPILE
    return b.numPts();
#else
    return 1;
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE Box at (Box const& b, long offset) noexcept
{
#if AMREX_DEVICE_COMPILE
    auto len = b.length3d();
    long k = offset / (len[0]*len[1]);
    long j = (offset - k*(len[0]*len[1])) / len[0];
    long i = (offset - k*(len[0]*len[1])) - j*len[0];
    IntVect iv{AMREX_D_DECL(static_cast<int>(i),
                            static_cast<int>(j),
                            static_cast<int>(k))};
    iv += b.smallEnd();
    return Box(iv,iv,b.type());
#else
    return b;
#endif
}

template <typename T>
struct range_impl
{
    AMREX_GPU_HOST_DEVICE
    explicit range_impl (T const& b) noexcept : m_b(b), m_n(range_detail::size(b)) {}

#ifdef AMREX_USE_DPCPP
    range_impl (T const& b, long gid, long grange) noexcept
        : m_b(b), m_n(range_detail::size(b)), m_gid(gid), m_grange(grange) {}
#endif

    struct iterator
    {
        AMREX_GPU_HOST_DEVICE
        iterator (T const& b, long i, long s) noexcept : mi_b(b), mi_i(i), mi_s(s)  {}

        AMREX_GPU_HOST_DEVICE
        void operator++ () noexcept { mi_i += mi_s; }

        AMREX_GPU_HOST_DEVICE
        bool operator!= (iterator const& rhs) const noexcept { return mi_i < rhs.mi_i; }

        AMREX_GPU_HOST_DEVICE
        T operator* () const noexcept { return range_detail::at(mi_b,mi_i); }

    private:
        T const& mi_b;
        long mi_i;
        long mi_s;
    };

    AMREX_GPU_HOST_DEVICE
    iterator begin () const noexcept {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return iterator(m_b, blockDim.x*blockIdx.x+threadIdx.x, blockDim.x*gridDim.x);
#elif defined (__SYCL_DEVICE_ONLY__)
        return iterator(m_b, m_gid, m_grange);
#else
        return iterator(m_b,0,1);
#endif
    }

    AMREX_GPU_HOST_DEVICE
    iterator end () const noexcept { return iterator(m_b,m_n,0); }

private:
    T m_b;
    long m_n;
#ifdef AMREX_USE_DPCPP
    long m_gid;
    long m_grange;
#endif
};
}

#ifdef AMREX_USE_DPCPP
template <typename T>
range_detail::range_impl<T> Range (T const& b, long gid, long grange) noexcept {
    return range_detail::range_impl<T>(b,gid,grange);
}
#endif

template <typename T>
AMREX_GPU_HOST_DEVICE
range_detail::range_impl<T> Range (T const& b) noexcept { return range_detail::range_impl<T>(b); }

}}

#endif
