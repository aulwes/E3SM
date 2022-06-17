#ifndef OCN_YAKL
#define OCN_YAKL

#include <iostream>
#ifdef YAKL_ARCH_CUDA
#include <cuda_runtime.h>
#elif defined(YAKL_ARCH_HIP)
#include <hip/hip_runtime.h>
#endif
#include "mpas_ocn_yakl_types.hxx"

#define SIGN(a,b) ( (b) >= 0 ) ? std::abs(a) : -std::abs(a)

#define YAKL_LOCAL(lvar, var) auto & lvar = *(var)
#define YAKL_LOCAL_NS(ns, var) auto & var = *(ns::var)

typedef struct
{
    int shape[10];
    void * ptr;
} ocn_yakl_type;

template <typename T>
std::ostream & operator << (std::ostream & os, const std::vector<T> & vec)
{
    for(auto elem : vec)
    {
        os<<elem<< " ";
    }
    return os;
}

template<typename R, int Mem, int Dims>
std::vector<int> getBounds(yakl::Array<R,Dims,Mem,yakl::styleFortran> & arr)
{
    std::vector<int>       bnds;

    auto rnk = arr.get_rank();
    for ( int i = 0; i < rnk; ++i )
    {
        bnds.push_back(arr.extent(i));
    }
    return bnds;
}


template <typename R, typename...T>
yakl::Array<R,sizeof...(T),yakl::memDefault,yakl::styleFortran> *
yakl_create_array(const char * name, T...dims)
{
    typedef yakl::Array<R,sizeof...(T),yakl::memDefault,yakl::styleFortran> ret_type;

    ret_type * w_var_p = new ret_type(name, dims...);

    return w_var_p;
}

template <typename R, typename...T>
yakl::Array<R,sizeof...(T),yakl::memDefault,yakl::styleFortran> *
yakl_wrap_array(const char * name, R * var_p, T...dims)
{
    typedef yakl::Array<R,sizeof...(T),yakl::memHost,yakl::styleFortran> arr_type;
    typedef yakl::Array<R,sizeof...(T),yakl::memDefault,yakl::styleFortran> ret_type;

    arr_type h_var("h_var_p", var_p, dims...);
    ret_type * w_var_p = new ret_type(name, dims...);

    h_var.deep_copy_to(*w_var_p);

    return w_var_p;
}


template <typename R, int N>
void
yakl_update_host(yakl::Array<R,N,yakl::memDefault,yakl::styleFortran> * d_var,
                 R * h_var_p,yakl::yakl_stream_t stream = 0 )
{
    typedef yakl::Array<R,N,yakl::memHost,yakl::styleFortran>  host_type;

    auto & rw_var = *d_var;
    host_type   h_var("h_var_p", h_var_p, getBounds(rw_var));
    rw_var.deep_copy_to(h_var, stream);
}

/*
template <typename R, int N>
void
yakl_update_host_from_host(yakl::Array<R,N,yakl::memHost,yakl::styleFortran> & d_var,
                 R * h_var_p,yakl::yakl_stream_t stream = 0 )
{
    typedef yakl::Array<R,N,yakl::memHost,yakl::styleFortran>  host_type;
    host_type   h_var("h_var_p", h_var_p, getBounds(d_var));
    yakl::memcpy_host_to_host( h_var.myData , d_var.myData , d_var.totElems(), 0);
}
*/

template <typename R, int N>
void
yakl_update_device(yakl::Array<R,N,yakl::memDefault,yakl::styleFortran> * d_var,
                 R * h_var_p,yakl::yakl_stream_t stream = 0 )
{
    typedef yakl::Array<R,N,yakl::memHost,yakl::styleFortran>  host_type;

    auto & rw_var = *d_var;
    host_type   h_var("h_var_p", h_var_p, getBounds(rw_var));
    h_var.deep_copy_to(rw_var, stream);
}

template <typename...T>
yakl::Array<real,sizeof...(T),yakl::memDefault,yakl::styleFortran> *
yakl_create_real(const char * name, T...Dims)
{
    return yakl_create_array<real>(name, Dims...);
}

template <typename...T>
yakl::Array<int,sizeof...(T),yakl::memDefault,yakl::styleFortran> *
yakl_create_int(const char * name, T...Dims)
{
    return yakl_create_array<int>(name, Dims...);
}

void yakl_stream_create(yakl::yakl_stream_t * stream);
void yakl_stream_destroy(yakl::yakl_stream_t stream);

void yakl_event_create(yakl::yakl_event_t * event);
void yakl_event_destroy(yakl::yakl_event_t event);

void yakl_record_event(yakl::yakl_event_t, yakl::yakl_stream_t);
void yakl_stream_wait(yakl::yakl_stream_t, yakl::yakl_event_t);

#endif
