#include <iostream>
#ifdef MPAS_GPTL_TIMERS
#include "gptl.h"
#endif
#include "mpas_ocn_yakl_types.hxx"
#include "mpas_ocn_yakl_c.hxx"

#include <chrono>
#include <ctime>
#include <cmath>
#include <algorithm>    // std::count_if
#include<numeric>

// Timer from https://gist.github.com/mcleary/b0bf4fa88830ff7c882d
class Timer
{
public:
    void start()
    {
        m_StartTime = std::chrono::system_clock::now();
        m_bRunning = true;
    }
    
    void stop()
    {
        m_EndTime = std::chrono::system_clock::now();
        m_bRunning = false;
    }
    
    double elapsedMilliseconds()
    {
        std::chrono::time_point<std::chrono::system_clock> endTime;
        
        if(m_bRunning)
        {
            endTime = std::chrono::system_clock::now();
        }
        else
        {
            endTime = m_EndTime;
        }
        
        return std::chrono::duration_cast<std::chrono::milliseconds>(endTime - m_StartTime).count();
    }
    
    double elapsedSeconds()
    {
        return elapsedMilliseconds() / 1000.0;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> m_StartTime;
    std::chrono::time_point<std::chrono::system_clock> m_EndTime;
    bool                                               m_bRunning = false;
};

namespace {
std::vector<double>     uptimes;

Timer   timer;

}; // namespace

extern "C"
void init_yakl()
{
#ifdef MPAS_GPTL_TIMERS
    //GPTLinitialize();
#endif
    yakl::init();
}

extern "C"
void finalize_yakl()
{
/*
    int cnt = std::count_if(uptimes.begin(), uptimes.end(), [](double v){return v > 1;});
    auto res = std::find_if(uptimes.begin(), uptimes.end(), [](double v){return v > 1;});
    auto total = std::accumulate(uptimes.begin(), uptimes.end(), 0.0);
    double avg = total / uptimes.size();
    const auto [mint,maxt] = std::minmax_element(uptimes.begin(), uptimes.end());
    
    std::cerr << "   Total,min,max,avg update device time = " << total << " " << *mint << " " << *maxt << " " << avg << " for " << uptimes.size() << " updates." << std::endl;
    if ( cnt > 0 )
    {
        std::cerr << " Update cnt exceeding 1 sec = " << cnt << std::endl;
        std::cerr << "     time = " << *res << std::endl;
    }

    cnt = std::count_if(uptimes.begin(), uptimes.end(), [](double v){return v > .5;});
    if ( cnt > 0 ) std::cerr << " Update cnt exceeding .5 sec = " << cnt << std::endl;
*/
    yakl::finalize();
}

extern "C"
void yakl_fence() {yakl::fence();}

extern "C"
void update_host_d2d(d_double_2d_t * w_var1, real * var_p)
{
    d_double_2d_t & rw_var = *w_var1;
    h_double_2d_t   h_var("h_var_p", var_p, getBounds(rw_var));
    w_var1->deep_copy_to(h_var);
}

extern "C"
void update_host_d1d(d_double_1d_t * w_var1, real * var_p)
{
    d_double_1d_t & rw_var = *w_var1;
    h_double_1d_t   h_var("h_var_p", var_p, getBounds(rw_var));
    w_var1->deep_copy_to(h_var);
}

extern "C"
void update_device_d1d(d_double_1d_t * w_var1, real * var_p)
{
    static int cnt = 0;
    
    //if ( cnt > 0 ) return;
    //++cnt;

    h_double_1d_t   h_var("h_var_p", var_p, w_var1->extent(0));
    timer.start();
    yakl::timer_start("update_device_d1d");
    h_var.deep_copy_to(*w_var1);
    yakl::timer_stop("update_device_d1d");
    timer.stop();
    uptimes.push_back(timer.elapsedSeconds());
}

extern "C"
void update_device_d2d(d_double_2d_t * w_var1, real * var_p)
{
    static int cnt = 0;
    
    //if ( cnt > 0 ) return;
    //++cnt;

    h_double_2d_t   h_var("h_var_p", var_p, w_var1->extent(0),w_var1->extent(1));
    //std::cerr << "   update_device_d2d: extent = " << w_var1->extent(0) << " " << w_var1->extent(1) << std::endl;
timer.start();
    yakl::timer_start("update_device_d2d");
    h_var.deep_copy_to(*w_var1);
    yakl::timer_stop("update_device_d2d");
    timer.stop();
    uptimes.push_back(timer.elapsedSeconds());
}

extern "C"
void unwrap_1d(d_double_1d_t * wv) {delete wv;}

extern "C"
void unwrap_2d(d_double_2d_t * wv) {delete wv;}

extern "C"
d_double_1d_t * yakl_create_d1d(int n) {return yakl_create_array<real>("var_d1d", n);}

extern "C"
d_double_2d_t * yakl_create_d2d(int n, int m) {return yakl_create_array<real>("var_d2d", n,m);}

extern "C"
void yakl_delete(d_double_1d_t * w_var1) {delete w_var1;}

extern "C"
void yakl_delete_d2d(d_double_2d_t * w_var1) {delete w_var1;}

extern "C"
d_double_1d_t * wrap_array_1d(real * var_p, int n)
{
    d_double_1d_t * w_var_p = yakl_wrap_array("var_p", var_p, n);

    return w_var_p;
}

extern "C"
d_double_2d_t * wrap_array_2d(real * var_p, int n, int m)
{
    d_double_2d_t * w_var_p = yakl_wrap_array("var_p", var_p, n, m);

    return w_var_p;
}

extern "C"
d_double_3d_t * wrap_array_3d(real * var_p, int n, int m, int k)
{
    auto w_var_p = yakl_wrap_array("var_p", var_p, n, m, k);

    return w_var_p;
}

void yakl_stream_create(yakl::yakl_stream_t * stream)
  {
    #ifdef YAKL_ARCH_CUDA
    cudaStreamCreate(stream);
    #elif YAKL_ARCH_HIP
    hipStreamCreate(stream);
    #else
    stream = nullptr;
    #endif
  }

void yakl_stream_destroy(yakl::yakl_stream_t stream)
  {
    #ifdef YAKL_ARCH_CUDA
    cudaStreamDestroy(stream);
    #elif YAKL_ARCH_HIP
    hipStreamDestroy(stream);
    #else
    stream = NULL;
    #endif
  }

void yakl_event_create(yakl::yakl_event_t * event)
  {
    #ifdef YAKL_ARCH_CUDA
    cudaEventCreate(event);
    #elif YAKL_ARCH_HIP
    hipEventCreate(event);
    #else
    stream = nullptr;
    #endif
  }

void yakl_event_destroy(yakl::yakl_event_t event)
  {
    #ifdef YAKL_ARCH_CUDA
    cudaEventDestroy(event);
    #elif YAKL_ARCH_HIP
    hipEventDestroy(event);
    #else
    stream = NULL;
    #endif
  }

void yakl_record_event(yakl::yakl_event_t event, yakl::yakl_stream_t stream)
{
  #ifdef YAKL_ARCH_CUDA
    cudaEventRecord(event, stream);
  #elif defined(YAKL_ARCH_HIP)
    hipEventRecord(event, stream);
  #endif
}

void yakl_stream_wait(yakl::yakl_stream_t stream, yakl::yakl_event_t event)
{
  #ifdef YAKL_ARCH_CUDA
    cudaStreamWaitEvent(stream, event, 0);
  #elif defined(YAKL_ARCH_HIP)
    hipStreamWaitEvent(stream, event, 0);
  #endif
}

