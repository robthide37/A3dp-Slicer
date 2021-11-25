#ifndef slic3r_StopableJob_hpp_
#define slic3r_StopableJob_hpp_

#include <memory>
#include <thread>
#include <mutex>
#include <assert.h>

namespace Slic3r::GUI {

// inspired by Job.hpp
// Mechanism to stack processing and do only last one
// Ability to stop processing developer must add check is_stopping() into process()
template<typename TIn> class StopableJob
{
    std::mutex           m_mutex;
    std::unique_ptr<TIn> m_input_next = nullptr;

    std::thread m_thread;
    // when running == true, than check of new input in thread.
    bool m_running = false;
    // faster interupt inside func, developer must add StopCondifion call
    bool m_stop = false;
public:
    /// <summary>
    /// Stop and join thread
    /// </summary>
    virtual ~StopableJob();

    /// <summary>
    /// Restart processing of input
    /// </summary>
    /// <param name="input">Data needed to process</param>
    void run(std::unique_ptr<TIn> input);    

    /// <summary>
    /// Check if job processing input now
    /// </summary>
    /// <returns>True when run now otherwise False</returns>
    bool is_running(); // const; -- mutex

    /// <summary>
    /// Check if actual job is stopping now
    /// </summary>
    /// <returns>True when at stop process otherwise False</returns>
    bool is_stoping(); // const; -- mutex

    /// <summary>
    /// set flag to stop processing
    /// </summary>
    void stop();

    /// <summary>
    /// Free thread resources by join thread
    /// Be Carefull, it is blocking until join
    /// Suggest to call stop() before join
    /// Thread join could throw exception
    /// https://en.cppreference.com/w/cpp/thread/thread/join
    /// </summary>
    void join();
protected:

    /// <summary>
    /// Thread job of processing input data
    /// Note: Use check function is_stoping(), when true than interupt processing
    /// </summary>
    /// <param name="input">input data to process</param>
    virtual void process(std::unique_ptr<TIn> input) = 0;
};

//////
// Implementation 
//////
template<typename TIn> 
StopableJob<TIn>::~StopableJob()
{
    stop();
    try {
        // thread join could throw exception
        // https://en.cppreference.com/w/cpp/thread/thread/join
        join();
    } catch (std::system_error err) {}
}

template<typename TIn> 
void StopableJob<TIn>::run(std::unique_ptr<TIn> input)
{
    if (input == nullptr) return;
    {
        // keep access to next input
        std::lock_guard lg(m_mutex);
        if (m_running) {
            // when runnig
            m_stop       = true;
            m_input_next = std::move(input);
            return; // on end of run will be used new input
        }
        m_running = true;
        m_stop    = false;
    }
    if (m_thread.joinable()) m_thread.join();
    // at this moment is not running --> stoped
    assert(m_input_next == nullptr);
    try { // Execute the job
        m_thread = std::thread(
            [this](std::unique_ptr<TIn> input) {
                do {
                    process(std::move(input));

                    std::lock_guard lg(m_mutex);
                    m_stop = false;
                    // this is not while (end)condition because of lock guard
                    if (m_input_next == nullptr) {
                        m_running = false;
                        return;
                    }
                    input        = std::move(m_input_next);
                    m_input_next = nullptr;
                } while (true);
            },
            std::move(input));
    } catch (std::exception &) {}
}

template<typename TIn> 
bool StopableJob<TIn>::is_running()
{
    std::lock_guard lg(m_mutex);
    return m_running;
}


template<typename TIn>
bool StopableJob<TIn>::is_stoping()
{
    std::lock_guard lg(m_mutex);
    return m_stop;
}

template<typename TIn>
void StopableJob<TIn>::stop()
{
    std::lock_guard lg(m_mutex);
    if (!m_running) return;
    m_input_next = nullptr;
    m_stop       = true;
}

template<typename TIn>
void StopableJob<TIn>::join()
{
    if (m_thread.joinable()) m_thread.join();
    assert(!m_running);
}

} // namespace Slic3r::GUI

#endif // slic3r_StopableJob_hpp_
