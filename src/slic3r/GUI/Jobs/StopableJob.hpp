#ifndef slic3r_StopableJob_hpp_
#define slic3r_StopableJob_hpp_

#include <memory>
#include <thread>
#include <mutex>

namespace Slic3r::GUI {

// inspired by Job.hpp
// All public function can be call only from UI thread
// Mechanism to stack processing and do only last one
// Ability to stop processing developer must add check into m_func
using StopCondition = std::function<bool(void)>;
template<typename TIn> class StopableJob
{
    std::mutex           m_mutex;
    std::unique_ptr<TIn> m_input_next = nullptr;

    std::thread m_thread;
    // when running == true, than check of new input in thread.
    bool m_running = false;
    // faster interupt inside func, developer must add StopCondifion call
    bool m_stop = false;
    
    using Func = std::function<void(std::unique_ptr<TIn>, StopCondition)>;
    Func m_func;

public:
    StopableJob(Func func) : m_func(func) {}
    virtual ~StopableJob() { 
        stop();
        try {
            // thread join could throw exception
            // https://en.cppreference.com/w/cpp/thread/thread/join
            join();
        } catch (std::system_error err) {}
    }

    void run(std::unique_ptr<TIn> input)
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
                        m_func(std::move(input), [this]() { return is_stoping(); });

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

    bool is_running()
    {
        std::lock_guard lg(m_mutex);
        return m_running;
    }
    bool is_stoping()
    {
        std::lock_guard lg(m_mutex);
        return m_stop;
    }
    void stop() { 
        std::lock_guard lg(m_mutex);
        if (!m_running) return;
        m_input_next = nullptr;
        m_stop = true; 
    }
    // Be Carefull, blocking until join
    // call stop when you not sure
    void join(int timeout_ms = 0)
    {
        if (m_thread.joinable()) m_thread.join();
        assert(!m_running);
    }
};

} // namespace Slic3r::GUI

#endif // slic3r_StopableJob_hpp_
