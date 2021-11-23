#ifndef slic3r_EmbossJob_hpp_
#define slic3r_EmbossJob_hpp_

#include <memory>
#include <thread>
#include <mutex>
#include "Job.hpp"

#include "libslic3r/Emboss.hpp"

namespace Slic3r {
class ModelVolume;
}

namespace Slic3r::GUI {

// inspired by Job.hpp
// All public function can be call only from UI thread
// Mechanism to stack processing and do only last one
template<typename TIn> class ReRunJob
{
    std::mutex           m_next_mutex;
    std::unique_ptr<TIn> m_input_next = nullptr;

    std::thread m_thread;
    // indicate checking of new input.
    std::atomic<bool> m_running{false};

    using Func = std::function<void(std::unique_ptr<TIn>)>;
    Func m_func;

public:
    ReRunJob(Func func) : m_func(func) {}
    virtual ~ReRunJob() { join(); }

    void re_run(std::unique_ptr<TIn> input)
    {
        if (input == nullptr) return;
        {
            // keep access to next input
            std::lock_guard lg(m_next_mutex);
            if (is_running()) {
                // when runnig
                m_input_next = std::move(input);
                return; // on end of run will be used new input
            }
            m_running.store(true);
        }
        if (m_thread.joinable()) m_thread.join();
        // at this moment is not running --> stoped
        assert(m_input_next == nullptr);
        try { // Execute the job
            m_thread = std::thread(
                [this](std::unique_ptr<TIn> input) {
                    do {
                        m_func(std::move(input));

                        std::lock_guard lg(m_next_mutex);
                        // it is not in while condition because of lock guard
                        if (m_input_next == nullptr) {
                            m_running.store(false);
                            return;
                        }
                        input        = std::move(m_input_next);
                        m_input_next = nullptr;
                    } while (true);
                },
                std::move(input));
        } catch (std::exception &) {}
    }

    bool is_running() const { return m_running.load(); }
    void join()
    {
        if (m_thread.joinable()) m_thread.join();
        assert(!is_running());
    }
};

struct EmbossData
{
    // Pointer on Data of font (glyph shapes)
    std::shared_ptr<Emboss::Font> font;
    // font item is not used for create object
    TextConfiguration text_configuration;
    // new volume name created from text
    std::string volume_name;
    // when volume_ptr == nullptr than new volume will be created
    ModelVolume *volume_ptr;
    // when volume_ptr == nullptr && object_idx < 0 than new object will be created
    int object_idx;
};

class EmbossJob : public ReRunJob<EmbossData>
{
public:
    EmbossJob();
};

} // namespace Slic3r::GUI

#endif // slic3r_EmbossJob_hpp_
