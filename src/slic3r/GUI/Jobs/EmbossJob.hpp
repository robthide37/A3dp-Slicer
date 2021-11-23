#ifndef slic3r_EmbossJob_hpp_
#define slic3r_EmbossJob_hpp_

#include <memory>
#include <thread>
#include <mutex>
#include "Job.hpp"

#include "libslic3r/Emboss.hpp"

namespace Slic3r {
class ModelVolume;
class ModelObject;
}

namespace Slic3r::GUI {

// thread process of conversion from text to model 
// [optionaly] union with model
class EmbossJob : protected Job
{
public:
    struct Data
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
        int          object_idx;
    };
    EmbossJob();    
    ~EmbossJob();// cancel work and join

    void restart(const Data &data);

    // do not allow Job::start
    bool join(int timeout_ms = 0) { return Job::join(timeout_ms); };
    bool is_running() const { return Job::is_running(); }
    void cancel() { Job::cancel(); }

protected:
    // Launched just before start(), a job can use it to prepare internals
    virtual void prepare() override;

    // The method where the actual work of the job should be defined.
    virtual void process() override;

    // Launched when the job is finished. It refreshes the 3Dscene by def.
    virtual void finalize() override;

private:
    std::optional<Data> m_data;
    std::optional<indexed_triangle_set> m_result;

    std::mutex          m_mutex; // protect next data
    std::optional<Data> m_data_next;
    std::thread         m_restart_thread;

    // TODO: move to objec list utils
    void select_volume(ModelVolume * volume);


    class Progress : public ProgressIndicator
    {
        // Inherited via ProgressIndicator
        virtual void set_range(int range) override {}
        virtual void set_cancel_callback(CancelFn = CancelFn()) override {}
        virtual void set_progress(int pr) override {}
        virtual void set_status_text(const char *) override {}
        virtual int  get_range() const override { return 100; }
    };
};
} // namespace Slic3r::GUI

#endif // slic3r_EmbossJob_hpp_
