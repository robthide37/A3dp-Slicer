#ifndef PLATERJOB_HPP
#define PLATERJOB_HPP

#include "BusyCursorJob.hpp"

#include "slic3r/GUI/GUI.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/I18N.hpp"
#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"

namespace Slic3r { namespace GUI {

class Plater;

template<class WorkerSubclass>
class PlaterWorker: public Worker {
    WorkerSubclass m_w;

    class PlaterJob : public Job {
        std::unique_ptr<Job> m_job;
        Plater *m_plater;

    public:
        void process(Ctl &c) override
        {
            CursorSetterRAII busycursor{c};
            m_job->process(c);
        }

        void finalize(bool canceled, std::exception_ptr &eptr) override
        {
            m_job->finalize(canceled, eptr);

            if (eptr) try {
                    std::rethrow_exception(eptr);
                }  catch (std::exception &e) {
                    show_error(m_plater, _L("An unexpected error occured: ") + e.what());
                    eptr = nullptr;
                }
        }

        PlaterJob(std::unique_ptr<Job> j)
            : m_job{std::move(j)}, m_plater{wxGetApp().plater()}
        {
            // TODO: decide if disabling slice button during UI job is what we want.
            //        if (m_plater)
            //            m_plater->sidebar().enable_buttons(false);
        }

        ~PlaterJob() override
        {
            // TODO: decide if disabling slice button during UI job is what we want.

            // Reload scene ensures that the slice button gets properly
            // enabled or disabled after the job finishes, depending on the
            // state of slicing. This might be an overkill but works for now.
            //        if (m_plater)
            //            m_plater->canvas3D()->reload_scene(false);
        }
    };

public:

    template<class ... WorkerArgs>
    PlaterWorker(WorkerArgs &&...args) : m_w{std::forward<WorkerArgs>(args)...} {}

    // Always package the job argument into a PlaterJob
    bool start_next(std::unique_ptr<Job> job) override
    {
        return m_w.start_next(std::make_unique<PlaterJob>(std::move(job)));
    }

    bool is_idle() const override { return m_w.is_idle(); }
    void cancel() override { m_w.cancel(); }
    void cancel_all() override { m_w.cancel_all(); }
    void process_events() override { m_w.process_events(); }
};

}} // namespace Slic3r::GUI

#endif // PLATERJOB_HPP
