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

template<class JobSubclass>
class PlaterJob : public Job {
    JobSubclass m_job;
    Plater *m_plater;

public:

    void process(Ctl &c) override
    {
        CursorSetterRAII busycursor{c};
        m_job.process(c);
    }

    void finalize(bool canceled, std::exception_ptr &eptr) override
    {
        m_job.finalize(canceled, eptr);

        if (eptr) try {
            std::rethrow_exception(eptr);
        }  catch (std::exception &e) {
            show_error(m_plater, _L("An unexpected error occured: ") + e.what());
            eptr = nullptr;
        }
    }

    template<class... Args>
    PlaterJob(Args &&...args)
        : m_job(std::forward<Args>(args)...), m_plater{wxGetApp().plater()}
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

template<class JobSubclass, class ... Args>
void replace_job(Plater& p, Args && ...args)
{
    replace_job(p.get_ui_job_worker(),
                std::make_unique<PlaterJob<JobSubclass>>(
                    std::forward<Args>(args)...));
}

}} // namespace Slic3r::GUI

#endif // PLATERJOB_HPP
