#include <exception>

#include "BoostThreadWorker.hpp"

namespace Slic3r { namespace GUI {

void BoostThreadWorker::WorkerMessage::deliver(BoostThreadWorker &runner)
{
    switch(MsgType(m_data.index())) {
    case Empty: break;
    case Status: {
        StatusInfo info = std::get<Status>(m_data);
        if (runner.get_pri()) {
            runner.get_pri()->set_progress(info.status);
            runner.get_pri()->set_status_text(info.msg.c_str());
        }
        break;
    }
    case Finalize: {
        JobEntry& entry = std::get<Finalize>(m_data);
        entry.job->finalize(entry.canceled, entry.eptr);

        // Unhandled exceptions are rethrown without mercy.
        if (entry.eptr)
            std::rethrow_exception(entry.eptr);

        break;
    }
    case MainThreadCall: {
        MainThreadCallData &calldata = std::get<MainThreadCall>(m_data);
        calldata.fn();
        calldata.promise.set_value();

        break;
    }
    }
}

void BoostThreadWorker::run()
{
    bool stop = false;
    while(!stop) {
        m_input_queue.consume_one_blk([this, &stop](JobEntry &e) {
            if (!e.job)
                stop = true;
            else {
                m_canceled.store(false);
                m_running.store(true);

                try {
                    e.job->process(*this);
                }  catch (...) {
                    e.eptr = std::current_exception();
                }

                e.canceled = m_canceled.load();
                m_output_queue.push(std::move(e)); // finalization message
                m_running.store(false);
            }
        });
    };
}

void BoostThreadWorker::update_status(int st, const std::string &msg)
{
    m_output_queue.push(st, msg);
}

std::future<void> BoostThreadWorker::call_on_main_thread(std::function<void ()> fn)
{
    MainThreadCallData cbdata{std::move(fn), {}};
    std::future<void> future = cbdata.promise.get_future();

    m_output_queue.push(std::move(cbdata));

    return future;
}

BoostThreadWorker::BoostThreadWorker(std::shared_ptr<ProgressIndicator> pri,
                                     boost::thread::attributes &attribs,
                                     const char *               name)
    : m_progress(std::move(pri)), m_name{name}
{
    if (m_progress)
        m_progress->set_cancel_callback([this](){ cancel(); });

    m_thread = create_thread(attribs, [this] { this->run(); });

    std::string nm{name};
    if (!nm.empty()) set_thread_name(m_thread, name);
}

constexpr int ABORT_WAIT_MAX_MS = 10000;

BoostThreadWorker::~BoostThreadWorker()
{
    try {
        cancel_all();
        m_input_queue.push(JobEntry{nullptr});
        join(ABORT_WAIT_MAX_MS);
    } catch(...) {
        BOOST_LOG_TRIVIAL(error)
            << "Could not join worker thread '" << m_name << "'";
    }
}

bool BoostThreadWorker::join(int timeout_ms)
{
    if (!m_thread.joinable())
        return true;

    if (timeout_ms <= 0) {
        m_thread.join();
    }
    else if (m_thread.try_join_for(boost::chrono::milliseconds(timeout_ms))) {
        return true;
    }
    else
        return false;

    return true;
}

void BoostThreadWorker::process_events()
{
    while (m_output_queue.consume_one([this](WorkerMessage &msg) {
        msg.deliver(*this);
    }));
}

bool BoostThreadWorker::start_next(std::unique_ptr<Job> job)
{
    if (job)
        m_input_queue.push(JobEntry{std::move(job)});

    return bool{job};
}

}} // namespace Slic3r::GUI
