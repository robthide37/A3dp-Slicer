#include "catch2/catch.hpp"

#include <chrono>
#include <thread>

#include "slic3r/GUI/Jobs/BoostThreadWorker.hpp"
#include "slic3r/GUI/Jobs/ProgressIndicator.hpp"

//#include <boost/thread/thread.hpp>

struct Progress: Slic3r::ProgressIndicator {
    int range = 100;
    int pr = 0;
    std::string statustxt;
    void set_range(int r) override { range = r; }
    void set_cancel_callback(CancelFn = CancelFn()) override {}
    void set_progress(int p) override { pr = p; }
    void set_status_text(const char *txt) override { statustxt = txt; }
    int  get_range() const override { return range; }
};

TEST_CASE("nullptr job should be ignored", "[Jobs]") {
    Slic3r::GUI::BoostThreadWorker worker{std::make_unique<Progress>()};
    worker.push(nullptr);

    REQUIRE(worker.is_idle());
}

TEST_CASE("State should not be idle while running a job", "[Jobs]") {
    using namespace Slic3r;
    using namespace Slic3r::GUI;
    BoostThreadWorker worker{std::make_unique<Progress>(), "worker_thread"};

    queue_job(worker, [&worker](Job::Ctl &ctl) {
        ctl.call_on_main_thread([&worker] {
            REQUIRE(!worker.is_idle());
        }).wait();
    });

    worker.wait_for_idle();

    REQUIRE(worker.is_idle());
}

TEST_CASE("Status messages should be received by the main thread during job execution", "[Jobs]") {
    using namespace Slic3r;
    using namespace Slic3r::GUI;
    auto pri = std::make_shared<Progress>();
    BoostThreadWorker worker{pri};

    queue_job(worker, [](Job::Ctl &ctl){
        for (int s = 0; s <= 100; ++s) {
            ctl.update_status(s, "Running");
        }
    });

    worker.wait_for_idle();

    REQUIRE(pri->pr == 100);
    REQUIRE(pri->statustxt == "Running");
}

TEST_CASE("Cancellation should be recognized be the worker", "[Jobs]") {
    using namespace Slic3r;
    using namespace Slic3r::GUI;

    auto pri = std::make_shared<Progress>();
    BoostThreadWorker worker{pri};

    queue_job(
        worker,
        [](Job::Ctl &ctl) {
            for (int s = 0; s <= 100; ++s) {
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
                ctl.update_status(s, "Running");
                if (ctl.was_canceled()) break;
            }
        },
        [](bool cancelled, std::exception_ptr &) { // finalize
            REQUIRE(cancelled == true);
        });

    std::this_thread::sleep_for(std::chrono::milliseconds(1));
    worker.cancel();

    worker.wait_for_current_job();

    REQUIRE(pri->pr != 100);
}

TEST_CASE("cancel_all should remove all pending jobs", "[Jobs]") {
    using namespace Slic3r;
    using namespace Slic3r::GUI;

    auto pri = std::make_shared<Progress>();
    BoostThreadWorker worker{pri};

    std::array<bool, 4> jobres = {false};

    queue_job(worker, [&jobres](Job::Ctl &) {
        jobres[0] = true;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    });
    queue_job(worker, [&jobres](Job::Ctl &) {
        jobres[1] = true;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    });
    queue_job(worker, [&jobres](Job::Ctl &) {
        jobres[2] = true;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    });
    queue_job(worker, [&jobres](Job::Ctl &) {
        jobres[3] = true;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    });

    std::this_thread::sleep_for(std::chrono::microseconds(500));
    worker.cancel_all();

    REQUIRE(jobres[0] == true);
    REQUIRE(jobres[1] == false);
    REQUIRE(jobres[2] == false);
    REQUIRE(jobres[3] == false);
}

TEST_CASE("Exception should be properly forwarded to finalize()", "[Jobs]") {
    using namespace Slic3r;
    using namespace Slic3r::GUI;

    auto pri = std::make_shared<Progress>();
    BoostThreadWorker worker{pri};

    queue_job(
        worker, [](Job::Ctl &) { throw std::runtime_error("test"); },
        [](bool /*canceled*/, std::exception_ptr &eptr) {
            REQUIRE(eptr != nullptr);
            try {
                std::rethrow_exception(eptr);
            } catch (std::runtime_error &e) {
                REQUIRE(std::string(e.what()) == "test");
            }

            eptr = nullptr;
        });

    worker.wait_for_idle();
    REQUIRE(worker.is_idle());
}
