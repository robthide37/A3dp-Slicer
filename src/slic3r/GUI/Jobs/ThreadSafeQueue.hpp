#ifndef THREADSAFEQUEUE_HPP
#define THREADSAFEQUEUE_HPP

#include <type_traits>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>

namespace Slic3r { namespace GUI {

// A thread safe queue for one producer and one consumer. Use consume_one_blk
// to block on an empty queue.
template<class T,
         template<class, class...> class Container = std::deque,
         class... ContainerArgs>
class ThreadSafeQueueSPSC
{
    std::queue<T, Container<T, ContainerArgs...>> m_queue;
    mutable std::mutex m_mutex;
    std::condition_variable m_cond_var;
public:

    // Consume one element, block if the queue is empty.
    template<class Fn> void consume_one_blk(Fn &&fn, std::atomic<bool> *pop_flag = nullptr)
    {
        static_assert(!std::is_reference_v<T>, "");
        static_assert(std::is_default_constructible_v<T>, "");
        static_assert(std::is_move_assignable_v<T> || std::is_copy_assignable_v<T>, "");

        T el;
        {
            std::unique_lock lk{m_mutex};
            m_cond_var.wait(lk, [this]{ return !m_queue.empty(); });

            if constexpr (std::is_move_assignable_v<T>)
                el = std::move(m_queue.front());
            else
                el = m_queue.front();

            m_queue.pop();

            if (pop_flag) // The optional atomic is set before the lock us unlocked
                pop_flag->store(true);
        }

        fn(el);
    }

    // Consume one element, return true if consumed, false if queue was empty.
    template<class Fn> bool consume_one(Fn &&fn)
    {
        T el;
        {
            std::unique_lock lk{m_mutex};
            if (!m_queue.empty()) {
                if constexpr (std::is_move_assignable_v<T>)
                    el = std::move(m_queue.front());
                else
                    el = m_queue.front();

                m_queue.pop();
            } else
                return false;
        }

        fn(el);

        return true;
    }

    // Push element into the queue.
    template<class...TArgs> void push(TArgs&&...el)
    {
        std::lock_guard lk{m_mutex};
        m_queue.emplace(std::forward<TArgs>(el)...);
        m_cond_var.notify_one();
    }

    bool empty() const
    {
        std::lock_guard lk{m_mutex};
        return m_queue.empty();
    }

    size_t size() const
    {
        std::lock_guard lk{m_mutex};
        return m_queue.size();
    }

    void clear()
    {
        std::lock_guard lk{m_mutex};
        while (!m_queue.empty()) m_queue.pop();
    }
};

}} // namespace Slic3r::GUI

#endif // THREADSAFEQUEUE_HPP
