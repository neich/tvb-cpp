//
// Created by imartin on 30-Jul-22.
//

#ifndef TVB_CPP_THREADPOOL_H
#define TVB_CPP_THREADPOOL_H

#include <thread>
#include <condition_variable>
#include <functional>
#include <queue>
#include <optional>

namespace tvb {

    template<typename R>
    class ThreadPool {
    public:
        explicit ThreadPool(int mt = 1) : max_threads(mt) {}


        void start() {
            threads.resize(max_threads);
            for (uint32_t i = 0; i < max_threads; i++) {
                threads.at(i) = std::thread(&ThreadPool<R>::ThreadLoop, this);
            }
        }

        void queue_job(const std::function<R()> &job) {
            {
                std::unique_lock<std::mutex> lock(queue_mutex);
                jobs.push(job);
                pending++;
            }
            mutex_condition.notify_one();
        }

        void stop() {
            {
                std::unique_lock<std::mutex> lock(queue_mutex);
                should_terminate = true;
            }
            mutex_condition.notify_all();
            for (std::thread &active_thread: threads) {
                active_thread.join();
            }
            threads.clear();
        }

        bool finished() {
            return pending == 0 && !has_results();
        }

        std::optional<R> get_result() {
            bool is_result = false;
            R r;
            {
                std::unique_lock<std::mutex> lock(queue_mutex);
                if (!results.empty()) {
                    is_result = true;
                    r = results.back();
                    results.pop_back();
                }
            }
            if (is_result)
                return r;
            else
                return {};
        }

        bool has_results() {
            bool has_results;
            {
                std::unique_lock<std::mutex> lock(queue_mutex);
                has_results = !results.empty();
            }
            return has_results;
        }


    private:
        void ThreadLoop() {
            while (true) {
                std::function<R()> job;
                {
                    std::unique_lock<std::mutex> lock(queue_mutex);
                    mutex_condition.wait(lock, [this] {
                        return !jobs.empty() || should_terminate;
                    });
                    if (should_terminate) {
                        return;
                    }
                    job = jobs.front();
                    jobs.pop();
                }
                R r = job();
                {
                    std::unique_lock<std::mutex> lock(queue_mutex);
                    results.push_back(r);
                    pending--;
                }
            }
        }

        int max_threads;
        bool should_terminate = false;           // Tells threads to stop looking for jobs
        std::mutex queue_mutex;                  // Prevents data races to the job queue
        std::condition_variable mutex_condition; // Allows threads to wait on new jobs or termination
        std::vector<std::thread> threads;
        std::queue<std::function<R()>> jobs;
        std::vector<R> results;
        int pending = 0;

    };
}

#endif //TVB_CPP_THREADPOOL_H
