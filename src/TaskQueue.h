// taskqueue.h

#ifndef TASKQUEUE_H

#define TASKQUEUE_H

#include <mutex>
#include <thread>
#include <vector>
#include <queue>
#include <future>
#include <condition_variable>

template <class F> class TaskQueue {
public:
  TaskQueue(int n=0) {
    abort = false;
    nactive = 0;
    if (n==0)
      n = std::thread::hardware_concurrency();
    if (n<2)
      n = 2;
    for (int k=0; k<n; k++)
      thrs.push_back(std::thread(&TaskQueue::worker, this));
  }
  ~TaskQueue() {
    abort = true;
    cond.notify_all();
    for (auto &t: thrs)
      t.join();
  }
  void post(F &request) {
    std::lock_guard<std::mutex> lock(mut);
    queue.push(std::move(request));
    cond.notify_one();
  }
  void wait() {
    std::unique_lock<std::mutex> lock(mut);
    backcond.wait(lock, [&]() {
                          return nactive==0 && queue.empty();
                        });
  }
  void worker() {
    std::unique_lock<std::mutex> lock(mut);
    while (!abort) {
      cond.wait(lock, [&]() {
                        return abort || !queue.empty();
                      });
      if (abort)
        break;
      F req(std::move(queue.front()));
      queue.pop();
      nactive ++;
      lock.unlock();
      req();
      lock.lock();
      nactive --;
      backcond.notify_one();
    }
    while (!queue.empty())
      queue.pop();
  }
    
private:
  std::condition_variable_any cond;
  std::condition_variable_any backcond;
  std::mutex mut;
  std::vector<std::thread> thrs;
  std::queue<F> queue;
  int nactive;
  bool abort;
};

#endif
