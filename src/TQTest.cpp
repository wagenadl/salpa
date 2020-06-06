#include <iostream>
#include "TaskQueue.h"
#include <future>
#include <chrono>
#include <thread>

int main() {
  TaskQueue<std::packaged_task<void()>> tq(8);
  for (int i=0; i<10; i++) {
    for (int j=0; j<10; j++) {
      std::cerr << "posting " << i << " + " << j << "\n";
      std::packaged_task<void()> task([i,j]() {
                std::this_thread::sleep_for(std::chrono::seconds(1));
                std::cerr << i << "+" << j << "=" << i+j << "\n";
              });
      tq.post(task);
    }
    std::cerr << "waiting\n";
    tq.wait();
  }
  return 0;
}
