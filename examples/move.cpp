#include "fuse/fuse.h"
#include "fuse/structs/leaftree.h"
#include "fuse/structs/hashtable.h"
#include "parlay/parallel.h"
#include <iostream>

int main() {

  int n = 1000000;

  fuse::tlf_leaftree_map<int,int> a;
  fuse::tlf_hashtable_map<int,int> b(n);
  fuse::tlf_atomic<int> c(0);

  // insert [0..n) into a in parallel
  parlay::parallel_for(0, n, [&] (long i) {
    fuse::atomic_region([&] { a.insert(i, i); });});

  // First n iterations move [0..n) from a to b in parallel.
  // Second n interations check that only in one at a time
  // If operations are not atomic then it could be in both
  // since we insert to b then remove from a.
  parlay::parallel_for(0, 2 * n, [&] (long i) {
    if (i < n) {                               
      fuse::atomic_region([&] {
        auto x = a.find(i);
        if (x.has_value()) {
          b.insert(i, *x);
          a.remove(i);
        }});
    } else {
      fuse::atomic_region([&] {
        if (a.find(i-n).has_value() == b.find(i-n).has_value())
          c.store(c.load() + 1); });
    }
  });

  if (c.load() == 0) std::cout << "I'm good" << std::endl;
  else std::cout << "Yikes!, failed atomicity: " << c.load() << std::endl;
}
