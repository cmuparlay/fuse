
# Fuse : A library for fusing lock based code into transactions

See:

Guy E. Blelloch, Zachary Kend and Yuanhao Wei\
TLF: Transactional Lock Fusion\
SPAA 2025

## Running the Experiments from the Paper

To test and run benchmarks:

```
mkdir build
cd build
cmake ..
cd benchmarks/transactions
make -j <numprocs that won't kill your machine>
./runtrans
```

This will write the timings to a file in the directory `timings`.  Alternatively
you can run individual benchmarks, e.g.:

```
./btree_fuse_tlf -n 100000,10000000 -u 0,5,50 -z 0,.99 -trans 1,16
```

will run the fuse version of btree on all combinations of tree sizes
of 100K and 10M, update rates of 0, 5 and 50%, zipfian parameters 0
and .99, and transaction sizes of 1 and 16 operations, reporting the
individual times and the geometric mean.

## Library Interface

The library can either used as a traditional stm incorporating
the TLF data structures or used to build your own TLF structures
The interface is hearer only.

```
  #include "fuse/fuse.h"
```

If used as an a traditional STM the relevant interface is:

```
namespace fuse {
  template <typename T>
  T* New(args...);     // use like 'new T(args...)'

  template <typename T>
  void Retire(T* x);  //  use like 'delete x'

  // f is a functor with no argument and returning T
  template <typename T, typename F>
  T atomic_region(F f);

  template <typename T, typename F>
  T atomic_read_only(F f);

  template <typename T>
  struct tlf_atomic<T> { // use like std::atomic
    atomic(T v);            // constructor
    T load();               // load contents
    void store(T val);      // store contents
    T operator=(T val);     // assignment
  }
}
```

Any allocation or freeing within an atomic region should be performed
with `New` and `Retire`, and any loads and stores of shared state
should use `tlf_atomic`.  In addition we supply a set of data
structures.  They can be used by including the relevant header file
For example:

```
#include "fuse/structs/btree.h"
fuse::tlf_btree_map<std::string, int> map1;
fuse::tlf_btree_map<std::string, int> map2;
fuse::tlf_atomic<int> counter;
...
int result = fuse::atomic_region([&] {
  std::string x = ...
  int cnt = counter.load();
  if (map1.find(x).has_value() &&
      !map2.find(x).has_value()) {
  map2.insert(x, cnt);
  else counter.store(cnt + 1);
  return cnt;
});
...
```

In designing ones own TLF concurrent data structures based on
optimistic locking, the following additional interface is useful:

```
namespace fuse {
  template <typename T, typename F>
  T with_epoch(F f);

  template <typename T>
  struct atomic<T> { // use like std::atomic
    atomic(T v);            // constructor
    T load();               // load contents
    void store(T val);      // store contents
    T operator=(T val);     // assignment
  }

  struct shared_mutex; // use like std::shared_mutex
}
```

Note that this `atomic` is not the same as `tlf_atomic`.  It is designed
for implementing concurrent data strcutures.  In particular if not inside
of a lock it is unprotected and will not ensure serializability (see the paper).
Indeed `tlf_atomic` is implemented with `atomic` and locks as follows:

```
template <typename T>
struct tlf_atomic {
  fuse::atomic<T> v;
  fuse::shared_mutex lock;

  T load() {
      std::shared_lock lck(lock);
      return v.load(); }

  void store(T x) {
      std::unique_lock lck(lock);
      v.store(x);}};
```

The `with_epoch` funtion is for memory protection, which is
implemented with epoch-based reclamation.  It should be wrapped around
any concurren functions if the memory it accesses can be freed
concurrently.    `atomic_region` and `atomic_read_only` include their own
'with-epoch`, so you don't need to uses this inside atomic regions.

See `structures/tlf_leaftree/ordered_map.h` for an example.  You can run this
by making as described above and using, e.g.:
```
./leaftree_fuse -n 100000,10000000 -u 0,5,50 -z 0,.99 -trans 1,16
```

All other structures in `structures/` are carried over from `verlib`
and are therefore written in a slightly different style (e.g. using
verlib namespace, and verlib style locks).
