
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

The interface is hearer only and can be included with:

```
  #include "fuse/fuse.h"
```

The interface supports the following:

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
  struct atomic<T> { // use like std::atomic
    atomic(T v);            // constructor
    T load();               // load contents
    void store(T val);      // store contents
    T operator=(T val);     // assignment
  }
  
  struct shared_mutex; // use like std::shared_mutex
}
```

See `structures/tlf_leaftree/ordered_map.h` for an example.  You can run this
by making as described above and using, e.g.:
```
./leaftree_fuse -n 100000,10000000 -u 0,5,50 -z 0,.99 -trans 1,16
```

All other structures in `structures/` are carried over from `verlib`
and are therefore written in a slightly different style (e.g. using
verlib namespace, and verlib style locks).
