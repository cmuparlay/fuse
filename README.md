
# Fuse : A library for fusing lock based code into transactions

See:

Guy E. Blelloch, Zachary Kent and Yuanhao Wei\
TLF: Transactional Lock Fusion\
SPAA 2025

## Running the Experiments from the Paper

To test and run benchmarks:

```
mkdir build
cd build
cmake ..
cd benchmarks/transactions
make -j <at most num cores on your machine>
./runtrans
```

This will write the timings to a file in the directory `timings`.  Alternatively
you can run individual benchmarks, e.g.:

```
./btree_mv_tlf -n 100000,10000000 -u 0,5,50 -z 0,.99 -trans 1,16
```

will run the multiversion TLF version of btree on all combinations of
tree sizes of 100K and 10M, update rates of 0, 5 and 50%, zipfian
parameters 0 and .99, and transaction sizes of 1 and 16 operations,
reporting the individual times and the geometric mean.

The naming convention is `<struct>_<stm_system>_<variant>` where

* `<struct>` is one of `arttree`, `btree`, `hash_block`, `skiplist`,
`leaftree`, `avltree`, `treap` or `list`, and
* `<stm_system>` is one of `mv` (multiversion) or `2plsf`,
* and `<variant>` is one of
  * `stm` : traditional stm
  * `tlf` : TLF with coarse grained locks (original locks in the optimistic locking code)
  * `tlf_stm` : TLF with fine grained locks (locks for each write)

Furthermore we have `<struct>` for each structure, which is just the lock-based
concurrent data structure and `<struct>_versioned` which is the original
data structure with multiversioning (using verlib).

## Library Interface

The library can either used as a traditional Software Transactional
Memory (STM) incorporating the TLF data structures or used to build
your own Optimistic Locking structures using TLF.
The interface is header only.

```
  #include "fuse/fuse.h"
```

If used as an a traditional STM, the relevant interface is:

```
namespace fuse {
  // f is a functor with no argument and returning T
  template <typename T, typename F>
  T atomic_region(F f);

  // can be used if f only does loads
  template <typename T, typename F>
  T atomic_read_only(F f);

  // use like std::atomic
  template <typename T>
  struct tlf_atomic<T> { 
    atomic(T v);            // constructor
    T load();               // load contents
    void store(T val);      // store contents
    T operator=(T val);     // assignment
  }

  // use like 'new T(args...)'
  template <typename T>
  T* New(args...);     

  // use like 'delete x'
  template <typename T>
  void Retire(T* x);  
}
```

Any loads and stores of shared state should use `tlf_atomic` and any
allocation or freeing within an atomic region should be performed with
`New` and `Retire`.  In addition, we supply a set of TLF data
structures.  These are all based on the
[Verlib](https://github.com/cmuparlay/verlib) concurrent data
structures, which all use optimistic locking.  The TLF structures can
be used by including the relevant header file.  For example:

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
Other structures include `arttree`, `hashtable`, `leaftree`, `treap`, `skiplist`, and `ordered_list`.
Another more detailed example is given [here](examples/move.cpp).

In designing your own optimistic locking data structures with TLF, the
following additional interface can be used:

```
namespace fuse {
  // use like std::atomic
  template <typename T>
  struct atomic<T> { 
    atomic(T v);            // constructor
    T load();               // load contents
    void store(T val);      // store contents
    T operator=(T val);     // assignment
  }

  // use like std::shared_mutex
  struct shared_mutex; 

  template <typename T, typename F>
  T with_epoch(F f);
}
```

Note that this `atomic` is not the same as `tlf_atomic`.  It is designed
for implementing concurrent data structures.  In particular if not inside
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

The `with_epoch` function is for memory protection, which is
implemented with epoch-based reclamation.  It should be wrapped around
any concurrent functions if the memory it accesses can be freed
concurrently.    `atomic_region` and `atomic_read_only` include their own
'with-epoch`, so you don't need to uses this inside atomic regions.

See [`tlf_leaftree`](include/structures/tlf_leaftree/ordered_map.h) for an example.
This is basically the example from the paper.
You can run this
by making as described above and using, e.g.:
```
./leaftree_fuse -n 100000,10000000 -u 0,5,50 -z 0,.99 -trans 1,16
```

All other structures in `structures/` are carried over from `verlib`
and are therefore written in a slightly different style (e.g. using
verlib namespace, and verlib style locks).
