#ifndef FUSE_SKIPLIST_H
#define FUSE_SKIPLIST_H
  
namespace fuse {
#include "../../structures/skiplist/ordered_map.h"
  template <typename K,
            typename V,
            typename Compare = std::less<K>>
  using tlf_skiplist_map = verlib::skiplist<K,V,Compare>;
}

#endif
