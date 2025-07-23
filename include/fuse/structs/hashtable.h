#ifndef FUSE_HASHTABLE_H
#define FUSE_HASHTABLE_H

namespace fuse {
#include "../../structures/hash_block/unordered_map.h"
  template <typename K,
            typename V,
            class Hash = std::hash<K>,
            class KeyEqual = std::equal_to<K>>
  using tlf_hashtable_map = verlib::hash_block<K,V,Hash,KeyEqual>;
}
#endif
