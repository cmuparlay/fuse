#ifndef FUSE_BTREE_H
#define FUSE_BTREE_H

namespace fuse {
#include "../../../structures/btree/ordered_map.h"
  template <typename K,
            typename V,
            typename Compare = std::less<K>>
  using tlf_btree_map = verlib::btree<K,V,Compare>;
}

#endif
