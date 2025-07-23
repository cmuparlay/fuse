#ifndef FUSE_TREAP_H
#define FUSE_TREAP_H
  
namespace fuse {
#include "../../structures/treap/ordered_map.h"
  template <typename K,
            typename V,
            typename Compare = std::less<K>>
  using tlf_treap_map = verlib::treap<K,V,Compare>;
}

#endif
