#ifndef FUSE_ARTTREE_H
#define FUSE_ARTTREE_H
  
namespace fuse {
#include "../../structures/arttree/ordered_map.h"
  template <typename K,
          typename V,
          typename String = verlib::int_string<K>>
  using tlf_arttree_map = verlib::arttree<K,V,String>;
}

#endif
