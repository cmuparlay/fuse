#ifndef FUSE_LIST_H
#define FUSE_LIST_H
  
namespace fuse {
#include "../../structures/list/ordered_map.h"
  template <typename K,
            typename V,
            typename Compare = std::less<K>>
  using tlf_ordered_list_map = verlib::list<K,V,Compare>;
}

#endif
