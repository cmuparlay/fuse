#include "ordered_map.h"

template <typename K,
	  typename V,
	  typename Less = std::less<K>>
using ordered_map = fuse::tlf_leaftree_map<K, V, Less>;
