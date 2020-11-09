#ifndef PARALLEL_HPP
#define PARALLEL_HPP

#include <algorithm>
#include <future>

namespace par {

/* Transform each element between beg and end with f and return output as std::vector. */
template <class RandomIt, class Func>
auto for_each(size_t min_len, RandomIt beg, RandomIt end, Func f) {
  using T = decltype(f(*beg));
  
  auto len = static_cast<size_t>(end - beg); // assert end > beg
  
  if (len < min_len) {
    std::vector<T> out(len);
    std::transform(beg, end, out.begin(), f);
    return out;
  }
  
  RandomIt mid = beg + len / 2;
  auto handle = std::async(std::launch::async, for_each<RandomIt, Func>,
                           min_len, mid, end, f);
  
  std::vector<T> head = for_each(min_len, beg, mid, f);
  std::vector<T> tail = handle.get();
  
  head.insert(head.end(), tail.begin(), tail.end());
  return head;
}

}  // namespace par

#endif  // PARALLEL_HPP