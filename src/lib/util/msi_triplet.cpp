#include "msi_triplet.hpp"

#include "serialization.hpp"

namespace triplet {

constexpr size_t size_triplet =
    sizeof(uint32_t) + sizeof(double) + sizeof(float);
constexpr size_t size_uint = sizeof(uint64_t);

/* Binary search in file inspired by std::lower_bound. */
uint64_t lower_bound(std::istream& stream, double val) {
  stream.seekg(0); // reset stream

  uint64_t len;  // number of triplets
  Serialization::read_uint64(stream, &len);

  uint64_t step;
  uint64_t it;
  uint64_t first = 0;
  triplet t;

  while (len > 0) {
    it = first;
    step = len / 2;
    it += step;

    stream.seekg(size_uint + it * size_triplet);
    read_triplet(stream, &t);

    if (t.mz < val) {
      first = it + 1;
      len -= step + 1;
    } else {
      len = step;
    }
  }

  return first;
}

/* Retrieve all triplets within m/z range. */
std::vector<triplet> get_range(std::istream& stream, double begin, double end) {
  if (begin >= end) {
    return {};
  }

  auto first = lower_bound(stream, begin);
  auto last = lower_bound(stream, end);

  size_t n_triplets = last - first;
  std::vector<triplet> ts{n_triplets};

  stream.seekg(size_uint + size_triplet * first);
  for (size_t i = 0; i < n_triplets; ++i) {
    read_triplet(stream, &ts[i]);
  }

  return ts;
}

bool read_triplet(std::istream& stream, triplet* t) {
  Serialization::read_uint32(stream, &t->index);
  Serialization::read_double(stream, &t->mz);
  Serialization::read_float(stream, &t->height);
  return stream.good();
}

bool write_triplet(std::ostream& stream, const triplet& t) {
  Serialization::write_uint32(stream, t.index);
  Serialization::write_double(stream, t.mz);
  Serialization::write_float(stream, t.height);
  return stream.good();
}

bool read_triplets(std::istream& stream, std::vector<triplet>* v) {
  return Serialization::read_vector(stream, v, read_triplet);
}

bool write_triplets(std::ostream& stream, const std::vector<triplet>& v) {
  return Serialization::write_vector(stream, v, write_triplet);
}

}  // namespace triplet
