#ifndef MSI_TRIPLET
#define MSI_TRIPLET

#include <iostream>
#include <vector>

namespace triplet {

/* An MSI data point composed of a spectrum id (index), m/z location, and
 * intensity (height). */
struct triplet {
  uint32_t index;
  double mz;
  float height;
};

/* Take all triplets within m/z range. Assumes file is sorted by m/z */
std::vector<triplet> get_range(std::istream& stream, double begin, double end);

/* Read/write triplet to/from binary stream. */
bool read_triplet(std::istream&, triplet*);
bool write_triplet(std::ostream&, const triplet&);

/* Read/write a vector of triplets to/from binary stream. */
bool read_triplets(std::istream&, std::vector<triplet>*);
bool write_triplets(std::ostream&, const std::vector<triplet>&);

}  // namespace triplet

#endif  // MSI_TRIPLET