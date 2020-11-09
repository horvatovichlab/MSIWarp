#ifndef MSI_SBD_HPP
#define MSI_SBD_HPP

#include <cstdint>
#include <iostream>
#include <vector>

namespace sbd {

// "<BQB"
struct header {
    uint8_t id;
    uint64_t n_spectra;
    uint8_t mz_size;
};

// "<QLfHH"
struct meta {
    uint64_t offset;
    uint32_t n_points;
    float tic;
    uint16_t x;
    uint16_t y;
};

//
struct spectrum {
  std::vector<double> mzs;
  std::vector<float> heights;
};

/* read/write header from/to binary stream */
bool read_header(std::istream&, header*);
bool write_header(std::ostream&, const header&);

/* read/write meta from/to binary stream */
bool read_meta(std::istream&, meta*);
bool write_meta(std::ostream&, const meta&);

/* read/write vector of meta items from/to binary stream */
bool read_meta_items(std::istream&, std::vector<meta>*);
bool write_meta_items(std::ostream&, const std::vector<meta>&);

/* read/write spectrum from/to binary stream */
bool read_spectrum(std::istream&, spectrum*, size_t);
bool write_spectrum(std::ostream&, spectrum*);


} // namespace sbd

#endif // MSI_SBD_HPP


