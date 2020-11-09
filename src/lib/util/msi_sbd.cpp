#include "msi_sbd.hpp"

#include "serialization.hpp"

namespace sbd {

constexpr size_t header_size =
    sizeof(uint8_t) + sizeof(uint64_t) + sizeof(uint8_t);
    
constexpr size_t meta_size = sizeof(uint64_t) + sizeof(uint32_t) +
                             sizeof(float) + sizeof(uint16_t) +
                             sizeof(uint16_t);

bool read_header(std::istream& stream, header* h) {
  Serialization::read_uint8(stream, &h->id);
  Serialization::read_uint64(stream, &h->n_spectra);
  Serialization::read_uint8(stream, &h->mz_size);
  return stream.good();
}

bool write_header(std::ostream& stream, const header& h) {
  Serialization::write_uint8(stream, h.id);
  Serialization::write_uint64(stream, h.n_spectra);
  Serialization::write_uint8(stream, h.mz_size);
  return stream.good();
}

bool read_meta(std::istream& stream, meta* m) {
  Serialization::read_uint64(stream, &m->offset);
  Serialization::read_uint32(stream, &m->n_points);
  Serialization::read_float(stream, &m->tic);
  Serialization::read_uint16(stream, &m->x);
  Serialization::read_uint16(stream, &m->y);
  return stream.good();
}

bool write_meta(std::ostream& stream, const meta& m) {
  Serialization::write_uint64(stream, m.offset);
  Serialization::write_uint32(stream, m.n_points);
  Serialization::write_float(stream, m.tic);
  Serialization::write_uint16(stream, m.x);
  Serialization::write_uint16(stream, m.y);
  return stream.good();
}

bool read_meta_items(std::istream& stream, std::vector<meta>* v) {
  for (auto& m : *v) {
      read_meta(stream, &m);
  }
  return stream.good();
}

bool write_meta_items(std::ostream& stream, const std::vector<meta>& v) {
    return Serialization::write_vector(stream, v, write_meta);
}

bool read_spectrum(std::istream& stream, spectrum* s, size_t n) {
  s->mzs.resize(n);
  s->heights.resize(n);

  for(auto& mz : s->mzs) {
    Serialization::read_double(stream, &mz);
  }

  for(auto& h : s->heights) {
    Serialization::read_float(stream, &h);
  }

  return stream.good();
}

bool write_spectrum(std::ostream& stream, spectrum*) {
  return stream.good();
}

} // namespace sbd
