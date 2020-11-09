#ifndef UTILS_SERIALIZATION_HPP
#define UTILS_SERIALIZATION_HPP

#include <iostream>
#include <vector>

// This namespace contains necessary functions to serialize commonly used types
// into a binary stream using the little endian byte order.
namespace Serialization {

// Write/read a boolean as a single byte to/from the stream.
bool read_bool(std::istream& stream, bool* value);
bool write_bool(std::ostream& stream, bool value);

// Write/read a string to/from the stream.
bool read_string(std::istream& stream, std::string* value);
bool write_string(std::ostream& stream, const std::string& value);

// Write/read a single byte to/from the stream.
bool read_uint8(std::istream& stream, uint8_t* value);
bool write_uint8(std::ostream& stream, uint8_t value);

// Write/read an uint16 to/from the stream.
bool read_uint16(std::istream& stream, uint16_t* value);
bool write_uint16(std::ostream& stream, uint16_t value);

// Write/read an uint32 to/from the stream.
bool read_uint32(std::istream& stream, uint32_t* value);
bool write_uint32(std::ostream& stream, uint32_t value);

// Write/read an uint64 to/from the stream.
bool read_uint64(std::istream& stream, uint64_t* value);
bool write_uint64(std::ostream& stream, uint64_t value);

// Write/read a single byte to/from the stream.
bool read_int8(std::istream& stream, int8_t* value);
bool write_int8(std::ostream& stream, int8_t value);

// Write/read an int16 to/from the stream.
bool read_int16(std::istream& stream, int16_t* value);
bool write_int16(std::ostream& stream, int16_t value);

// Write/read an int32 to/from the stream.
bool read_int32(std::istream& stream, int32_t* value);
bool write_int32(std::ostream& stream, int32_t value);

// Write/read an int64 to/from the stream.
bool read_int64(std::istream& stream, int64_t* value);
bool write_int64(std::ostream& stream, int64_t value);

// Write/read a single precision floating point value to/from the stream.
bool read_float(std::istream& stream, float* value);
bool write_float(std::ostream& stream, float value);

// Write/read a double precision floating point value to/from the stream.
bool read_double(std::istream& stream, double* value);
bool write_double(std::ostream& stream, double value);

// Write/read a generic vector to/from the stream.
template <class T, class Func>
bool read_vector(std::istream& stream, std::vector<T>* value, Func f) {
  uint64_t num_elements = 0;
  read_uint64(stream, &num_elements);
  *value = std::vector<T>(num_elements);
  for (size_t i = 0; i < num_elements; ++i) {
    f(stream, &(*value)[i]);
  }
  return stream.good();
}

template <class T, class Func>
bool write_vector(std::ostream& stream, const std::vector<T>& value, Func f) {
  uint64_t num_elements = value.size();
  write_uint64(stream, num_elements);
  for (size_t i = 0; i < num_elements; ++i) {
    f(stream, value[i]);
  }
  return stream.good();
}




}  // namespace Serialization

#endif /* UTILS_SERIALIZATION_HPP */