#include <cstring>
#include <iostream>
#include <limits>

#include "serialization.hpp"

// In order for the serialization to work properly, we must verify that the
// floating point sizes are correct.
static_assert(
    sizeof(float) == sizeof(uint32_t),
    "error: single precision floating point size must be equal 4 bytes");
static_assert(
    sizeof(double) == sizeof(uint64_t),
    "error: double precision floating point size must be equal to 8 bytes");
static_assert(std::numeric_limits<double>::is_iec559,
              "error: floating point type is not IEEE754");

bool Serialization::read_bool(std::istream& stream, bool* value) {
  uint8_t read_bool = 0;
  read_uint8(stream, &read_bool);
  if (read_bool == 1) {
    *value = true;
  } else {
    *value = false;
  }
  return stream.good();
}

bool Serialization::write_bool(std::ostream& stream, bool value) {
  if (value) {
    write_uint8(stream, 1);
  } else {
    write_uint8(stream, 0);
  }
  return stream.good();
}

bool Serialization::read_string(std::istream& stream, std::string* value) {
  uint64_t size = 0;
  read_uint64(stream, &size);
  *value = "";
  for (size_t i = 0; i < size; ++i) {
    uint8_t read_value = 0;
    read_uint8(stream, &read_value);
    value->push_back(read_value);
  }
  return stream.good();
}

bool Serialization::write_string(std::ostream& stream,
                                 const std::string& value) {
  write_uint64(stream, value.length());
  for (size_t i = 0; i < value.length(); ++i) {
    write_uint8(stream, value[i]);
  }
  return stream.good();
}

bool Serialization::read_uint8(std::istream& stream, uint8_t* value) {
  stream.read(reinterpret_cast<char*>(value), sizeof(uint8_t));
  return stream.good();
}

bool Serialization::write_uint8(std::ostream& stream, uint8_t value) {
  stream.write(reinterpret_cast<const char*>(&value), sizeof(uint8_t));
  return stream.good();
}

bool Serialization::read_uint16(std::istream& stream, uint16_t* value) {
  // Read bytes from the stream.
  uint16_t read_value;
  stream.read(reinterpret_cast<char*>(&read_value), sizeof(uint16_t));

  // Interpret the bytes as little endian stream.
  uint8_t* bytes = reinterpret_cast<uint8_t*>(value);
  bytes[1] = read_value >> 8;
  bytes[0] = read_value >> 0;
  return stream.good();
}

bool Serialization::write_uint16(std::ostream& stream, uint16_t value) {
  // The data stream is encoded in little endian.
  uint8_t bytes[2];
  bytes[1] = value >> 8;
  bytes[0] = value >> 0;
  stream.write(reinterpret_cast<const char*>(&bytes), sizeof(uint16_t));
  return stream.good();
}

bool Serialization::read_uint32(std::istream& stream, uint32_t* value) {
  // Read bytes from the stream.
  uint32_t read_value;
  stream.read(reinterpret_cast<char*>(&read_value), sizeof(uint32_t));

  // Interpret the bytes as little endian stream.
  uint8_t* bytes = reinterpret_cast<uint8_t*>(value);
  bytes[3] = read_value >> 24;
  bytes[2] = read_value >> 16;
  bytes[1] = read_value >> 8;
  bytes[0] = read_value >> 0;
  return stream.good();
}

bool Serialization::write_uint32(std::ostream& stream, uint32_t value) {
  // The data stream is encoded in little endian.
  uint8_t bytes[4];
  bytes[3] = value >> 24;
  bytes[2] = value >> 16;
  bytes[1] = value >> 8;
  bytes[0] = value >> 0;
  stream.write(reinterpret_cast<const char*>(&bytes), sizeof(uint32_t));
  return stream.good();
}

bool Serialization::read_uint64(std::istream& stream, uint64_t* value) {
  // Read bytes from the stream.
  uint64_t read_value;
  stream.read(reinterpret_cast<char*>(&read_value), sizeof(uint64_t));

  // Interpret the bytes as little endian stream.
  uint8_t* bytes = reinterpret_cast<uint8_t*>(value);
  bytes[7] = read_value >> 56;
  bytes[6] = read_value >> 48;
  bytes[5] = read_value >> 40;
  bytes[4] = read_value >> 32;
  bytes[3] = read_value >> 24;
  bytes[2] = read_value >> 16;
  bytes[1] = read_value >> 8;
  bytes[0] = read_value >> 0;
  return stream.good();
}

bool Serialization::write_uint64(std::ostream& stream, uint64_t value) {
  // The data stream is encoded in little endian.
  uint8_t bytes[8];
  bytes[7] = value >> 56;
  bytes[6] = value >> 48;
  bytes[5] = value >> 40;
  bytes[4] = value >> 32;
  bytes[3] = value >> 24;
  bytes[2] = value >> 16;
  bytes[1] = value >> 8;
  bytes[0] = value >> 0;
  stream.write(reinterpret_cast<const char*>(&bytes), sizeof(uint64_t));
  return stream.good();
}

bool Serialization::read_int8(std::istream& stream, int8_t* value) {
  uint8_t raw_value;
  if (!read_uint8(stream, &raw_value)) {
    return false;
  }
  std::memcpy(value, &raw_value, sizeof(raw_value));
  return true;
}

bool Serialization::write_int8(std::ostream& stream, int8_t value) {
  uint8_t raw_value;
  std::memcpy(&raw_value, &value, sizeof(value));
  return write_uint8(stream, raw_value);
}

bool Serialization::read_int16(std::istream& stream, int16_t* value) {
  uint16_t raw_value;
  if (!read_uint16(stream, &raw_value)) {
    return false;
  }
  std::memcpy(value, &raw_value, sizeof(raw_value));
  return true;
}

bool Serialization::write_int16(std::ostream& stream, int16_t value) {
  uint16_t raw_value;
  std::memcpy(&raw_value, &value, sizeof(value));
  return write_uint16(stream, raw_value);
}

bool Serialization::read_int32(std::istream& stream, int32_t* value) {
  uint32_t raw_value;
  if (!read_uint32(stream, &raw_value)) {
    return false;
  }
  std::memcpy(value, &raw_value, sizeof(raw_value));
  return true;
}

bool Serialization::write_int32(std::ostream& stream, int32_t value) {
  uint32_t raw_value;
  std::memcpy(&raw_value, &value, sizeof(value));
  return write_uint32(stream, raw_value);
}

bool Serialization::read_int64(std::istream& stream, int64_t* value) {
  uint64_t raw_value;
  if (!read_uint64(stream, &raw_value)) {
    return false;
  }
  std::memcpy(value, &raw_value, sizeof(raw_value));
  return true;
}

bool Serialization::write_int64(std::ostream& stream, int64_t value) {
  uint64_t raw_value;
  std::memcpy(&raw_value, &value, sizeof(value));
  return write_uint64(stream, raw_value);
}

bool Serialization::read_float(std::istream& stream, float* value) {
  uint32_t raw_value;
  if (!read_uint32(stream, &raw_value)) {
    return false;
  }
  std::memcpy(value, &raw_value, sizeof(raw_value));
  return true;
}

bool Serialization::write_float(std::ostream& stream, float value) {
  uint32_t raw_value;
  std::memcpy(&raw_value, &value, sizeof(value));
  return write_uint32(stream, raw_value);
}

bool Serialization::read_double(std::istream& stream, double* value) {
  uint64_t raw_value;
  if (!read_uint64(stream, &raw_value)) {
    return false;
  }
  std::memcpy(value, &raw_value, sizeof(raw_value));
  return true;
}

bool Serialization::write_double(std::ostream& stream, double value) {
  uint64_t raw_value;
  std::memcpy(&raw_value, &value, sizeof(value));
  return write_uint64(stream, raw_value);
}