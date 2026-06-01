#include <arpa/inet.h>  // Pour htonl/htons (Linux/macOS)

#include <cstring>
#include <string>

#include "../gauge/GaugeField.h"
// Endian swap for doubles (Standard ILDG)
inline void swap_endian_64(void* ptr) {
    uint64_t* val = reinterpret_cast<uint64_t*>(ptr);
    *val = __builtin_bswap64(*val);
}

std::string generate_ildg_xml(int lx, int ly, int lz, int lt, int precision = 64);
void save_ildg(const GaugeField& field, const Geometry& geo, const std::string& run_name, const std::string& run_dir);
void read_ildg(GaugeField& field, const Geometry& geo, const std::string& run_name, const std::string& run_dir);
