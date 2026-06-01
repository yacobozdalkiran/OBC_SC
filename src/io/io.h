#ifndef ECMC_MPI_IO_H
#define ECMC_MPI_IO_H

#include <random>
#include <vector>

#include "../ecmc/ecmc.h"
#include "params.h"

namespace io {
// Output
void save_plaquette(const std::vector<double>& data, const std::string& filename,
                    const std::string& dirpath, int precision);
void save_event_nb(const std::vector<size_t>& event_nb, const std::string& filename,
                   const std::string& dirpath);
void save_event_nb(const std::vector<size_t>& event_nb, const std::vector<size_t>& lift_nb,
                   const std::vector<size_t> rev_nb, const std::vector<double> lambda,
                   const std::string& filename, const std::string& dirpath);
void save_seed(std::mt19937_64& rng, const std::string& filename, const std::string& dirpath);
void save_seed(std::vector<std::mt19937_64>& rng, const std::string& filename,
               const std::string& dirpath);
void save_params(const RunParamsHbCB& rp, const std::string& filename, const std::string& dirpath);
void save_params(const RunParamsECMC& rp, const std::string& filename, const std::string& dirpath);
void add_shift(int shift, const std::string& filename, const std::string& dirpath);
void add_finished(const std::string& filename, const std::string& dirpath);
void save_state(const LocalChainState& state, const std::string& filename,
                const std::string& dirpath);
// Input
std::string trim(const std::string& s);
void load_params(const std::string& filename, RunParamsECMC& rp);
void load_params(const std::string& filename, RunParamsHbCB& rp);
void load_state(LocalChainState& state, const std::string& filename, const std::string& dirpath);

// Utilitaries
std::string format_double(double val, int precision);

// Load parameters from input file
bool read_params(RunParamsHbCB& params, const std::string& input);
bool read_params(RunParamsECMC& params, const std::string& input);
}  // namespace io

// Printing
void print_parameters(const RunParamsHbCB& rp);
void print_parameters(const RunParamsECMC& rp);
void print_time(long elapsed);

#endif  // ECMC_MPI_IO_H
