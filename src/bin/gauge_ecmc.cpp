#include <omp.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "../ecmc/ecmc.h"
#include "../gauge/GaugeField.h"
#include "../io/ildg.h"
#include "../io/io.h"
#include "../observables/observables.h"

namespace fs = std::filesystem;

void generate_ecmc_cb(const RunParamsECMC& rp, bool existing) {
    //========================Objects initialization====================
    // MPI
    // Lattice creation + RNG
    Geometry geo(rp.T, rp.L);
    GaugeField field(geo);
    std::mt19937_64 rng(rp.seed);
    if (!rp.cold_start) {
        field.hot_start(geo, rng);
    }

    // Chain state
    LocalChainState state{};
    Distributions d(rp.ecmc_params);

    if (existing) {
        read_ildg(field, geo, rp.run_name, rp.run_dir);
        io::load_state(state, rp.run_name, rp.run_dir);

        fs::path state_path = fs::path(rp.run_dir) / rp.run_name / (rp.run_name + "_seed") /
                              (rp.run_name + "_seed.txt");
        std::ifstream ifs(state_path);
        if (ifs.is_open()) {
            ifs >> rng;
        } else {
            // Optionnel : Alerte si on attend un checkpoint mais qu'il manque un morceau
            std::cerr << "Seed file missing, using initial seed." << std::endl;
        }
    }

    // Params ECMC
    ECMCParams ep = rp.ecmc_params;

    // Measure vectors
    std::vector<double> plaquette;
    plaquette.reserve(rp.save_each / rp.N_plaquette + 2);

    std::vector<size_t> event_nb;
    std::vector<size_t> lift_nb;
    std::vector<double> lambda;
    std::vector<size_t> rev_nb;
    lift_nb.reserve(rp.save_each / rp.N_plaquette + 2);
    event_nb.reserve(rp.save_each / rp.N_plaquette + 2);
    lambda.reserve(rp.save_each / rp.N_plaquette + 2);
    rev_nb.reserve(rp.save_each / rp.N_plaquette + 2);

    io::save_params(rp, rp.run_name, rp.run_dir);

    // Print params
    print_parameters(rp);

    //==============================ECMC Checkboard===========================
    // Thermalisation
    // Skipped if existing (for successive jobs in slurm)

    if (!existing) {
        std::cout << "\n\n===========================================\n";
        std::cout << "Thermalisation : " << rp.N_therm << " configurations \n";
        std::cout << "===========================================\n";

        for (int i = 0; i < rp.N_therm; i++) {
            ecmc::sample_persistant_norev(state, d, field, geo, ep, rng);

            // Plaquette measure (not saved for thermalization)
            if (i % rp.N_plaquette == 0) {
                double p = mean_plaquette_weighted(field, geo, rp.T_min, rp.T_max);
                std::cout << "\n====== Configuration " << i << " ======\n";
                std::cout << "(Therm) Sample " << i / rp.N_plaquette << ", <P> = " << p << "\n";
            }
        }
    }

    // Sampling
    std::cout << "\n\n===========================================\n";
    std::cout << "Sampling : " << rp.N_samples / rp.N_plaquette << " <P> samples\n "
              << "===========================================\n";

    for (int i = 0; i < rp.N_samples; i++) {
        ecmc::sample_persistant_norev(state, d, field, geo, ep, rng);

        // Plaquette measure
        if ((i % rp.N_plaquette == 0) and (i > 0 or !existing)) {
            double p = mean_plaquette_weighted(field, geo, rp.T_min, rp.T_max);
            std::cout << "\n====== Configuration " << i << " ======\n";
            std::cout << "Sample " << i / rp.N_plaquette << ", <P> = " << p << "\n";
            plaquette.emplace_back(p);
            unsigned long local_lifts = state.lift_counter;
            unsigned long local_events = state.event_counter;
            unsigned long local_rev = state.rev_counter;
            double local_lambda = (rp.N_plaquette * state.theta_sample) / (double)local_lifts;

            lambda.emplace_back(local_lambda);
            lift_nb.emplace_back(local_lifts);
            event_nb.emplace_back(local_events);
            rev_nb.emplace_back(local_rev);

            std::cout << ">>> Average \u03bb = " << local_lambda << ", Events : " << local_events
                      << " (Reversibility : " << (double)local_rev / (double)local_lifts * 100
                      << "%)\n";
            // Event counter reinitialized
            state.event_counter = 0;
            state.lift_counter = 0;
            state.rev_counter = 0;
        }

        // Save conf/seed/chain state/obs
        if (i > 0 and i % rp.save_each == 0) {
            std::cout << "\n\n==========================================\n";
            // Write the output
            int precision = 10;
            io::save_plaquette(plaquette, rp.run_name, rp.run_dir, precision);
            io::add_shift(i, rp.run_name, rp.run_dir);
            io::save_event_nb(event_nb, lift_nb, rev_nb, lambda, rp.run_name, rp.run_dir);
            // Save conf
            save_ildg(field, geo, rp.run_name, rp.run_dir);
            // Save seeds
            io::save_seed(rng, rp.run_name, rp.run_dir);
            // Save chain state
            io::save_state(state, rp.run_name, rp.run_dir);
            std::cout << "==========================================\n";
            // Clear the measures
            plaquette.clear();
            plaquette.reserve(rp.save_each / rp.N_plaquette + 2);
            event_nb.clear();
            event_nb.reserve(rp.save_each / rp.N_plaquette + 2);
            lift_nb.clear();
            lift_nb.reserve(rp.save_each / rp.N_plaquette + 2);
            rev_nb.clear();
            rev_nb.reserve(rp.save_each / rp.N_plaquette + 2);
            lambda.clear();
            lambda.reserve(rp.save_each / rp.N_plaquette + 2);
        }
    }

    //===========================Output======================================

    // Save conf/seed/chain state/obs
    std::cout << "\n\n==========================================\n";
    // Write the output
    int precision = 10;
    io::save_plaquette(plaquette, rp.run_name, rp.run_dir, precision);
    io::add_shift(rp.N_samples, rp.run_name, rp.run_dir);
    io::add_finished(rp.run_name, rp.run_dir);
    io::save_event_nb(event_nb, lift_nb, rev_nb, lambda, rp.run_name, rp.run_dir);
    // Save conf
    save_ildg(field, geo, rp.run_name, rp.run_dir);
    // Save seeds
    io::save_seed(rng, rp.run_name, rp.run_dir);
    // Save chain state
    io::save_state(state, rp.run_name, rp.run_dir);
    std::cout << "==========================================\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.txt>" << std::endl;
        return 1;
    }

    // Charging the parameters of the run
    RunParamsECMC params;
    bool existing = io::read_params(params, argv[1]);

    // Measuring time
    auto start_time = std::chrono::high_resolution_clock::now();

    generate_ecmc_cb(params, existing);

    auto end_time = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> total_time = end_time - start_time;
    std::cout << std::fixed << std::setprecision(4);
    print_time(total_time.count());
    // End MPI
}
