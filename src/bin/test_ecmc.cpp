#include <chrono>
#include <iostream>
#include <print>

#include "../ecmc/ecmc.h"
#include "../gauge/GaugeField.h"
#include "../observables/observables.h"

int main() {
    std::println("Test ECMC");
    int T, L, seed, n_sweeps, T_min, T_max;
    double beta;

    // 1. Gather all inputs first
    std::cout << "T : ";
    std::cin >> T;
    std::cout << "L : ";
    std::cin >> L;
    std::print("Seed : ");
    std::cin >> seed;
    std::print("Number of sweeps : ");
    std::cin >> n_sweeps;
    std::print("Beta : ");
    std::cin >> beta;
    std::print("T_min : ");
    std::cin >> T_min;
    std::print("T_max : ");
    std::cin >> T_max;
    std::println("");

    // 2. Start total computation timer
    auto total_start = std::chrono::high_resolution_clock::now();

    // 3. Initialization
    auto init_start = std::chrono::high_resolution_clock::now();
    std::print("Geometry initialisation... ");
    Geometry geo(T, L);
    std::print("Done !\n");
    std::print("Configuration initialisation... ");
    GaugeField field(geo);
    std::print("Done !\n");
    auto init_end = std::chrono::high_resolution_clock::now();
    std::println(
        "Initialization time: {} ms",
        std::chrono::duration_cast<std::chrono::milliseconds>(init_end - init_start).count());

    std::mt19937_64 rng(seed);
    double plaquette{};
    ECMCParams params{.beta = beta,
                      .N_samples = n_sweeps,
                      .param_theta_sample = 600,
                      .param_theta_refresh= 200,
                      .poisson = false,
                      .epsilon_set = 0.15};
    LocalChainState state{};
    Distributions d{params};

    // 4. Simulation Loop
    for (int sweeps = 0; sweeps < n_sweeps; sweeps++) {
        std::println("===========================");
        std::println("Configuration {}", sweeps);

        auto sweep_start = std::chrono::high_resolution_clock::now();
        ecmc::sample_persistant_norev(state, d, field, geo, params, rng);
        auto sweep_end = std::chrono::high_resolution_clock::now();
        std::println(
            "Sweep time: {} ms",
            std::chrono::duration_cast<std::chrono::milliseconds>(sweep_end - sweep_start).count());

        auto plaq_start = std::chrono::high_resolution_clock::now();
        plaquette = mean_plaquette_weighted(field, geo, T_min, T_max);
        auto plaq_end = std::chrono::high_resolution_clock::now();
        std::println("<P> = {}", plaquette);
        std::println(
            "Plaquette measurement time: {} ms",
            std::chrono::duration_cast<std::chrono::milliseconds>(plaq_end - plaq_start).count());
    }

    auto total_end = std::chrono::high_resolution_clock::now();
    std::println("===========================");
    std::println("Total computation time: {} s",
                 std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start).count());

    return 0;
}
