//
// Created by ozdalkiran-l on 1/13/26.
//

#ifndef INC_4D_MPI_ECMC_MPI_H
#define INC_4D_MPI_ECMC_MPI_H

#include <array>
#include <random>

#include "../gauge/GaugeField.h"
#include "../io/params.h"

// Saves the state of a local Event-Chain
struct LocalChainState {
    size_t site;
    int mu;
    int epsilon;
    SU3 R;
    bool initialized = false;

    // Budgets de Poisson persistants
    double theta_parcouru_refresh = 0.0;
    // Budgets de Poisson persistants
    double theta_sample = 0.0;
    double theta_refresh= 0.0;
    // Compteur pour le set de matrices
    size_t set_counter;
    // Compteur de lifts
    size_t event_counter = 0;
    size_t lift_counter = 0;
    size_t rev_counter = 0;
};

// Exp distributions if poisson
struct Distributions {
    std::uniform_int_distribution<int> random_dir;
    std::uniform_int_distribution<int> random_eps;
    std::exponential_distribution<double> dist_sample;
    std::exponential_distribution<double> dist_refresh;
    // Constructeur pour initialiser les taux
    Distributions(const ECMCParams& p)
        : random_dir(0, 3),
          random_eps(0, 1),
          dist_sample(1.0 / p.param_theta_sample),
          dist_refresh(1.0 / p.param_theta_refresh)
    {};
};

namespace ecmc {
void compute_list_staples(const GaugeField& field, const Geometry& geo, size_t site, int mu,
                          std::array<SU3, 6>& list_staple, std::array<double, 6>& mask_staple);
#pragma omp declare simd
inline void solve_reject_fast(double A, double B, double& gamma, double& reject, int epsilon) {
    // Utilisation de ternaires pour éviter les sauts (branches)
    B = (epsilon == -1) ? -B : B;

    double R = std::sqrt(A * A + B * B);
    double invR = 1.0 / R;
    double period = 2.0 * R;

    double discarded_number = std::floor(gamma / period);
    gamma -= discarded_number * period;

    double phi = std::atan2(-A, B);
    phi += (phi < 0.0) ? (2.0 * M_PI) : 0.0;

    double alpha;
    double p1 = R - A;

    // Le compilateur Intel transforme ce bloc en "masking" SIMD
    if (phi < (M_PI * 0.5) || phi > (M_PI * 1.5)) {
        alpha = (gamma > p1) ? (gamma - p1) * invR - 1.0 : (gamma + A) * invR;
    } else {
        alpha = gamma * invR - 1.0;
    }

    // Clamp (std::clamp est vectorisable en C++17, ou version manuelle)
    alpha = (alpha > 1.0) ? 1.0 : ((alpha < -1.0) ? -1.0 : alpha);

    double theta = phi + std::asin(alpha);

    // Normalisation 2*PI sans if/else complexes
    theta += (theta < 0.0) ? (2.0 * M_PI) : 0.0;
    theta -= (theta >= 2.0 * M_PI) ? (2.0 * M_PI) : 0.0;

    reject = theta + 2.0 * M_PI * discarded_number;
}
void solve_reject(double A, double B, double& gamma, double& reject, int epsilon);
void compute_reject_angles_fast(const GaugeField& field, size_t site, int mu,
                                const std::array<SU3, 6>& list_staple,
                                const std::array<double, 6>& mask_staple, const SU3& R, int epsilon,
                                const double& beta, std::array<double, 6>& reject_angles,
                                std::mt19937_64& rng);
size_t selectVariable_norev(const std::array<double, 3>& probas, std::mt19937_64& rng);
double compute_ds(const SU3& Pi, const SU3& R_mat);
std::pair<std::pair<size_t, int>, int> lift_improved_fast_norev(const GaugeField& field,
                                                                const Geometry& geo, size_t site,
                                                                int mu, int j, SU3& R,
                                                                std::mt19937_64& rng);
void update(GaugeField& field, size_t site, int mu, double theta, int epsilon, const SU3& R);
size_t random_site(const Geometry& geo, std::mt19937_64& rng);
void sample_persistant_norev(LocalChainState& state, Distributions& d, GaugeField& field,
                             const Geometry& geo, const ECMCParams& params, std::mt19937_64& rng);
}  // namespace ecmc

#endif  // INC_4D_MPI_ECMC_MPI_H
