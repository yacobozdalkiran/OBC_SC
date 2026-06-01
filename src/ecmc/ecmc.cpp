//
// Created by ozdalkiran-l on 1/13/26.
//

#include "ecmc.h"

#include <limits>

#include "../su3/utils.h"

// Computes the list of the 6 staples around a gauge link
void ecmc::compute_list_staples(const GaugeField& field, const Geometry& geo, size_t site, int mu,
                                std::array<SU3, 6>& list_staple,
                                std::array<double, 6>& mask_staples) {
    size_t link_idx = site * 4 + mu;
    for (auto& staple : list_staple) {
        staple.setZero();
    }
    mask_staples.fill(0.0);
    // Forward staples
    for (size_t i = geo.fwd_start[link_idx]; i < geo.fwd_start[link_idx + 1]; ++i) {
        const auto& s = geo.fwd_staples_opt[i];
        list_staple[s.j_local] = s.coeff * (field.view_link_const_off(s.off0) *
                                            field.view_link_const_off(s.off1).adjoint() *
                                            field.view_link_const_off(s.off2).adjoint());
        mask_staples[s.j_local] = 1.0;
    }

    // Backward staples
    for (size_t i = geo.bwd_start[link_idx]; i < geo.bwd_start[link_idx + 1]; ++i) {
        const auto& s = geo.bwd_staples_opt[i];
        list_staple[s.j_local] += s.coeff * (field.view_link_const_off(s.off0).adjoint() *
                                             field.view_link_const_off(s.off1).adjoint() *
                                             field.view_link_const_off(s.off2));
        mask_staples[s.j_local] = 1.0;
    }
}

void ecmc::compute_reject_angles_fast(const GaugeField& field, size_t site, int mu,
                                      const std::array<SU3, 6>& list_staple,
                                      const std::array<double, 6>& mask_staple, const SU3& R,
                                      int epsilon, const double& beta,
                                      std::array<double, 6>& reject_angles, std::mt19937_64& rng) {
    static std::uniform_real_distribution<double> unif01_g(0.0, 1.0);
    const double beta_red = -(beta / 3.0);
    const SU3 T = R.adjoint() * field.view_link_const(site, mu);
    const double invalid_reject = std::numeric_limits<double>::max();
    // 1. Pré-génération des gamma
    double gammas[6];
    for (int i = 0; i < 6; ++i) {
        gammas[i] = -std::log(1.0 - unif01_g(rng));
    }

// 2. Boucle de calcul parallèle (SIMD)
// Avec Intel oneAPI, ce pragma force le compilateur à utiliser SVML pour log/atan2/asin
#pragma omp simd
    for (int i = 0; i < 6; i++) {
        // Eigen peut être utilisé à l'intérieur de omp simd si les expressions sont simples
        // Sinon, on accède directement aux données pour garantir la vectorisation

        // Calcul ligne 0
        std::complex<double> m00 = T(0, 0) * list_staple[i](0, 0) + T(0, 1) * list_staple[i](1, 0) +
                                   T(0, 2) * list_staple[i](2, 0);
        std::complex<double> m01 = T(0, 0) * list_staple[i](0, 1) + T(0, 1) * list_staple[i](1, 1) +
                                   T(0, 2) * list_staple[i](2, 1);
        std::complex<double> m02 = T(0, 0) * list_staple[i](0, 2) + T(0, 1) * list_staple[i](1, 2) +
                                   T(0, 2) * list_staple[i](2, 2);

        // Calcul ligne 1
        std::complex<double> m10 = T(1, 0) * list_staple[i](0, 0) + T(1, 1) * list_staple[i](1, 0) +
                                   T(1, 2) * list_staple[i](2, 0);
        std::complex<double> m11 = T(1, 0) * list_staple[i](0, 1) + T(1, 1) * list_staple[i](1, 1) +
                                   T(1, 2) * list_staple[i](2, 1);
        std::complex<double> m12 = T(1, 0) * list_staple[i](0, 2) + T(1, 1) * list_staple[i](1, 2) +
                                   T(1, 2) * list_staple[i](2, 2);

        std::complex<double> P00 = m00 * R(0, 0) + m01 * R(1, 0) + m02 * R(2, 0);
        std::complex<double> P11 = m10 * R(0, 1) + m11 * R(1, 1) + m12 * R(2, 1);

        const double valid = mask_staple[i];
        const double invalid = 1.0 - valid;
        double A = (P00.real() + P11.real()) * beta_red * valid;
        double B = (P11.imag() - P00.imag()) * beta_red * valid + invalid;
        // Appel de la version inline vectorisée
        double reject = 0.0;
        solve_reject_fast(A, B, gammas[i], reject, epsilon);
        reject_angles[i] = valid * reject + invalid * invalid_reject;
    }
}

size_t ecmc::selectVariable_norev(const std::array<double, 3>& probas, std::mt19937_64& rng) {
    static std::uniform_real_distribution<double> unif01(0.0, 1.0);
    double r = unif01(rng);

    if (r < probas[0]) return 0;
    if (r < probas[0] + probas[1]) return 1;
    return 2;
}

// Optimised computation of ImTr(lambda_3*R_mat.adjoint()*Pi*R_mat)
double ecmc::compute_ds(const SU3& Pi, const SU3& R_mat) {
    // Calcule Im( (R.adj * Pi * R)_00 - (R.adj * Pi * R)_11 )
    // On ne calcule que les colonnes 0 et 1 de (Pi * R)
    // Puis le produit scalaire avec les lignes de R.adjoint
    SU3 M = Pi * R_mat;
    Complex res = 0;
    for (int k = 0; k < 3; ++k) {
        res += std::conj(R_mat(k, 0)) * M(k, 0);
        res -= std::conj(R_mat(k, 1)) * M(k, 1);
    }
    return res.imag();
};

std::pair<std::pair<size_t, int>, int> ecmc::lift_improved_fast_norev(const GaugeField& field,
                                                                      const Geometry& geo,
                                                                      size_t site, int mu, int j,
                                                                      SU3& R,
                                                                      std::mt19937_64& rng) {
    // Choose a link with same probas, no reversibility
    std::array<double, 3> probas = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    size_t index_lift = selectVariable_norev(probas, rng);
    int epsilon = 1;
    if (j % 2 == 0) {
        // Forward plaquette
        if (index_lift == 0) {
            // Change R
            SU3 R2 = field.view_link_const(site, mu).adjoint() * R;
            R = R2;
            // Change epsilon
            epsilon = -1;
        }
        if (index_lift == 1) {
            auto coord_u4 = geo.get_link_staple(site, mu, j, 2);
            SU3 R3 = field.view_link_const(coord_u4.first, coord_u4.second).adjoint() * R;
            R = R3;
            // No change for epsilon
        }
        // If index_lift==2, no change for epsilon or R
    } else {
        // Backward plaquette
        auto coord_u7 = geo.get_link_staple(site, mu, j, 2);
        if (index_lift == 0) {
            auto coord_u6 = geo.get_link_staple(site, mu, j, 1);
            SU3 R5 = field.view_link_const(coord_u6.first, coord_u6.second).adjoint() *
                     field.view_link_const(coord_u7.first, coord_u7.second) * R;
            R = R5;
            // No change for epsilon
        }
        if (index_lift == 1) {
            SU3 R6 = field.view_link_const(coord_u7.first, coord_u7.second) * R;
            R = R6;
            // No change for epsilon
        }
        if (index_lift == 2) {
            SU3 R7 = field.view_link_const(coord_u7.first, coord_u7.second) * R;
            R = R7;
            epsilon = -1;
        }
    }
    return std::make_pair(geo.get_link_staple(site, mu, j, index_lift), epsilon);
}
// Updates the gauge field with XY embedding
void ecmc::update(GaugeField& field, size_t site, int mu, double theta, int epsilon, const SU3& R) {
    const SU3& Uold = field.view_link_const(site, mu);
    SU3 T = R.adjoint() * Uold;
    // 3. Application de el_3 de manière "sparse"
    // el_3 = diag(exp(i*xi), exp(-i*xi), 1)
    double xi = epsilon * theta;
    double cxi, sxi;
    sincos(xi, &sxi, &cxi);  // Utilise sincos si dispo (Linux/GNU), sinon cos et sin
    std::complex<double> PhasePlus(cxi, sxi);
    std::complex<double> PhaseMoins(cxi, -sxi);
    T.row(0) *= PhasePlus;
    T.row(1) *= PhaseMoins;
    field.view_link(site, mu) = R * T;
    field.projection_su3(site, mu);
}

void ecmc::sample_persistant_norev(LocalChainState& state, Distributions& d, GaugeField& field,
                                   const Geometry& geo, const ECMCParams& params,
                                   std::mt19937_64& rng) {
    // Constantes et Distributions
    const double beta = params.beta;
    const bool poisson = params.poisson;

    static std::uniform_int_distribution<int> random_spatial(0, geo.L - 1);
    static std::uniform_int_distribution<int> random_temp(0, geo.T - 1);

    // Initialisation de l'état de la chaîne si nécessaire
    if (!state.initialized) {
        int x_site = random_spatial(rng);
        int y_site = random_spatial(rng);
        int z_site = random_spatial(rng);
        int t_site = random_temp(rng);
        int mu_link = d.random_dir(rng);
        while (t_site == geo.T - 1 && mu_link == 3) {
            x_site = random_spatial(rng);
            y_site = random_spatial(rng);
            z_site = random_spatial(rng);
            t_site = random_temp(rng);
            mu_link = d.random_dir(rng);
        }
        state.site = geo.index(x_site, y_site, z_site, t_site);
        state.mu = mu_link;
        state.epsilon = 2 * d.random_eps(rng) - 1;
        state.R = random_su3(rng);
        state.theta_refresh_site =
            poisson ? d.dist_refresh_site(rng) : params.param_theta_refresh_site;
        state.theta_refresh_R = poisson ? d.dist_refresh_R(rng) : params.param_theta_refresh_R;
        state.theta_parcouru_refresh_site = 0.0;
        state.theta_parcouru_refresh_R = 0.0;
        state.set_counter = 0;
        state.event_counter = 0;
        state.initialized = true;
    }

    // Initalisation de l'état de la chaîne (persistant)
    size_t site_current = state.site;
    int mu_current = state.mu;
    int epsilon_current = state.epsilon;
    SU3 R = state.R;
    size_t set_counter = state.set_counter;
    size_t event_counter = state.event_counter;
    size_t lift_counter = state.lift_counter;
    size_t rev_counter = state.rev_counter;

    // Budget d'angle
    double theta_sample = poisson ? d.dist_sample(rng) : params.param_theta_sample;
    double theta_refresh_site = state.theta_refresh_site;
    double theta_refresh_R = state.theta_refresh_R;
    double theta_parcouru_sample = 0.0;
    double theta_parcouru_refresh_site = state.theta_parcouru_refresh_site;
    double theta_parcouru_refresh_R = state.theta_parcouru_refresh_R;

    // Buffers de travail
    std::array<double, 6> reject_angles;
    std::array<SU3, 6> list_staple;
    std::array<double, 6> mask_staple;

    while (true) {
        compute_list_staples(field, geo, site_current, mu_current, list_staple, mask_staple);
        compute_reject_angles_fast(field, site_current, mu_current, list_staple, mask_staple, R,
                                   epsilon_current, beta, reject_angles, rng);

        int j = 0;
        double theta_reject = reject_angles[0];
        for (int k = 1; k < 6; ++k) {
            if (reject_angles[k] < theta_reject) {
                theta_reject = reject_angles[k];
                j = k;
            }
        }

        // Distances aux frontières
        double dist_to_sample = theta_sample - theta_parcouru_sample;
        double dist_to_refresh_site = theta_refresh_site - theta_parcouru_refresh_site;
        double dist_to_refresh_R = theta_refresh_R - theta_parcouru_refresh_R;

        // Premier événement
        double theta_step =
            std::min({theta_reject, dist_to_sample, dist_to_refresh_site, dist_to_refresh_R});

        if (theta_step == dist_to_sample) {
            // --- EVENT: SAMPLE ---
            update(field, site_current, mu_current, dist_to_sample, epsilon_current, R);
            event_counter++;
            // --- SAUVEGARDE DE L'ÉTAT AVANT LE RETOUR ---
            state.site = site_current;
            state.mu = mu_current;
            state.epsilon = epsilon_current;
            state.R = R;
            state.theta_parcouru_refresh_site = theta_parcouru_refresh_site + dist_to_sample;
            state.theta_parcouru_refresh_R = theta_parcouru_refresh_R + dist_to_sample;
            state.theta_sample = theta_sample;
            state.theta_refresh_site = theta_refresh_site;
            state.theta_refresh_R = theta_refresh_R;
            state.set_counter = set_counter;
            state.event_counter = event_counter;
            state.lift_counter = lift_counter;
            state.rev_counter = rev_counter;
            return;
        } else if (theta_step == dist_to_refresh_site) {
            // --- EVENT: REFRESH SITE ---
            update(field, site_current, mu_current, dist_to_refresh_site, epsilon_current, R);
            event_counter++;

            theta_parcouru_sample += dist_to_refresh_site;
            theta_parcouru_refresh_R += dist_to_refresh_site;
            theta_parcouru_refresh_site = 0.0;
            if (poisson) theta_refresh_site = d.dist_refresh_site(rng);

            // Safety OBC : no refresh to a forbidden link
            int x_site = random_spatial(rng);
            int y_site = random_spatial(rng);
            int z_site = random_spatial(rng);
            int t_site = random_temp(rng);
            int mu_link = d.random_dir(rng);
            while (t_site == geo.T - 1 && mu_link == 3) {
                x_site = random_spatial(rng);
                y_site = random_spatial(rng);
                z_site = random_spatial(rng);
                t_site = random_temp(rng);
                mu_link = d.random_dir(rng);
            }

            site_current = geo.index(x_site, y_site, z_site, t_site);
            mu_current = mu_link;
            epsilon_current = 2 * d.random_eps(rng) - 1;
        } else if (theta_step == dist_to_refresh_R) {
            // --- EVENT: REFRESH R ---
            update(field, site_current, mu_current, dist_to_refresh_R, epsilon_current, R);
            event_counter++;

            theta_parcouru_sample += dist_to_refresh_R;
            theta_parcouru_refresh_site += dist_to_refresh_R;
            theta_parcouru_refresh_R = 0.0;
            if (poisson) theta_refresh_R = d.dist_refresh_R(rng);

            R = random_su3(rng);
        } else {
            // --- EVENT: LIFT ---
            update(field, site_current, mu_current, theta_reject, epsilon_current, R);
            event_counter++;

            theta_parcouru_sample += theta_reject;
            theta_parcouru_refresh_site += theta_reject;
            theta_parcouru_refresh_R += theta_reject;

            auto l = lift_improved_fast_norev(field, geo, site_current, mu_current, j, R, rng);
            // On lifte
            set_counter++;
            lift_counter++;
            rev_counter +=
                (l.first.first == site_current and l.first.second == mu_current and l.second == -1)
                    ? 1
                    : 0;
            site_current = l.first.first;
            mu_current = l.first.second;
            epsilon_current = epsilon_current * l.second;
            proj_su3(R);
        }
    }
}
