//
// Created by ozdalkiran-l on 1/13/26.
//

#include "ecmc.h"

#include "../su3/utils.h"

// Computes the list of the 6 staples around a gauge link
void ecmc::compute_list_staples(const GaugeField& field, const Geometry& geo, size_t site, int mu,
                                std::array<SU3, 6>& list_staple) {
    size_t index = 0;
    size_t x = site;                        // x
    size_t xmu = geo.get_neigh(x, mu, up);  // x+mu
    for (int nu = 0; nu < 4; nu++) {
        if (nu == mu) {
            continue;
        }
        // Staple forward
        size_t xnu = geo.get_neigh(x, nu, up);  // x+nu
        const auto& U0 = field.view_link_const(xmu, nu);
        const auto& U1 = field.view_link_const(xnu, mu);
        const auto& U2 = field.view_link_const(x, nu);
        list_staple[index] = U0 * (U2 * U1).adjoint();

        // Staple backward
        size_t xmunu = geo.get_neigh(xmu, nu, down);  // x+mu-nu
        size_t xmnu = geo.get_neigh(x, nu, down);     // x-nu
        auto V0 = field.view_link_const(xmunu, nu);
        auto V1 = field.view_link_const(xmnu, mu);
        auto V2 = field.view_link_const(xmnu, nu);
        list_staple[index + 1] = (V1 * V0).adjoint() * V2;
        index += 2;
    }
}

// Generates the 6 reject angles for a link
void ecmc::compute_reject_angles(const GaugeField& field, const Geometry& geo, size_t site, int mu,
                                 const std::array<SU3, 6>& list_staple, const SU3& R, int epsilon,
                                 const double& beta, std::array<double, 6>& reject_angles,
                                 std::mt19937_64& rng) {
    static std::uniform_real_distribution<double> unif01_g(0.0, 1.0);
    SU3 T = R.adjoint() * field.view_link_const(site, mu);
    const double beta_red = -(beta / 3.0);
    for (int i = 0; i < 6; i++) {
        double gamma = -std::log(unif01_g(rng));
        auto M_row0 = T.row(0) * list_staple[i];
        auto M_row1 = T.row(1) * list_staple[i];
        // P(0,0) = M_row0 * R_col0
        std::complex<double> P00 = M_row0(0) * R(0, 0) + M_row0(1) * R(1, 0) + M_row0(2) * R(2, 0);
        // P(1,1) = M_row1 * R_col1
        std::complex<double> P11 = M_row1(0) * R(0, 1) + M_row1(1) * R(1, 1) + M_row1(2) * R(2, 1);
        double A = (P00.real() + P11.real()) * beta_red;
        double B = (-P00.imag() + P11.imag()) * beta_red;
        solve_reject_fast(A, B, gamma, reject_angles[i], epsilon);
    }
}

void ecmc::compute_reject_angles_fast(const GaugeField& field, const Geometry& geo, size_t site, int mu,
                                      const std::array<SU3, 6>& list_staple, const SU3& R,
                                      int epsilon, const double& beta,
                                      std::array<double, 6>& reject_angles, std::mt19937_64& rng) {
    static std::uniform_real_distribution<double> unif01_g(0.0, 1.0);
    const double beta_red = -(beta / 3.0);
    const SU3 T = R.adjoint() * field.view_link_const(site, mu);
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

        double A = (P00.real() + P11.real()) * beta_red;
        double B = (P11.imag() - P00.imag()) * beta_red;
        // Appel de la version inline vectorisée
        solve_reject_fast(A, B, gammas[i], reject_angles[i], epsilon);
    }
}

// Selects an index between 0 and probas.size()-1 using the tower of probability method
// static dist to avoid initialization cost
size_t ecmc::selectVariable(const std::array<double, 4>& probas, std::mt19937_64& rng) {
    static std::uniform_real_distribution<double> unif01(0.0, 1.0);
    double r = unif01(rng);

    if (r < probas[0]) return 0;
    if (r < probas[0] + probas[1]) return 1;
    if (r < probas[0] + probas[1] + probas[2]) return 2;
    return 3;
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
// Optimized version of lift_improved
std::pair<std::pair<size_t, int>, int> ecmc::lift_improved_fast(const GaugeField& field,
                                                                const Geometry& geo, size_t site,
                                                                int mu, int j, SU3& R,
                                                                std::mt19937_64& rng) {
    std::array<std::pair<size_t, int>, 4>
        links_plaquette_j;  // We add the current link to get the plaquette
    links_plaquette_j[0] = std::make_pair(site, mu);
    links_plaquette_j[1] = geo.get_link_staple(site, mu, j, 0);
    links_plaquette_j[2] = geo.get_link_staple(site, mu, j, 1);
    links_plaquette_j[3] = geo.get_link_staple(site, mu, j, 2);

    SU3 U0 = field.view_link_const(site, mu);
    SU3 U1 = field.view_link_const(links_plaquette_j[1].first, links_plaquette_j[1].second);
    SU3 U2 = field.view_link_const(links_plaquette_j[2].first, links_plaquette_j[2].second);
    SU3 U3 = field.view_link_const(links_plaquette_j[3].first, links_plaquette_j[3].second);
    std::array<double, 4> probas{};
    std::array<double, 4> abs_dS{};
    double sum = 0.0;
    std::array<int, 4> sign_dS{};
    std::array<SU3, 4> P{};

    if (j % 2 == 0) {  // Forward plaquette
        SU3 U01 = U0 * U1;
        SU3 U32 = U3 * U2;
        P[0] = U01 * U32.adjoint();
        P[1] = U1 * U32.adjoint() * U0;
        P[2] = U2 * U01.adjoint() * U3;
        P[3] = P[0].adjoint();
    } else {                // Backward plaquette
        SU3 U21 = U2 * U1;  // 1 mult + adjoint
        SU3 T = U0 * U21.adjoint();
        P[0] = T * U3;
        P[1] = U1 * U0.adjoint() * U3.adjoint() * U2;
        P[3] = U3 * T;
        P[2] = P[3].adjoint();
    }
    for (size_t i = 0; i < 4; i++) {
        probas[i] = compute_ds(P[i], R);  // Less matmul
        sign_dS[i] = dsign(probas[i]);
        probas[i] = abs(probas[i]);
        abs_dS[i] = probas[i];
        sum += probas[i];
    }

    for (size_t i = 0; i < 4; i++) {
        probas[i] /= sum;
    }

    size_t index_lift = selectVariable(probas, rng);
    int new_epsilon = -sign_dS[index_lift];
    return make_pair(links_plaquette_j[index_lift], new_epsilon);
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

// Returns a random non frozen site
size_t ecmc::random_site(const Geometry& geo, std::mt19937_64& rng) {
    static std::uniform_int_distribution<size_t> random_coord(0, geo.V - 1);
    return random_coord(rng);
}

void ecmc::sample_persistant(LocalChainState& state, Distributions& d, GaugeField& field,
                             const Geometry& geo, const ECMCParams& params,
                             std::mt19937_64& rng) {
    // Constantes et Distributions
    const double beta = params.beta;
    const bool poisson = params.poisson;

    // Initialisation de l'état de la chaîne si nécessaire
    if (!state.initialized) {
        state.site = random_site(geo, rng);
        state.mu = d.random_dir(rng);
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

    // Buffer de matrices (Optimisation : Statique pour éviter l'allocation)
    // static std::vector<SU3> set_matrices(101);
    // ecmc_set(params.epsilon_set, set_matrices, rng);

    // Buffers de travail
    std::array<double, 6> reject_angles;
    std::array<SU3, 6> list_staple;

    while (true) {
        compute_list_staples(field, geo, site_current, mu_current, list_staple);
        compute_reject_angles_fast(field, geo, site_current, mu_current, list_staple, R, epsilon_current,
                                   beta, reject_angles, rng);

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

            site_current = random_site(geo, rng);
            mu_current = d.random_dir(rng);
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

            auto l = lift_improved_fast(field, geo, site_current, mu_current, j, R, rng);
            set_counter++;
            lift_counter++;
            rev_counter += (l.first.first == site_current and l.first.second == mu_current) ? 1 : 0;
            site_current = l.first.first;
            mu_current = l.first.second;
            epsilon_current = l.second;
        }
    }
}

void ecmc::sample_persistant_norev(LocalChainState& state, Distributions& d, GaugeField& field,
                                   const Geometry& geo, const ECMCParams& params,
                                   std::mt19937_64& rng) {
    // Constantes et Distributions
    const double beta = params.beta;
    const bool poisson = params.poisson;

    // Initialisation de l'état de la chaîne si nécessaire
    if (!state.initialized) {
        state.site = random_site(geo, rng);
        state.mu = d.random_dir(rng);
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

    while (true) {
        compute_list_staples(field, geo, site_current, mu_current, list_staple);
        compute_reject_angles_fast(field, geo, site_current, mu_current, list_staple, R, epsilon_current,
                                   beta, reject_angles, rng);

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

            site_current = random_site(geo, rng);
            mu_current = d.random_dir(rng);
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
