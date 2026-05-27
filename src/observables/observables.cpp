//
// Created by ozdalkiran-l on 1/8/26.
//
#include "observables.h"

#include <iostream>

// Returns the mean spatial plaquette at a time slice t0 (average over the spatial plaquettes at
// time t0)
double mean_plaquette_spatial_t0(const GaugeField& field, const Geometry& geo, int t0) {
    double sum = 0.0;
    int count = 0;
    for (int z = 0; z < geo.L; z++) {
        for (int y = 0; y < geo.L; y++) {
            for (int x = 0; x < geo.L; x++) {
                size_t site = geo.index(x, y, z, t0);
                for (int mu = 0; mu < 3; mu++) {
                    for (int nu = mu + 1; nu < 3; nu++) {
                        auto U0 = field.view_link_const(site, mu);
                        auto U1 = field.view_link_const(geo.get_neigh(site, mu, up), nu);
                        auto U2 = field.view_link_const(geo.get_neigh(site, nu, up), mu).adjoint();
                        auto U3 = field.view_link_const(site, nu).adjoint();
                        sum += (U0 * U1 * U2 * U3).trace().real() / 3.0;
                        ++count;
                    }
                }
            }
        }
    }
    return sum / count;
}

// Returns the mean temporal plaquettes at time t0 (only possible if t0 < T-1)
double mean_plaquette_temporal_t0(const GaugeField& field, const Geometry& geo, int t0) {
    if (t0 > geo.T - 2) {
        std::cerr << "Invalid time slice t0 for temporal plaquettes\n";
        return 0.0;
    }
    int count = 0;
    double sum = 0.0;
    for (int z = 0; z < geo.L; z++) {
        for (int y = 0; y < geo.L; y++) {
            for (int x = 0; x < geo.L; x++) {
                size_t site = geo.index(x, y, z, t0);
                for (int mu = 0; mu < 3; mu++) {
                    int nu = 3;
                    auto U0 = field.view_link_const(site, mu);
                    auto U1 = field.view_link_const(geo.get_neigh(site, mu, up), nu);
                    auto U2 = field.view_link_const(geo.get_neigh(site, nu, up), mu).adjoint();
                    auto U3 = field.view_link_const(site, nu).adjoint();
                    sum += (U0 * U1 * U2 * U3).trace().real() / 3.0;
                    ++count;
                }
            }
        }
    }
    return sum / count;
}

// Returns the temporal+spatial weighted mean plaquette at time t0
double mean_plaquette_weighted_t0(const GaugeField& field, const Geometry& geo, int t0) {
    double res = 0.0;
    double count = 0.0;
    res += ((t0 == 0 or t0 == geo.T - 1) ? 0.5 : 1.0) * mean_plaquette_spatial_t0(field, geo, t0);
    count += ((t0 == 0 or t0 == geo.T - 1) ? 0.5 : 1.0);
    if (t0 < geo.T - 1) {
        res += mean_plaquette_temporal_t0(field, geo, t0);
        count += 1;
    }
    return res / count;
}

// Returns the weighted mean plaquette between T_min, T_max with T_min<=T_max
double mean_plaquette_weighted(const GaugeField& field, const Geometry& geo, int T_min, int T_max) {
    if (T_max < T_min) {
        std::cerr << "Invalid time slice T_min > T_max\n";
        return -1;
    }
    double total_weighted_sum = 0.0;
    double total_weight = 0.0;
    for (int t = T_min; t <= T_max; t++) {
        // Spatial part
        double w_s = ((t == 0 or t == geo.T - 1) ? 0.5 : 1.0);
        total_weighted_sum += w_s * mean_plaquette_spatial_t0(field, geo, t);
        total_weight += w_s;

        // Temporal part (connects t to t+1)
        if (t < geo.T - 1 and t < T_max) {
            total_weighted_sum += mean_plaquette_temporal_t0(field, geo, t);
            total_weight += 1.0;
        }
    }
    return total_weighted_sum / total_weight;
}
