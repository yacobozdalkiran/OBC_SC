#include <iostream>
#include <print>

#include "../gauge/GaugeField.h"
#include "../observables/observables.h"

int main(int argc, char* argv[]) {
    std::cout << "Configuration initialization\n";
    int T = 8;
    int L = 4;
    Geometry geo(T, L);
    GaugeField field(geo);

    int x = 2;
    int y = 3;
    int z = 3;
    int t = 6;
    int mu = 3;
    std::println("Site ({},{},{},{}) Mu = {}", x, y, z, t, mu);
    std::cout << field.view_link_const(geo.index(x, y, z, t), mu) << "\n";
    std::cout << "Hot start\n";
    int seed = 123;
    std::mt19937_64 rng(seed);
    field.hot_start(geo, rng);
    std::println("Site ({},{},{},{}) Mu = {}", x, y, z, t, mu);
    std::cout << field.view_link_const(geo.index(x, y, z, t), mu) << "\n";
    for (int i = 0; i < T; i++) {
        double pt = 0.0;
        double p = mean_plaquette_spatial_t0(field, geo, i);
        if (i < T - 1) pt = mean_plaquette_temporal_t0(field, geo, i);
        std::println("Plaquette sp at time {} = {}", i, p);
        if (i < T - 1) std::println("Plaquette temp at time {} = {}", i, pt);
    }
    std::println("Total mean plaquette = {}", mean_plaquette_weighted(field, geo, 0, T-1));
    field.cold_start(geo);
    std::println("Cold start");
    std::println("Site ({},{},{},{}) Mu = {}", x, y, z, t, mu);
    std::cout << field.view_link_const(geo.index(x, y, z, t), mu) << '\n';
    for (int i = 0; i < T; i++) {
        double pt = 0.0;
        double p = mean_plaquette_spatial_t0(field, geo, i);
        if (i < T - 1) pt = mean_plaquette_temporal_t0(field, geo, i);
        std::println("Plaquette sp at time {} = {}", i, p);
        if (i < T - 1) std::println("Plaquette temp at time {} = {}", i, pt);
    }
    std::println("Total mean plaquette = {}", mean_plaquette_weighted(field, geo, 0, T-1));
    return 0;
}
