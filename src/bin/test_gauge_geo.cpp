#include <iostream>
#include <print>

#include "../gauge/GaugeField.h"

int main(int argc, char* argv[]) {
    std::cout <<"Configuration initialization\n";
    int T = 8;
    int L = 4;
    Geometry geo(T,L);
    GaugeField field(geo);
    int x = 2;
    int y = 3;
    int z = 3;
    int t = 7;
    int mu = 3;
    std::println("Site ({},{},{},{}) Mu = {}", x,y,z,t,mu);
    std::cout << field.view_link_const(geo.index(x,y,z,t), mu) << "\n";
    std::cout << "Hot start\n";
    int seed = 123;
    std::mt19937_64 rng(seed);
    field.hot_start(geo, rng);
    std::println("Site ({},{},{},{}) Mu = {}", x,y,z,t,mu);
    std::cout << field.view_link_const(geo.index(x,y,z,t), mu) << "\n";
    field.cold_start(geo);
    std::println("Cold start");
    std::println("Site ({},{},{},{}) Mu = {}", x,y,z,t,mu);
    std::cout << field.view_link_const(geo.index(x,y,z,t),mu) << '\n';
    return 0;
}
