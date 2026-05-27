#include "Geometry.h"

#include <cstdint>
#include <print>

Geometry::Geometry(int T_, int L_) {
    T = T_;
    L = L_;
    V = T * L * L * L;
    neighbors.resize(V * 8, SIZE_MAX);
    for (int t = 0; t < T; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
                    size_t site_idx = index(x, y, z, t);

                    // Direction X (mu=0)
                    neighbors[index_neigh(site_idx, 0, up)] = index((x + 1) % L, y, z, t);
                    neighbors[index_neigh(site_idx, 0, down)] = index((x - 1 + L) % L, y, z, t);
                    // Direction Y (mu=1)
                    neighbors[index_neigh(site_idx, 1, up)] = index(x, (y + 1) % L, z, t);
                    neighbors[index_neigh(site_idx, 1, down)] = index(x, (y - 1 + L) % L, z, t);
                    // Direction Z (mu=2)
                    neighbors[index_neigh(site_idx, 2, up)] = index(x, y, (z + 1) % L, t);
                    neighbors[index_neigh(site_idx, 2, down)] = index(x, y, (z - 1 + L) % L, t);
                    // Direction T (mu=3) with OBC
                    if (t < T - 1) neighbors[index_neigh(site_idx, 3, up)] = index(x, y, z, t + 1);
                    if (t > 0) neighbors[index_neigh(site_idx, 3, down)] = index(x, y, z, t - 1);
                }
            }
        }
    }


    links_staples.resize(V * 4 * 6 * 3, std::make_pair(SIZE_MAX, -1));
    // 4 links per site, 6 staples per link, 3 links per staple
    for (int t = 0; t < T; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
                    size_t site = index(x, y, z, t);  // x
                    for (int mu = 0; mu < 4; mu++) {
                        int j = 0;
                        for (int nu = 0; nu < 4; nu++) {
                            if (nu == mu) continue;
                            if (mu == 3 and t == T - 1) continue;

                            size_t xmu = get_neigh(site, mu, up);     // x+mu
                            size_t xnu = get_neigh(site, nu, up);     // x+nu
                            size_t xmunu = get_neigh(xmu, nu, down);  // x+mu-nu
                            size_t xmnu = get_neigh(site, nu, down);  // x-nu
                            // We fill only the valid staples (ie not containing a link fixed to Id)
                            if (not(t == T - 1 and mu == 3) and not(t == T - 1 and nu == 3)) {
                                links_staples[index_staples(site, mu, j, 0)] = {xmu, nu};
                                links_staples[index_staples(site, mu, j, 1)] = {xnu, mu};
                                links_staples[index_staples(site, mu, j, 2)] = {site, nu};
                            }

                            if (not(t == 0 and nu == 0) and not(t == T - 1 and mu == 3)) {
                                links_staples[index_staples(site, mu, j + 1, 0)] = {xmunu, nu};
                                links_staples[index_staples(site, mu, j + 1, 1)] = {xmnu, mu};
                                links_staples[index_staples(site, mu, j + 1, 2)] = {xmnu, nu};
                            }
                            j += 2;
                        }
                    }
                }
            }
        }
    }


    staple_valid.resize(V * 4 * 6, false);
    // 4 links per site, 6 staples per link
    for (int t = 0; t < T; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
                    size_t site = index(x, y, z, t);  // x
                    for (int mu = 0; mu < 4; mu++) {
                        int j = 0;
                        for (int nu = 0; nu < 4; nu++) {
                            if (nu == mu) continue;
                            if (mu == 3 and t == T - 1) continue;
                            // We mark the valid staples (ie not containing a link fixed to Id)
                            if (not(t == T - 1 and mu == 3) and not(t == T - 1 and nu == 3)) {
                                staple_valid[index_staple_valid(site, mu, j)] = true;
                            }

                            if (not(t == 0 and nu == 0) and not(t == T - 1 and mu == 3)) {
                                staple_valid[index_staple_valid(site, mu, j + 1)] = true;
                            }
                            j += 2;
                        }
                    }
                }
            }
        }
    }


    staple_coeff.resize(V * 4 * 6, 0);
    // 4 links per site, 6 staples per link
    for (int t = 0; t < T; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
                    size_t site = index(x, y, z, t);  // x
                    for (int mu = 0; mu < 4; mu++) {
                        int j = 0;
                        for (int nu = 0; nu < 4; nu++) {
                            if (nu == mu) continue;
                            if (mu == 3 and t == T - 1) continue;
                            // We mark the valid staples (ie not containing a link fixed to Id)
                            // Forward staple
                            if (not(t == T - 1 and mu == 3) and not(t == T - 1 and nu == 3)) {
                                //Spatial plaquette of a border spatial link
                                if ((t == T - 1 and mu != 3 and nu !=3) or (t==0 and mu !=3 and nu != 3))
                                    staple_coeff[index_staple_valid(site, mu, j)] = 0.5;
                                else
                                    staple_coeff[index_staple_valid(site, mu, j)] = 1.0;
                            }
                            // Backward staple
                            if (not(t == 0 and nu == 0) and not(t == T - 1 and mu == 3)) {
                                //Spatial plaquette of a border spatial link
                                if ((t == T - 1 and mu != 3 and nu !=3) or (t==0 and mu !=3 and nu != 3))
                                    staple_coeff[index_staple_valid(site, mu, j+1)] = 0.5;
                                else
                                    staple_coeff[index_staple_valid(site, mu, j+1)] = 1.0;
                            }
                            j += 2;
                        }
                    }
                }
            }
        }
    }

}
