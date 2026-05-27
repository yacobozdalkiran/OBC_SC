#ifndef GEOMETRYCB_H
#define GEOMETRYCB_H

#include <vector>

#include "types.h"

class Geometry {
public:
    int T;
    int L;
    size_t V;

private:
    std::vector<size_t> neighbors;
    std::vector<std::pair<size_t, int>> links_staples;
    std::vector<bool> staple_valid;

public:
    explicit Geometry(int T_, int L_);

    // Index function for links
    [[nodiscard]] size_t index(int x, int y, int z, int t) const {
        return ((static_cast<size_t>(t) * L + z) * L + y) * L + x;
    }

    // Index of a site in neighbor vector
    [[nodiscard]] static size_t index_neigh(size_t site, int mu, dir d) {
        return site * 8 + mu * 2 + d;
    }

    // Index of a link in links_staples
    [[nodiscard]] static size_t index_staples(size_t site, int mu, int i_staple, int i_link) {
        return site * 3 * 6 * 4 + mu * 3 * 6 + i_staple * 3 + i_link;
    }

    // Index of a staple in staple_valid
    [[nodiscard]] static size_t index_staple_valid(size_t site, int mu, int i_staple) {
        return site * 3 * 6 * 4 + mu * 3 * 6 + i_staple;
    }

    // Get a neighbor
    [[nodiscard]] size_t get_neigh(size_t site, int mu, dir d) const {
        return neighbors[index_neigh(site, mu, d)];
    }

    // Returns the coords of the link of index i_link of the staple of index i_staple of <site,mu>
    [[nodiscard]] std::pair<size_t, int> get_link_staple(size_t site, int mu, int i_staple,
                                                         int i_link) const {
        return links_staples[index_staples(site, mu, i_staple, i_link)];
    }
};

#endif
