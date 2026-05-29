//
// Created by ozdalkiran-l on 1/8/26.
//

#include "GaugeField.h"

#include "../su3/utils.h"

// Initialises the gauge field with random SU3 matrices
void GaugeField::hot_start(const Geometry& geo, std::mt19937_64& rng) {
    for (int t = 0; t < T; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; mu++) {
                        // OBC condition
                        if (!(t == T - 1 && mu == 3)) view_link(site, mu) = random_su3(rng);
                    }
                }
            }
        }
    }
}

// Initialises the gauge field with identity matrices
void GaugeField::cold_start(const Geometry& geo) {
    for (int t = 0; t < T; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; mu++) {
                        // OBC condition
                        if (!(t == T - 1 && mu == 3)) view_link(site, mu) = SU3::Identity();
                    }
                }
            }
        }
    }
}

// Projects a link on SU3 using Gramm-Schmidt
void GaugeField::projection_su3(size_t site, int mu) {
    auto U = view_link(site, mu);

    SU3 temp = U;

    Eigen::Vector3cd c0 = temp.col(0);
    c0.normalize();

    Eigen::Vector3cd c1 = temp.col(1);
    c1 -= c0 * c0.dot(c1);
    c1.normalize();

    Eigen::Vector3cd c2;
    c2(0) = std::conj(c0(1) * c1(2) - c0(2) * c1(1));
    c2(1) = std::conj(c0(2) * c1(0) - c0(0) * c1(2));
    c2(2) = std::conj(c0(0) * c1(1) - c0(1) * c1(0));

    temp.col(0) = c0;
    temp.col(1) = c1;
    temp.col(2) = c2;

    U = temp;
}

// Projects the whole field on SU3
void GaugeField::project_field_su3() {
    for (size_t site = 0; site < V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            projection_su3(site, mu);
        }
    }
}

// Computes the sum of all staples of a site with correct OBC coeffs using optimized branchless logic
void GaugeField::compute_staple(const Geometry& geo, size_t site, int mu, SU3& staple) const {
    staple.setZero();
    size_t link_idx = site * 4 + mu;

    // Forward staples
    for (size_t i = geo.fwd_start[link_idx]; i < geo.fwd_start[link_idx + 1]; ++i) {
        const auto& s = geo.fwd_staples_opt[i];
        staple += s.coeff * (view_link_const_off(s.off0) * 
                            view_link_const_off(s.off1).adjoint() * 
                            view_link_const_off(s.off2).adjoint());
    }

    // Backward staples
    for (size_t i = geo.bwd_start[link_idx]; i < geo.bwd_start[link_idx + 1]; ++i) {
        const auto& s = geo.bwd_staples_opt[i];
        staple += s.coeff * (view_link_const_off(s.off0).adjoint() * 
                            view_link_const_off(s.off1).adjoint() * 
                            view_link_const_off(s.off2));
    }
}
