//
// Created by ozdalkiran-l on 1/8/26.
//

#ifndef INC_4D_MPI_GAUGEFIELD_H
#define INC_4D_MPI_GAUGEFIELD_H

#include <Eigen/Dense>
#include <complex>
#include <random>
#include <vector>

#include "../geometry/Geometry.h"

using Complex = std::complex<double>;
using SU3 = Eigen::Matrix<std::complex<double>, 3, 3, Eigen::RowMajor>;

class GaugeField {
    int L;
    int T;
    size_t V;

public:
    std::vector<Complex> links;

    // Initializes a squared cold gauge conf with halos
    explicit GaugeField(const Geometry& geo)
        : L(geo.L),
          T(geo.T),
          V(geo.V),
          links(V * 4 * 9, Complex(0.0, 0.0)) {
        for (size_t site = 0; site < V; site++) {
            for (int mu = 0; mu < 4; mu++) {
                view_link(site, mu) = SU3::Identity();
            }
        }
        //OBC : we fix U_t(x , y, z, T-1) = 0
        for (int z=0; z<L; z++){
            for (int y=0; y<L; y++){
                for (int x=0; x<L; x++){
                    size_t site_tT = geo.index(x,y,z,T-1);
                    view_link(site_tT, 3) = SU3::Zero();
                }
            }
        }
    }

    void hot_start(const Geometry& geo, std::mt19937_64& rng);
    void cold_start(const Geometry& geo);

    // Non const mapping of links to SU3 matrices
    Eigen::Map<SU3> view_link(size_t site, int mu) {
        return Eigen::Map<SU3>(&links[(site * 4 + mu) * 9]);
    }

    // Const mapping of links to SU3 matrices
    [[nodiscard]] Eigen::Map<const SU3> view_link_const(size_t site, int mu) const {
        return Eigen::Map<const SU3>(&links[(site * 4 + mu) * 9]);
    }

    void projection_su3(size_t site, int mu);
    void project_field_su3();
    void compute_staple(const Geometry &geo, size_t site, int mu, SU3 &staple) const; 
};

#endif  // INC_4D_MPI_GAUGEFIELD_H
