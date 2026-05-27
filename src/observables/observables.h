//
// Created by ozdalkiran-l on 1/8/26.
//

#ifndef INC_4D_MPI_OBSERVABLES_H
#define INC_4D_MPI_OBSERVABLES_H

#include "../gauge/GaugeField.h"

double mean_plaquette_spatial_t0(const GaugeField& field, const Geometry& geo, int t0);
double mean_plaquette_temporal_t0(const GaugeField& field, const Geometry& geo, int t0);
double mean_plaquette_weighted_t0(const GaugeField& field, const Geometry& geo, int t0);
double mean_plaquette_weighted(const GaugeField& field, const Geometry& geo, int T_min, int T_max);

#endif //INC_4D_MPI_OBSERVABLES_H
