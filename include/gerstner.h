#include "mesh.h"

#ifndef GERSTNER_H
#define GERSTNER_H

typedef struct {
    double wavelength;
    double amplitude;
    double speed;
    Vector3 direction;
    double steepness;
} GerstnerWave;




void apply_water_simulation(Mesh* mesh, 
        GerstnerWave* waves_x, int num_waves_x, 
        GerstnerWave* waves_z, int num_waves_z, 
        double t);

void reset_water_simulation(Mesh* mesh);

#endif