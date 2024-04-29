#include "gerstner.h"
#include "math.h"
#include <stdlib.h>

double rand_double(double min, double max) {
    return min + ((double)rand() / (double)RAND_MAX) * (max - min);
}

void apply_water_simulation(Mesh* mesh, 
        GerstnerWave* waves_x, int num_waves_x, 
        GerstnerWave* waves_z, int num_waves_z, 
        double t) {
    // Seed the random number generator
    srand((unsigned int)t);

    // Accumulate a slowly varying random phase shift
    //want to get away from this
    // static double random_phase_shift = 0.0;
    // random_phase_shift += 0.1 * rand_double(-1.0, 1.0);

    for (int i = 0; i < mesh->num_vertices; i++) {
        Vertex* v = &mesh->vertices[i];

        // Compute the height of the water surface at the vertex position
        double height = 0.0;

        // Add waves in the x direction
        for (int j = 0; j < num_waves_x; j++) {
            GerstnerWave* w = &waves_x[j];
            double wi = 2.0 * M_PI / w->wavelength;
            double phase = wi * (w->speed * t );
            height += w->amplitude * (cos(wi * (v->position.x * w->direction.x + v->position.z * w->direction.z) + phase) +
                                     w->steepness * (v->position.x * w->direction.x + v->position.z * w->direction.z) * sin(wi * (v->position.x * w->direction.x + v->position.z * w->direction.z) + phase));
        }

        // Add waves in the z direction
        for (int j = 0; j < num_waves_z; j++) {
            GerstnerWave* w = &waves_z[j];
            double wi = 2.0 * M_PI / w->wavelength;
            double phase = wi * (w->speed * t );
            height += w->amplitude * (cos(wi * (v->position.x * w->direction.x + v->position.z * w->direction.z) + phase) +
                                     w->steepness * (v->position.x * w->direction.x + v->position.z * w->direction.z) * sin(wi * (v->position.x * w->direction.x + v->position.z * w->direction.z) + phase));
        }

        // Compute the partial derivatives of the water surface
        double dHdx = 0.0, dHdy = 0.0, dHdt = 0.0;

        // Compute the partial derivatives for the x-direction waves
        for (int j = 0; j < num_waves_x; j++) {
            GerstnerWave* w = &waves_x[j];
            double wi = 2.0 * M_PI / w->wavelength;
            double phase = wi * (w->speed * t );
            dHdx += -wi * w->direction.x * (w->amplitude * sin(wi * (v->position.x * w->direction.x + v->position.z * w->direction.z) + phase) +
                                           w->steepness * (v->position.x * w->direction.x + v->position.z * w->direction.z) * cos(wi * (v->position.x * w->direction.x + v->position.z * w->direction.z) + phase));
        }

        // Compute the partial derivatives for the z-direction waves
        for (int j = 0; j < num_waves_z; j++) {
            GerstnerWave* w = &waves_z[j];
            double wi = 2.0 * M_PI / w->wavelength;
            double phase = wi * (w->speed * t );
            dHdy += -wi * w->direction.z * (w->amplitude * sin(wi * (v->position.x * w->direction.x + v->position.z * w->direction.z) + phase) +
                                           w->steepness * (v->position.x * w->direction.x + v->position.z * w->direction.z) * cos(wi * (v->position.x * w->direction.x + v->position.z * w->direction.z) + phase));
        }

        dHdt = 0.0;
        for (int j = 0; j < num_waves_x; j++) {
            GerstnerWave* w = &waves_x[j];
            double wi = 2.0 * M_PI / w->wavelength;
            double phase = wi * (w->speed * t );
            dHdt += -wi * w->speed * w->amplitude * sin(wi * (v->position.x * w->direction.x + v->position.z * w->direction.z) + phase);
        }
        for (int j = 0; j < num_waves_z; j++) {
            GerstnerWave* w = &waves_z[j];
            double wi = 2.0 * M_PI / w->wavelength;
            double phase = wi * (w->speed * t );
            dHdt += -wi * w->speed * w->amplitude * sin(wi * (v->position.x * w->direction.x + v->position.z * w->direction.z) + phase);
        }

        // Compute the binormal, tangent, and normal vectors
        Vector3 binormal = {1.0, dHdx, 0.0};
        Vector3 tangent = {0.0, dHdy, 1.0};
        Vector3 normal = vec3_normalized(vec3_cross_prod(binormal, tangent));

        // Update the vertex position along the normal vector
        v->position = vec3_add(v->position_static, vec3_scale(v->normal_static, height));
    }

    // Update the mesh normals
    compute_plane_normals(mesh);
    // compute_face_normals(mesh);

}

void reset_water_simulation(Mesh* mesh) {
    for (int i = 0; i < mesh->num_vertices; i++) {
        Vertex* v = &mesh->vertices[i];
        v->position = v->position_static;
    }
    compute_face_normals(mesh);
}