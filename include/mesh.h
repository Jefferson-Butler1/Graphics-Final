#ifndef MESH_H
#define MESH_H
#include "vector.h"
#include "lightmodel.h"


typedef struct {
    Vector3 position;
    Vector3 normal;
} Vertex;

typedef struct {
    // a, b, and c are indices of vertices in the mesh
    Vertex* a;
    Vertex* b;
    Vertex* c;
    Vector3 normal;
} Triangle;

typedef struct {

} cachedMesh;

typedef struct {
    int num_tris;
    Triangle* tris;
    int num_vertices;
    Vertex* vertices;

    Vector3 center;
    Vector3 scale;

    Vector3 bounding_box_max;
    Vector3 bounding_box_min;

    double transform[4][4];
    double inverse_transform[4][4];

    PhongMaterial material;
    double roughness;

} Mesh;

/**
 * @brief Loads triangles from an ascii .ply file into a given mesh struct.
 * Assumes that the mesh is triangulated and that it has precomputed vertex normals only.
 * 
 * @param mesh A pointer to the mesh to load the data into
 * @param filename The name of the file to be loaded
 */
void load_mesh_from_ply(Mesh* mesh, char* filename);

/**
 * @brief Frees all mesh vertices and faces from memory
 * 
 * @param mesh the mesh to be deleted
 */
void delete_mesh(Mesh mesh);

/**
 * @brief Translates a mesh's transform by its transform matrix.
 * 
 * @param mesh The mesh to be translated
 */
void translate_mesh(Mesh* mesh, Vector3 translation);

/**
 * @brief Scales a mesh's transform by the given scaling factors.
 * 
 * @param mesh The mesh to be scaled
 * @param scale The scale vector to use 
 */
void scale_mesh(Mesh* mesh, Vector3 scale);

/**
 * @brief Rotates a Mesh on the x axis
 * 
 * @param mesh the mesh to be rotated
 * @param degrees the degrees on the x axis to rotate the mesh by
 */

void rotate_mesh_x_degrees(Mesh* mesh, double degrees);
/**
 * @brief Rotates a Mesh on the y axis
 * 
 * @param mesh the mesh to be rotated
 * @param degrees the degrees on the y axis to rotate the mesh by
 */
void rotate_mesh_y_degrees(Mesh* mesh, double degrees);

/**
 * @brief Rotates a Mesh on the z axis
 * 
 * @param mesh the mesh to be rotated
 * @param degrees the degrees on the z axis to rotate the mesh by
 */
void rotate_mesh_z_degrees(Mesh* mesh, double degrees);

/**
 * @brief Applies a transform to a mesh
 * 
 * @param mesh the mesh to be transformed

 */
void apply_mesh_transform(Mesh* mesh);

/**
 * @brief Computes the bounding box for a given mesh
 * 
 * @param mesh 
 */
void compute_mesh_bounds(Mesh* mesh);

/**
 * @brief Computes the face normals for a given mesh
 * 
 * @param mesh 
 */

void compute_face_normals(Mesh* mesh);

/**
 * @brief Prints the mesh to the console
 * 
 * @param mesh 
 */
void print_mesh(Mesh mesh);

/**
 * @brief Returns a pointer to a mesh loaded from a file
 * 
 * @param filename 
 * @return Mesh* 
 */
Mesh* get_mesh(char* filename);
#endif