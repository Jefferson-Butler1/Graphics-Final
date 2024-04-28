#include "raytrace.h"
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "matrix.h"
#include "FPToolkit.h"
#include "colors.h"
#include "trig.h"
#include "lightmodel.h"

bool SHOW_WORLD_DIRECTION = true;
void invert_show_world_direction(){ SHOW_WORLD_DIRECTION = !SHOW_WORLD_DIRECTION; }
#define SHOW_MISSES 0

int MAX_BOUNCES = 6;


const double EPSILON = 0.000001;

bool intersect_triangle(double* t_out, Vector2* barycentric_out, double closest_t, Ray ray, Triangle triangle){
    //TODO: precompute triangle normal
    //I don't really understand this intuitively. It would be a good idea to go back to this to get a better grasp on it
    Vector3 edge_1 = vec3_sub(triangle.b->position, triangle.a->position);
    Vector3 edge_2 = vec3_sub(triangle.c->position, triangle.a->position);

    // Vector3 normal = vec3_cross_prod(edge_1, edge_2); //This may need to be normalized
    Vector3 normal = triangle.normal;
    Vector3 pvec = vec3_cross_prod(ray.direction, edge_2);
    double determinant = vec3_dot_prod(edge_1, pvec);
    if(fabs(determinant) < EPSILON) return false; //Ray is paralell to triangle
    
    double inverse_determinant = 1.0 / determinant;
    Vector3 vertex_to_origin = vec3_sub(ray.origin, triangle.a->position);

    double u = vec3_dot_prod(vertex_to_origin, pvec) * inverse_determinant;
    if(u < 0 || u > 1) return false;

    Vector3 edge_1_cross_prod = vec3_cross_prod(vertex_to_origin, edge_1);
    double v = vec3_dot_prod(ray.direction, edge_1_cross_prod) * inverse_determinant;

    if(v < 0 || u + v > 1) return false;

    double t = vec3_dot_prod(edge_2, edge_1_cross_prod) * inverse_determinant;
    if(t < closest_t && t > EPSILON) {
        if(t_out != NULL) *t_out = t;
        if(barycentric_out != NULL){
            barycentric_out->x = u;
            barycentric_out->y = v;
        }
        return true;
    }
    return false;
}

bool intersects_bounding_box(Mesh mesh, Ray ray){
    Vector3 box_min = mesh.bounding_box_min;
    Vector3 box_max = mesh.bounding_box_max;

    Vector3 direction_inv = {1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z};
    Vector3 tmin = {(box_min.x - ray.origin.x) * direction_inv.x, (box_min.y - ray.origin.y) * direction_inv.y, (box_min.z - ray.origin.z) * direction_inv.z};
    Vector3 tmax = {(box_max.x - ray.origin.x) * direction_inv.x, (box_max.y - ray.origin.y) * direction_inv.y, (box_max.z - ray.origin.z) * direction_inv.z};

    if (ray.direction.x == 0) {
        tmin.x = ray.origin.x < box_min.x || ray.origin.x > box_max.x ? -INFINITY : INFINITY;
        tmax.x = ray.origin.x < box_min.x || ray.origin.x > box_max.x ? INFINITY : -INFINITY;
    }
    if (ray.direction.y == 0) {
        tmin.y = ray.origin.y < box_min.y || ray.origin.y > box_max.y ? -INFINITY : INFINITY;
        tmax.y = ray.origin.y < box_min.y || ray.origin.y > box_max.y ? INFINITY : -INFINITY;
    }
    if (ray.direction.z == 0) {
        tmin.z = ray.origin.z < box_min.z || ray.origin.z > box_max.z ? -INFINITY : INFINITY;
        tmax.z = ray.origin.z < box_min.z || ray.origin.z > box_max.z ? INFINITY : -INFINITY;
    }

    double t_enter = fmax(fmin(tmin.x, tmax.x), fmax(fmin(tmin.y, tmax.y), fmin(tmin.z, tmax.z)));
    double t_exit = fmin(fmax(tmin.x, tmax.x), fmin(fmax(tmin.y, tmax.y), fmax(tmin.z, tmax.z)));

    return t_enter <= t_exit && t_exit >= 0;
}

void raytrace_scene(int width, int height, Camera cam, 
                    RaytracedParametricObject3D* objs, int num_objs, 
                    Mesh* meshes, int num_meshes, bool skipMeshes,
                    PhongLight* lights, int num_lights, 
                    int numBounces){
    double dwidth = (double)width;
    double dheight = (double)height;
    double film_extent = tan(to_radians(cam.half_fov_degrees));
    for(int y = 0; y < height; y++){
        for(int x = 0; x < width; x++){
            //TODO: make this work for different aspect ratios
            Vector3 pixel_camera_space = {
                ((x - (dwidth / 2)) / dwidth) * (film_extent * 2),
                ((y - (dheight / 2)) / dheight) * (film_extent * 2),
                1
            };

            Vector3 world_space_dir = vec3_sub(mat4_mult_point(pixel_camera_space, cam.inverse_view_matrix), cam.eye);
            Ray ray = {
                .origin=cam.eye,
                .direction=world_space_dir
            };
            if(SHOW_WORLD_DIRECTION){
                G_rgb(SPREAD_COL3(vec3_normalized(world_space_dir)));
                G_pixel(x, y);
            }
            RayHitInfo hit;
            if(raytrace(&hit, ray, numBounces ? numBounces : MAX_BOUNCES, 
                        objs, num_objs, 
                        meshes,num_meshes, false,
                        lights, num_lights)){
                G_rgb(SPREAD_COL3(hit.color));
                G_pixel(x, y);
            }
            else if (SHOW_MISSES){
                G_rgb(.2,.54,.54);
                G_pixel(x, y);
            }
    
        }
    }
}

bool raytrace  (RayHitInfo* out, Ray ray, int depth,
                RaytracedParametricObject3D* objs, int num_objs, 
                Mesh* meshes, int num_meshes, bool skipMeshes,
                PhongLight* lights, int num_lights){
    double closest_t = INFINITY;
    bool did_hit = false;
    if(depth == 0) return false;
    // RayHitInfo result = {
    //     .location={NAN, NAN, NAN},
    //     .normal={NAN, NAN, NAN},
    //     .color={BLACK}
    // };
    Vector3 tip = vec3_add(ray.origin, ray.direction);
    if(!skipMeshes && num_meshes > 0){
        for(int m = 0; m < num_meshes; m++){
            if(meshes[m].hidden) continue;
            Mesh mesh = meshes[m];
            // printf("Mesh %d\tintersects: %d\n", m, intersects_bounding_box(mesh, ray));
            if(!intersects_bounding_box(mesh, ray)) continue;
            for(int t = 0; t < mesh.num_tris; t++){
                Triangle triangle = mesh.tris[t];
                double t;
                Vector2 barycentric;
                if(intersect_triangle(&t, &barycentric, closest_t, ray, triangle)){
                    closest_t = t;
                    if(out == NULL) return true;
                    did_hit = true;

                    Vector3 obj_space_location = vec3_add(ray.origin, vec3_scale(ray.direction, t));
                    out->location = obj_space_location;
                    out->normal = triangle.normal;
                    Vector3 offset_location = vec3_add(out->location, vec3_scale(out->normal, 0.0000000001));
                    vec3_normalize(&out->normal);
                    
                    Vector3 reflection = vec3_reflection(
                        vec3_normalized(vec3_sub(out->location, ray.origin)),
                        out->normal);

                    Ray reflection_ray = {
                        .origin=offset_location,
                        .direction=reflection
                    };

                    RayHitInfo reflections;
                    bool hit = raytrace(&reflections, reflection_ray, depth - 1, objs, num_objs, meshes, num_meshes, false, lights, num_lights);
                    Color3 color = phong_lighting_eye(out->location, out->normal, ray.origin, mesh.material, lights, num_lights);

                    out->color = vec3_add(
                                    vec3_scale(color, mesh.roughness),
                                    vec3_scale(
                                        (hit ? 
                                        reflections.color : 
                                        SHOW_WORLD_DIRECTION ? 
                                            vec3_normalized(reflection_ray.direction) : 
                                            (Color3) {BLACK}
                                        ),
                                        1.0 - mesh.roughness)
                                    );
                }
            }
        }
    }

    for(int o = 0; o < num_objs; o++){
        
        RaytracedParametricObject3D object =  objs[o] ;
        
        Vector3 tsource = mat4_mult_point(ray.origin, object.inverse);
        Vector3 tdir = vec3_sub(mat4_mult_point(tip, object.inverse), tsource);

        // Now we need to foil the terms:
        //(tsrc * tsrc)
        // Vector3 tip_dir = ray.direction; //vec3_sub(ttip, tsource);
        // Vector3 tsource = ray.origin;
        Vector3 va = vec3_mult(tdir, tdir);
        // (tdir * tsrc) + (tsrc * tdir)
        Vector3 vb = vec3_scale(vec3_mult(tdir, tsource), 2);
        // (tdir * tdir)
        Vector3 vc = vec3_mult(tsource, tsource);
        //Add terms together
        double a = va.x + va.y + va.z;
        double b = vb.x + vb.y + vb.z;
        double c = vc.x + vc.y + vc.z - 1;

        //
        double under_sqrt = (b * b) - (4 * a * c);
        if(under_sqrt < 0) continue; // did not intersect with object
        double t_plus = (-b + sqrt(under_sqrt)) / (2 * a);
        double t_minus = (-b - sqrt(under_sqrt)) / (2 * a);

        double t;
        if(t_plus < 0 && t_minus < 0) continue; // object is behind camera
        if(t_plus > 0 && (t_plus < t_minus || t_minus < 0)) t = t_plus;
        else t = t_minus;

        if(t < closest_t && t > 0){
            closest_t = t;
            if(out == NULL) return true;
            did_hit = true;

            Vector3 obj_space_location = vec3_add(tsource, vec3_scale(tdir, t));
            out->location = vec3_add(ray.origin, vec3_scale(ray.direction, t));
            Vector3 obj_normal = vec3_scale(obj_space_location, 2);
            out->normal.x = obj_normal.x * object.inverse[0][0] + obj_normal.y * object.inverse[1][0] + obj_normal.z * object.inverse[2][0];
            out->normal.y = obj_normal.x * object.inverse[0][1] + obj_normal.y * object.inverse[1][1] + obj_normal.z * object.inverse[2][1];
            out->normal.z = obj_normal.x * object.inverse[0][2] + obj_normal.y * object.inverse[1][2] + obj_normal.z * object.inverse[2][2];

            Vector3 offset_location = vec3_add(out->location, vec3_scale(out->normal, 0.0000000001));

            vec3_normalize(&out->normal);
            
            Vector3 reflection = vec3_reflection(
                vec3_normalized(vec3_sub(out->location, ray.origin)),
                out->normal);

            Ray reflection_ray = {
                .origin=offset_location,
                .direction=reflection
            };

            RayHitInfo reflections;
            bool hit = raytrace(&reflections, reflection_ray, depth - 1, objs, num_objs, meshes, num_meshes, false, lights, num_lights);
            Color3 color = phong_lighting_eye(out->location, out->normal, ray.origin, object.material, lights, num_lights);
            out->color = vec3_add(
                            vec3_scale(color, object.roughness),
                            vec3_scale(
                                (hit ? 
                                reflections.color : 
                                // (Color3) {BLACK}
                                reflection_ray.direction
                                ),
                                1.0 - object.roughness)
                            );
        }
    }
    return did_hit;
}
