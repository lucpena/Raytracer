#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <omp.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "model.h"
#include "geometry.h"

// envmap is the SkyBox
int envmap_width, envmap_height;
std::vector<Vec3f> envmap;

// The duck model
Model duck("duck.obj");

struct Light {
    Light(const Vec3f &p, const float &i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

struct Material {
    Material(const float &r, const Vec4f &a, const Vec3f &color, const float &spec) : refractive_index(r), albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Material() : refractive_index(1), albedo(1, 0, 0, 0), diffuse_color(), specular_exponent() {}
    float refractive_index;
    Vec4f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
};


struct Sphere {
    Vec3f center;
    float radius;
    Material material;

    Sphere(const Vec3f &c, const float &r, const Material &m) : center(c), radius(r), material(m) {}

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const{
        Vec3f L = center - orig;
        float tca = L * dir;
        float d2 = L * L - tca * tca;
        if (d2 > radius * radius) return false;
        float thc = sqrt(radius * radius - d2);
              t0 = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true; 
    }
};

Vec3f reflect(const Vec3f &I, const Vec3f &N){
    return I - N * 2.f * (I * N);
}

// Using Snell's law to calculate the refraction
Vec3f refract(const Vec3f &I, const Vec3f &N, const float &refractive_index){
    float cosi = - std::max(-1.f, std::min(1.f, I * N));
    float etai = 1;
    float etat = refractive_index;

    Vec3f n = N;

    // If the ray is inside the object, swap the indices and invert the normal to get the correct result
    if(cosi < 0){
        cosi = -cosi;
        std::swap(etai, etat);
        n = -N;
    }

    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);

    return k < 0 ? Vec3f(0, 0, 0) : I * eta + n * (eta * cosi - sqrt(k));
}

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, Vec3f &hit, Vec3f &N, Material &material){
    float spheres_dist = std::numeric_limits<float>::max();
    for (size_t i = 0; i < spheres.size(); i++){
        float dist_i;

        if(spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist){
            spheres_dist = dist_i;
            hit = orig + dir * dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }

    // Checkerboard plane
    float checkerboard_dist = std::numeric_limits<float>::max();
    if(fabs(dir.y) > 1e-3){
        // The plane has equation y = -4
        float d = -(orig.y + 4)/dir.y;
        Vec3f pt = orig + dir * d;
        if(d > 0 && fabs(pt.x) < 10 && pt.z < -10 && pt.z > -30 && d < spheres_dist){
            checkerboard_dist = d;
            hit = pt;
            N = Vec3f(0, 1, 0);
            material.diffuse_color = (int(.5 * hit.x + 1000) + int(.5 * hit.z)) & 1 ? Vec3f(1, 1, 1) : Vec3f(1, .7, .3);
            material.diffuse_color = material.diffuse_color * .35;
        }
    }

    return std::min(spheres_dist, checkerboard_dist) <1000;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, const std::vector<Vec3f> &envmap, size_t depth = 0){
    Vec3f point, N;
    Material material;

    // Show the background color when there is no hits
    // Depth is for the reflection
    if (depth > 10 || !scene_intersect(orig, dir, spheres, point, N, material)) {

        int x_raw = ((int) ((atan2(dir.z, dir.x) / (2 * M_PI) + 0.5) * envmap_width) + (int) (1.f * envmap_width / (2 * M_PI))) % envmap_width;
        int y_raw = (int)(acos(dir.y) / M_PI * envmap_height);
        int x = std::max(0, std::min(x_raw, envmap_width  - 1));
        int y = std::max(0, std::min(y_raw, envmap_height - 1));
        return envmap[x + y * envmap_width];
    }

    Vec3f reflect_dir = reflect(dir, N).normalize();

    // Offset the original point to avoid occlusion by the object itself
    Vec3f reflect_orig = reflect_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
    Vec3f reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights, envmap, depth + 1);

    // Refraction
    Vec3f refract_dir = refract(dir, N, material.refractive_index).normalize();
    Vec3f refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f refract_color = cast_ray(refract_orig, refract_dir, spheres, lights, envmap, depth + 1);


    float diffuse_light_intensity  = 0;
    float specular_light_intensity = 0;

    for(size_t i = 0; i < lights.size(); i++){
        Vec3f light_dir = (lights[i].position - point).normalize();

        // For the Shadows
        float light_distance = (lights[i].position - point).norm();

        // Check if the point lies in the shadow of the lights[i]
        Vec3f shadow_orig = light_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        if(scene_intersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt - shadow_orig).norm() < light_distance) continue;

        diffuse_light_intensity += lights[i].intensity * std::max(0.f, light_dir * N);
        specular_light_intensity += powf(std::max(0.f, reflect(light_dir, N) * dir), material.specular_exponent) * lights[i].intensity;
    }

    return {
        material.diffuse_color * 
        diffuse_light_intensity * material.albedo[0] + 
        Vec3f(1., 1., 1.) * 
        specular_light_intensity * material.albedo[1] + 
        reflect_color * material.albedo[2] +
        refract_color * material.albedo[3]
    };
}

void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights, const std::vector<Model> &models, const std::vector<Vec3f> &envmap){
    const int width  = 1920 * 3;
    const int height = 1080 * 3;
    const int fov    = M_PI/3.;

    std::vector<Vec3f> framebuffer(width * height);
    
    // The heart of the Ray Tracing
    #pragma omp parallel for
    for(size_t j = 0; j < height; j++){
        for(size_t i = 0; i < width; i++){
            float dir_x =  (i + 0.5) -  width / 2.;
            float dir_y = -(j + 0.5) + height / 2.;
            float dir_z = -height/(2.*tan(fov/2.));
            framebuffer[i + j * width] = cast_ray(Vec3f(0, 0, 0), Vec3f(dir_x, dir_y, dir_z).normalize(), spheres, lights, envmap);
        }
    }

    // Saving the FrameBuffer to a file
    std::ofstream ofs("./out.ppm", std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";

    for(size_t i = 0; i < height * width; ++i){
        for(size_t j = 0; j < 3; j++){
            Vec3f &c = framebuffer[i];
            float max = std::max(c[0], std::max(c[1], c[2]));

            if (max > 1) c = c * (1. / max);

            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }

    ofs.close();
}

int main(){
    int n = -1;
    unsigned char *pixmap = stbi_load("envmap.jpg", &envmap_width, &envmap_height, &n, 0);

    if(!pixmap || 3 != n){
        std::cerr << "Error: can not load the enviroment map (SkyBox)." << std::endl;
        return -1;
    }

    envmap = std::vector<Vec3f>(envmap_width * envmap_height);

    for(int j = envmap_height - 1; j >= 0; j--){
        for(int i = 0; i < envmap_width; i++){
            envmap[i + j * envmap_width] = Vec3f(
                pixmap[(i + j * envmap_width) * 3 + 0],
                pixmap[(i + j * envmap_width) * 3 + 1],
                pixmap[(i + j * envmap_width) * 3 + 2]
            ) * (1 / 255.);
        }
    }

    stbi_image_free(pixmap);

    Material      ivory(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3),   50.);
    Material      glass(1.5, Vec4f(0.0,  0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8),  125.);
    Material red_rubber(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1),   10.);
    Material     mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1425.);

    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(Vec3f(  -3,    0, -16), 2,      ivory));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2,      glass));
    spheres.push_back(Sphere(Vec3f( 1.5, -0.5, -18), 3, red_rubber));
    spheres.push_back(Sphere(Vec3f(   7,    5, -18), 4,      ivory));
    spheres.push_back(Sphere(Vec3f(  -9,    4, -18), 2,     mirror));

    std::vector<Light> lights;
    lights.push_back(Light(Vec3f(-20, 20,  20), 1.5));
    lights.push_back(Light(Vec3f( 30, 50, -25), 1.8));
    lights.push_back(Light(Vec3f( 30, 20,  30), 1.7));

    std::vector<Model> models;
    models.push_back(Model("duck.obj"));

    std::cout << "Gerando imagem. Aguarde..." << std::endl;

    render(spheres, lights, models, envmap);

    std::cout << "Imagem gerada com sucesso! verificar out.ppm" << std::endl;
    return 0;
}
