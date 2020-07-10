#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <string>

#include "geometry.h"

class Model {

private:
    std::vector<Vec3f> verts;
    std::vector<Vec3i> faces;

public:
    Model(const char* filename);

    // Number of vertices
    int nverts() const;

    // Nuber of Triangles
    int nfaces() const;

    bool ray_triangle_intersect(const int &fi, const Vec3f &orig, const Vec3f &dir, float &tnear);

    // Coordinates of the vertex i
    const Vec3f &point(int i) const;
    Vec3f &point(int i);

    // Index of the vertex fot the triangle fi and local index li
    int vert(int fi, int li) const;

    // Bounding box for all the vertices, includinh isolated ones
    void get_bbox(Vec3f &min, Vec3f &max);

};

std::ostream& operator << (std::ostream& out, Model &m);

#endif // MODEL_H