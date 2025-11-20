#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "raytracer.h"
#include "hittable.h"
#include "aabb.h"

class triangle : public hittable {
public:
	// Stationary Triangle
	triangle(const point3& v0, const point3& v1, const point3& v2, shared_ptr<material> mat)
		: v0(v0), v1(v1), v2(v2), mat(mat)
	{
		// Build aabb from min/max of vertices, expand slightly to avoid degenerate boxes
		const double eps = 1e-4;
		double minx = std::fmin(std::fmin(v0.x(), v1.x()), v2.x());
		double miny = std::fmin(std::fmin(v0.y(), v1.y()), v2.y());
		double minz = std::fmin(std::fmin(v0.z(), v1.z()), v2.z());
		double maxx = std::fmax(std::fmax(v0.x(), v1.x()), v2.x());
		double maxy = std::fmax(std::fmax(v0.y(), v1.y()), v2.y());
		double maxz = std::fmax(std::fmax(v0.z(), v1.z()), v2.z());

		point3 bmin(minx - eps, miny - eps, minz - eps);
		point3 bmax(maxx - eps, maxy - eps, maxz - eps);
		bbox = aabb(bmin, bmax);
	}

	bool hit(const ray& r, interval ray_t, hit_record& rec) const override {
		// M�ller�Trumbore intersection
		const double EPS = 1e-8;
		auto edge1 = v1 - v0;
		auto edge2 = v2 - v0;

		auto h = cross(r.direction(), edge2);
		auto a = dot(edge1, h);

		if (std::fabs(a) < EPS) // ray parallel or nearly parallel to triangle
			return false;

		auto f = 1.0 / a;
		auto s = r.origin() - v0;
		auto u = f * dot(s, h);

		if (u < 0.0 || u > 1.0)
			return false;

		auto q = cross(s, edge1);
		auto v = f * dot(r.direction(), q);

		if (v < 0.0 || (u + v) > 1.0)
			return false;

		auto t = f * dot(edge2, q);

		if (!ray_t.surrounds(t))
			return false;

		//fill hit record
		rec.t = t;
		rec.p = r.at(t);
		rec.u = u;
		rec.v = v;
		vec3 outward_normal = unit_vector(cross(edge1, edge2));
		rec.set_face_normal(r, outward_normal);
		rec.mat = mat;

		return true;
	}

	aabb bounding_box() const override { return bbox; }

private:
	point3 v0;
	point3 v1;
	point3 v2;
	shared_ptr<material> mat;
	aabb bbox;
};

#endif