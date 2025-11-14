// RayTracer.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "vectorlib.hpp"
#include <vector>
#include <utility>
#include <algorithm>
#include <limits>
#include <cmath>
#include <thread>
#include <random>
#include <omp.h>

#define WIDTH 1000
#define HEIGHT 1000

constexpr float PI = 3.14159265358979323846f;

// Thread-local random float in [0,1)
inline float random_float() {
    static thread_local std::mt19937 generator(std::random_device{}());
    static thread_local std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
    return distribution(generator);
}

struct Rayon {
	Point origin;
	Direction direction;
};

struct IntersectionResult {
	float t;
	Vec3 color;
	std::string type;
	Direction normal; // Normal at intersection (added to avoid dynamic_cast in the hot loop)
};

struct SceneObject {
	virtual ~SceneObject() = default;
	virtual IntersectionResult intersect(const Rayon& r) const = 0;
};

struct Sphere : public SceneObject {
	float x, y, z; // Centre de la sphère
	float radius; // Rayon de la sphère
	Vec3 color; // Couleur de la sphère

	Sphere(float x_, float y_, float z_, float radius_, Vec3 color_)
		: x(x_), y(y_), z(z_), radius(radius_), color(color_) {
	}

	IntersectionResult intersect(const Rayon& r) const override {
		Direction oc = r.origin - Point{ x, y, z };
		float a = r.direction._dx * r.direction._dx + r.direction._dy * r.direction._dy + r.direction._dz * r.direction._dz;
		float b = 2.0f * (oc._dx * r.direction._dx + oc._dy * r.direction._dy + oc._dz * r.direction._dz);
		float c = oc._dx * oc._dx + oc._dy * oc._dy + oc._dz * oc._dz - radius * radius;
		float discriminant = b * b - 4 * a * c;
		if (discriminant > 0) {
			float t1 = (-b - sqrt(discriminant)) / (2.0f * a);
			float t2 = (-b + sqrt(discriminant)) / (2.0f * a);
			float t = (t1 > 0) ? t1 : ((t2 > 0) ? t2 : -1.0f);
			if (t > 0) {
				// compute normal at hit point
				float hx = r.origin._x + r.direction._dx * t;
				float hy = r.origin._y + r.direction._dy * t;
				float hz = r.origin._z + r.direction._dz * t;
				Direction normal = Direction(hx - x, hy - y, hz - z).normalize();
				return { t, color, "sphere", normal };
			}
		}
		return { -1.0f, Vec3{0,0,0}, "sphere", Direction{0,0,0} }; // Pas d'intersection
	}
};

struct Plane : public SceneObject {
	float x, y, z; // Un point du plan
	Direction normale; // La normale du plan (direction perpendiculaire)
	Vec3 color; // Couleur du plan

	Plane(float x_, float y_, float z_, Direction normale_, Vec3 color_)
		: x(x_), y(y_), z(z_), normale(normale_), color(color_) {
	}

	IntersectionResult intersect(const Rayon& r) const override {
		Direction n = normale.normalize();
		Point p0{ x, y, z };
		float denom = n._dx * r.direction._dx + n._dy * r.direction._dy + n._dz * r.direction._dz;
		if (fabs(denom) > 1e-6) {
			Direction p0l0 = p0 - r.origin;
			float t = (p0l0._dx * n._dx + p0l0._dy * n._dy + p0l0._dz * n._dz) / denom;
			if (t >= 0) {
				return { t, color, "plane", n };
			}
		}
		return { -1.0f, Vec3{0,0,0}, "plane", Direction{0,0,0} }; // Pas d'intersection
	}
};

struct Light : SceneObject {
	float x, y, z; // Position de la lumière
	Vec3 intensity; // Intensité lumineuse (RGB)

	Light(float x_, float y_, float z_, Vec3 inten) : x(x_), y(y_), z(z_), intensity(inten) {}

	IntersectionResult intersect(const Rayon& r) const override {
		return { -1.0f, Vec3{0,0,0}, "light", Direction{0,0,0} }; // Les lumières ne sont pas des objets à intersecter
	}
};

struct PixelRGB {
	float r, g, b;
	static PixelRGB Black() {
		return PixelRGB{ 0, 0, 0 };
	}
	static PixelRGB Red() {
		return PixelRGB{ 255, 0, 0 };
	}
	static PixelRGB Green() {
		return PixelRGB{ 0, 255, 0 };
	}
	static PixelRGB Blue() {
		return PixelRGB{ 0, 0, 255 };
	}

	//distance to color mapping function
	static PixelRGB distanceToColor(float distance, Vec3 objectColor) {
		if (distance < 0) return PixelRGB{ 0, 0, 0 }; // Pas d'intersection
		float intensity = std::max(0.0f, 1.0f - distance / 150.0f);
		return PixelRGB{ intensity * objectColor.unsafeIndex(0), intensity * objectColor.unsafeIndex(1), intensity * objectColor.unsafeIndex(2) }; // Object color based on distance
	}
};

// mkCameraRayon now prend des offsets sub-pixel (ox, oy) en [0,1] pour l'anti-aliasing
struct mkCameraRayon {
	float line, column;
	float ox, oy; // offset subpixel horizontal et vertical (0..1)
	mkCameraRayon(float l, float c, float ox_ = 0.5f, float oy_ = 0.5f) : line(l), column(c), ox(ox_), oy(oy_) {}

	Rayon operator()() {
		static const float fov = 80.0f;
		static const float aspectRatio = float(WIDTH) / float(HEIGHT);
		static const float factor = tan(fov / 2 * PI / 180);
		// Utiliser offsets sub-pixel pour calculer la position du rayon dans l'écran
		float px = (2 * ((column + ox) / WIDTH) - 1) * factor * aspectRatio;
		float py = (1 - 2 * ((line + oy) / HEIGHT)) * factor;
		Point origin{ 0, 0, 50 };
		Direction direction{ px, py, -1 };
		direction = direction.normalize();
		return Rayon{ origin, direction };
	}
};

std::pair<float, Vec3> raytrace(const Rayon& r, const std::vector<SceneObject*>& scene) {
	float closest_t = std::numeric_limits<float>::max();
	Vec3 hit_color = { 0,0,0 };
	for (const SceneObject* obj : scene) {
		IntersectionResult result = obj->intersect(r);
		float t = result.t;
		Vec3 color = result.color;
		if (t > 0 && t < closest_t) {
			closest_t = t;
			hit_color = color; // <-- Garde la bonne couleur
		}
	}
	return { closest_t == std::numeric_limits<float>::max() ? -1.0f : closest_t, hit_color };
};

//maximum pixel value in the image
float maximumPixel(const std::vector<PixelRGB>& image) {
	float maxPixel = 0.0f;
	for (const PixelRGB& pixel : image) {
		maxPixel = std::max(maxPixel, std::max(pixel.r, std::max(pixel.g, pixel.b)));
	}
	return maxPixel;
}

Vec3 computeDirectLighting(
	const Point& hitPoint,
	const Direction& normal,
	const Vec3& objectColor,
	const std::vector<SceneObject*>& scene,
	const Light* light)
{
	// Direction vers la lumière
	Direction toLight(
		light->x - hitPoint._x,
		light->y - hitPoint._y,
		light->z - hitPoint._z
	);
	float distanceToLight = toLight.length();
	toLight = toLight.normalize();

	// Ombre : lancer un rayon vers la lumière
	float epsilon = 1e-3f;
	Point shadowOrigin(
		hitPoint._x + normal._dx * epsilon,
		hitPoint._y + normal._dy * epsilon,
		hitPoint._z + normal._dz * epsilon
	);
	Rayon shadowRay{ shadowOrigin, toLight };

	// Vérifier s'il y a un objet entre le point et la lumière
	for (const SceneObject* obj : scene) {
		if (dynamic_cast<const Light*>(obj)) continue;
		IntersectionResult shadowResult = obj->intersect(shadowRay);
		if (shadowResult.t > 0 && shadowResult.t < distanceToLight) {
			// Dans l'ombre
			return Vec3{ 0, 0, 0 };
		}
	}

	// Éclairage Lambertien
	float NdotL = std::max(0.0f, normal._dx * toLight._dx + normal._dy * toLight._dy + normal._dz * toLight._dz);
	float attenuation = 1.0f / (0.1f * distanceToLight); // Atténuation avec la distance
	return Vec3{
		objectColor._x * light->intensity._x / 255.0f * NdotL * attenuation,
		objectColor._y * light->intensity._y / 255.0f * NdotL * attenuation,
		objectColor._z * light->intensity._z / 255.0f * NdotL * attenuation
	};
}

// Trace un seul rayon et renvoie la couleur (avant normalisation globale)
Vec3 traceRay(const Rayon& rayonDepuisLePixel, const std::vector<SceneObject*>& scene) {
	float closest_t = std::numeric_limits<float>::max();
	Vec3 hitColor = { 0,0,0 };
	Direction hitNormal(0, 0, 0);
	Point hitPoint(0, 0, 0);
	bool hit = false;
	std::string hitType;

	// Recherche de l'intersection la plus proche
	for (const SceneObject* obj : scene) {
		IntersectionResult result = obj->intersect(rayonDepuisLePixel);
		if (result.t > 0 && result.t < closest_t) {
			closest_t = result.t;
			hitColor = result.color;
			hitNormal = result.normal;
			hitType = result.type;
			hitPoint = Point(
				rayonDepuisLePixel.origin._x + rayonDepuisLePixel.direction._dx * closest_t,
				rayonDepuisLePixel.origin._y + rayonDepuisLePixel.direction._dy * closest_t,
				rayonDepuisLePixel.origin._z + rayonDepuisLePixel.direction._dz * closest_t
			);
			hit = true;
		}
	}

	if (!hit) return Vec3{ 0,0,0 };

	// Calcul de l'éclairage direct par les lumières de la scène
	Vec3 finalColor = { 0,0,0 };
	for (const SceneObject* obj : scene) {
		const Light* light = dynamic_cast<const Light*>(obj);
		if (light) {
			Vec3 contrib = computeDirectLighting(hitPoint, hitNormal, hitColor, scene, light);
			finalColor._x += contrib._x;
			finalColor._y += contrib._y;
			finalColor._z += contrib._z;
		}
	}

	return finalColor;
}

//int main()
//{
//	std::vector<SceneObject*> scene;
//
//	//scene.push_back(new Sphere{ 5, 0, -95, 30, {255,0,0 } }); // sphère rouge au fond
//	//scene.push_back(new Sphere{ 0, 0, -80, 20, {0,255,0 } }); // sphère verte au fond
//	scene.push_back(new Sphere{ -5, 0, -50, 10, {0,0,255 } }); // Petite sphère bleue devant
//	//scene.push_back(new Sphere{ -5, 0, -40, 5, {255,0,0 } }); // Petite sphère rouge devant
//	scene.push_back(new Sphere{ 50, 10, -60, 10, {50,100,200 } });	// Petite sphère bleue
//	scene.push_back(new Sphere{ -60, -20, -50, 8, {230,10,10 } }); // Petite sphère rouge
//	scene.push_back(new Sphere{ 30, -50, -70, 30, {0,255,255 } }); // Grande sphère cyan
//	scene.push_back(new Sphere{ -30, 35, -40, 5, {255,255,255 } }); // Petite sphère blanche
//	scene.push_back(new Sphere{ 20, 58, -80, 40, {255,255,0 } }); // Sphere jaune en haut
//	scene.push_back(new Plane{ 0, 60, 0, {0,1,0}, {255,0,0} }); // Mur plafond
//	scene.push_back(new Plane{ 0, -60, 0, {0,1,0}, {255,255,0} }); // Mur sol
//	scene.push_back(new Plane{ 60, 0, 0, {1,0,0}, {0,0,255} }); //Mur droite
//	scene.push_back(new Plane{ -60, 0, 0, {1,0,0}, {0,255,0} }); //Mur gauche
//	scene.push_back(new Plane{ 0, 0, -80, {0,0,1}, {255,255,255} }); //Mur fond
//	scene.push_back(new Light{ -40, 20, -10, {500,500,500} }); // Lumière gauche
//
//
//	std::vector<PixelRGB> image(WIDTH * HEIGHT, PixelRGB::Black());
//
//	// Samples per pixel for anti-aliasing
//	constexpr int SAMPLES_PER_PIXEL = 16; // e.g., 16 samples per pixel
//
//	// Boucle principale : pour chaque pixel, lancer SAMPLES_PER_PIXEL rayons et faire la moyenne
//	const int samplesPerPixel = SAMPLES_PER_PIXEL;
//
//	// Determine number of threads
//	unsigned int hwThreads = std::thread::hardware_concurrency();
//	if (hwThreads == 0) hwThreads = 4;
//	unsigned int numThreads = std::min<unsigned int>(hwThreads, HEIGHT);
//
//	// Render using OpenMP parallel for over scanlines
//	// Each iteration of the outer loop (line) is independent
//	#pragma omp parallel for schedule(static)
//	for (int line = 0; line < HEIGHT; ++line) {
//		for (int column = 0; column < WIDTH; ++column)
//		{
//			Vec3 accumColor{ 0,0,0 };
//
//			// Random sampling inside the pixel area
//			for (int s = 0; s < samplesPerPixel; ++s) {
//				float ox = random_float();
//				float oy = random_float();
//
//				Rayon rayon = mkCameraRayon((float)line, (float)column, ox, oy)();
//				Vec3 sampleColor = traceRay(rayon, scene);
//				accumColor._x += sampleColor._x;
//				accumColor._y += sampleColor._y;
//				accumColor._z += sampleColor._z;
//			}
//
//			// Moyenne des échantillons
//			Vec3 finalColor = {
//				accumColor._x / float(samplesPerPixel),
//				accumColor._y / float(samplesPerPixel),
//				accumColor._z / float(samplesPerPixel)
//			};
//
//			image[line * WIDTH + column] = PixelRGB{
//				std::min(finalColor._x, 255.0f),
//				std::min(finalColor._y, 255.0f),
//				std::min(finalColor._z, 255.0f)
//			};
//		}
//	}
//
//	float maxPixel = maximumPixel(image);
//	if (maxPixel == 0) maxPixel = 1; // Protection contre la division par zéro
//
//	std::vector<unsigned char> imgData(WIDTH * HEIGHT * 3);
//	for (int i = 0; i < WIDTH * HEIGHT; ++i) {
//		// Apply simple gamma correction (gamma=2) after normalization
//		float rnorm = (image[i].r / maxPixel);
//		float gnorm = (image[i].g / maxPixel);
//		float bnorm = (image[i].b / maxPixel);
//		rnorm = std::sqrt(std::max(0.0f, rnorm));
//		gnorm = std::sqrt(std::max(0.0f, gnorm));
//		bnorm = std::sqrt(std::max(0.0f, bnorm));
//
//		imgData[i * 3 + 0] = static_cast<unsigned char>(std::min(255.0f, rnorm * 255.0f));
//		imgData[i * 3 + 1] = static_cast<unsigned char>(std::min(255.0f, gnorm * 255.0f));
//		imgData[i * 3 + 2] = static_cast<unsigned char>(std::min(255.0f, bnorm * 255.0f));
//	}
//	int result = stbi_write_png("C:/Users/abasei/Documents/Gamagora/output.png", WIDTH, HEIGHT, 3, imgData.data(), WIDTH * 3);
//	if (result)
//		std::cout << "Image sauvegardee avec succes !" << std::endl;
//	else
//		std::cout << "Erreur lors de la sauvegarde de l'image !" << std::endl;
//
//	// Libération mémoire sommaire
//	for (SceneObject* obj : scene) delete obj;
//}