#include <iostream>

struct Vec3
{
    float _x, _y, _z;

    Vec3(float x, float y, float z) : _x(x), _y(y), _z(z)
    {
    }

    Vec3 operator+(const Vec3& other) const // Const pour que l'objet ne soit pas modifié
    {
        return { _x + other._x, _y + other._y, _z + other._z };
    }

    float dot(const Vec3& other) const
    {
        return _x * other._x + _y * other._y + _z * other._z;
    }

    float unsafeIndex(int i) const
    {
        switch (i)
        {
        case 0:
            return _x;
        case 1:
            return _y;
        case 2:
            return _z;
        default:
            throw "out of bound";
        }
    }

    float length() const
    {
        return std::sqrt(this->dot(*this));
    }
};

struct Direction
{
    float _dx, _dy, _dz;

    Direction(float dx, float dy, float dz) : _dx(dx), _dy(dy), _dz(dz)
    {
    }

    Direction normalize() const
    {
        float norm = std::sqrt(_dx * _dx + _dy * _dy + _dz * _dz);
        if (norm == 0)
        {
            return Direction{ 0, 0, 0 };
        }
        return Direction{ _dx / norm, _dy / norm, _dz / norm };
    }

    // Retourne le produit scalaire de deux directions.
    // Si le res est positif, les deux directions sont dans le même sens.
    // Si le res est négatif, les deux directions sont dans le sens opposé.
    // Si le res est nul, les deux directions sont orthogonales (perpendiculaires).
    float scalar(Direction d1, Direction d2)
    {
        return d1._dx * d2._dx + d1._dy * d2._dy + d1._dz * d2._dz;
    }

    float length() const
    {
        return std::sqrt(_dx * _dx + _dy * _dy + _dz * _dz);
	}
};

struct Point
{
    float _x, _y, _z;

    Point(float x, float y, float z) : _x(x), _y(y), _z(z)
    {
    }

    Point operator+(const Point& other) const
    {
        return { _x + other._x, _y + other._y, _z + other._z };
    }

    Direction operator-(const Point& other) const
    {
        return { _x - other._x, _y - other._y, _z - other._z };
    }

    Direction translate(const Direction& d) const
    {
        return { _x + d._dx, _y + d._dy, _z + d._dz };
    }

    float distance(const Point& other) const
    {
        Direction d = *this - other;
        return std::sqrt(d._dx * d._dx + d._dy * d._dy + d._dz * d._dz);
    }

    Point lerp(const Point& other, float t) const
    {
        return Point((1 - t) * _x + t * other._x, (1 - t) * _y + t * other._y, (1 - t) * _z + t * other._z);
    }

    void print() const
    {
        std::cout << "Point(" << _x << ", " << _y << ", " << _z << ")" << std::endl;
    }
};

struct QuantiteLumiere
{
    float _r, _g, _b;
};

struct Surface
{
    float _r, _g, _b;

    Surface(float r, float g, float b) : _r(r), _g(g), _b(b)
    {
    }

    Surface illuminate(const Surface& surface, float intensity) const
    {
        return Surface{
            std::min(_r * surface._r * intensity / 100.0f, 100.0f),
            std::min(_g * surface._g * intensity / 100.0f, 100.0f),
            std::min(_b * surface._b * intensity / 100.0f, 100.0f) };
    }
};

static QuantiteLumiere Black()
{
    return QuantiteLumiere{
        0,
        0,
        0 };
}
static QuantiteLumiere White()
{
    return QuantiteLumiere{
        100,
        100,
        100 };
}
static QuantiteLumiere Red()
{
    return QuantiteLumiere{
        100,
        0,
        0 };
}
static QuantiteLumiere Green()
{
    return QuantiteLumiere{
        0,
        100,
        0 };
}
static QuantiteLumiere Blue()
{
    return QuantiteLumiere{
        0,
        0,
        100 };
}

//int main()
//{
//    Vec3 v{ 0, 1, 2 };
//
//    v.unsafeIndex(2);
//
//    v.length();
//
//    Direction d{ 1, 2, 3 };
//    d.normalize();
//
//    //afficher d dans la console
//
//    std::cout << "Direction(" << d._dx << ", " << d._dy << ", " << d._dz << ")" << std::endl;
//};

//Pourquoi le code ne s'exécute pas dans visual studio ? "Sélectionner un élément de démarrage valide"

