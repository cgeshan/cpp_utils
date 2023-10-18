#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cmath>

template <typename T>
class Vector
{
private:
    T x, y, z;

public:
    // Constructors
    Vector()
    {
        x = 0;
        y = 0;
        z = 0;
    }

    Vector(T x, T y, T z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    // Public getters for x, y, and z
    T getX() const
    {
        return x;
    }

    T getY() const
    {
        return y;
    }

    T getZ() const
    {
        return z;
    }

    // Public setters for x, y, and z
    void setX(T x)
    {
        this->x = x;
    }

    void setY(T y)
    {
        this->y = y;
    }

    void setZ(T z)
    {
        this->z = z;
    }

    // Vector addition
    Vector operator+(const Vector &vec) const
    {
        return Vector(x + vec.x, y + vec.y, z + vec.z);
    }

    // Vector subtraction
    Vector operator-(const Vector &vec) const
    {
        return Vector(x - vec.x, y - vec.y, z - vec.z);
    }

    // Scalar multiplication
    Vector operator*(T scalar) const
    {
        return Vector(x * scalar, y * scalar, z * scalar);
    }

    // Dot product of vectors
    T dot(const Vector &vec) const
    {
        return x * vec.x + y * vec.y + z * vec.z;
    }

    // Cross product of vectors
    Vector cross(const Vector &vec) const
    {
        return Vector(
            y * vec.z - z * vec.y,
            z * vec.x - x * vec.z,
            x * vec.y - y * vec.x);
    }

    // Magnitude (length) of the vector
    T magnitude() const
    {
        return std::sqrt(x * x + y * y + z * z);
    }

    // Normalize the vector (convert to a unit vector)
    Vector normalize() const
    {
        T mag = magnitude();
        if (mag == 0)
            return *this; // Avoid division by zero
        return Vector(x / mag, y / mag, z / mag);
    }

    // Print the vector
    void print() const
    {
        std::cout << "(" << x << ", " << y << ", " << z << ")" << std::endl;
    }

    // Angle between two vectors in radians
    double angle(const Vector &vec) const
    {
        double magProduct = magnitude() * vec.magnitude();
        if (magProduct == 0)
            return 0.0; // Avoid division by zero
        return std::acos(dot(vec) / magProduct);
    }

    // Projection of the vector onto another vector
    Vector projection(const Vector &vec) const
    {
        double vecMagSquared = vec.magnitude() * vec.magnitude();
        if (vecMagSquared == 0)
            return Vector(); // Avoid division by zero
        double scaleFactor = dot(vec) / vecMagSquared;
        return vec * scaleFactor;
    }

    // Check if two vectors are orthogonal (perpendicular)
    bool isOrthogonal(const Vector &vec) const
    {
        return dot(vec) == 0;
    }

    // Check if two vectors are parallel
    bool isParallel(const Vector &vec) const
    {
        if (magnitude() == 0 || vec.magnitude() == 0)
            return true; // Handle zero vectors
        double angleRad = angle(vec);
        return angleRad == 0 || angleRad == M_PI;
    }

    // Euclidean distance between two vectors
    T distance(const Vector &vec) const
    {
        Vector diff = *this - vec;
        return diff.magnitude();
    }

    // Midpoint between two vectors
    Vector midpoint(const Vector &vec) const
    {
        return (*this + vec) * 0.5;
    }

    // Check if the vector is a zero vector
    bool isZero() const
    {
        return x == 0 && y == 0 && z == 0;
    }
};

#endif
