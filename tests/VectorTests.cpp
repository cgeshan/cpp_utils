#include "../src/cgVector.h"

void testVectorConstructors()
{
    // Test the default constructor
    Vector<int> v1;
    assert(v1.getX() == 0 && v1.getY() == 0 && v1.getZ() == 0);

    // Test the parameterized constructor
    Vector<double> v2(1.0, 2.0, 3.0);
    assert(v2.getX() == 1.0 && v2.getY() == 2.0 && v2.getZ() == 3.0);
}

void testVectorOperations()
{
    // Create some vectors for testing
    Vector<double> v1(1.0, 2.0, 3.0);
    Vector<double> v2(2.0, 3.0, 4.0);

    // Test vector addition
    Vector<double> result1 = v1 + v2;
    assert(result1.getX() == 3.0 && result1.getY() == 5.0 && result1.getZ() == 7.0);

    // Test vector subtraction
    Vector<double> result2 = v1 - v2;
    assert(result2.getX() == -1.0 && result2.getY() == -1.0 && result2.getZ() == -1.0);

    // Test vector scalar multiplication
    Vector<double> result3 = v1 * 2.0;
    assert(result3.getX() == 2.0 && result3.getY() == 4.0 && result3.getZ() == 6.0);

    // Test dot product
    double dotProduct = v1.dot(v2);
    assert(fabs(dotProduct - 20.0) < 1e-6);

    // Test cross product
    Vector<double> crossProduct = v1.cross(v2);
    assert(crossProduct.getX() == -1.0 && crossProduct.getY() == 2.0 && crossProduct.getZ() == -1.0);

    // Test magnitude
    double magnitude = v1.magnitude();
    assert(fabs(magnitude - 3.741657) < 1e-6);

    // Test vector normalization
    Vector<double> normalized = v1.normalize();
    double mag = normalized.magnitude();
    assert(fabs(mag - 1.0) < 1e-6);

    // Test projection
    Vector<double> projection = v1.projection(v2);
    assert(projection.getX() == 1.25 && projection.getY() == 2.0 && projection.getZ() == 2.75);

    // Test isOrthogonal
    bool orthogonal = v1.isOrthogonal(v2);
    assert(!orthogonal);

    // Test isParallel
    bool parallel = v1.isParallel(v2);
    assert(!parallel);

    // Test distance
    double distance = v1.distance(v2);
    assert(fabs(distance - 1.732051) < 1e-6);

    // Test midpoint
    Vector<double> mid = v1.midpoint(v2);
    assert(mid.getX() == 1.5 && mid.getY() == 2.5 && mid.getZ() == 3.5);

    // Test isZero
    Vector<int> v3(0, 0, 0);
    assert(v1.isZero() == false);
    assert(v3.isZero() == true);
}

int main()
{
    testVectorConstructors();
    testVectorOperations();

    std::cout << "All tests passed." << std::endl;

    return 0;
}
