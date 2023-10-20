#include "../src/cgVector.h"

void testVectorConstructors()
{
    std::cout << "Test 1: Default Constructor" << std::endl;
    // Test the default constructor
    Vector<int> v1;
    assert(v1.getX() == 0 && v1.getY() == 0 && v1.getZ() == 0);

    std::cout << "Test 2: Parameterized Constructor" << std::endl;
    // Test the parameterized constructor
    Vector<double> v2(1.0, 2.0, 3.0);
    assert(v2.getX() == 1.0 && v2.getY() == 2.0 && v2.getZ() == 3.0);
}

void testVectorOperations()
{
    std::cout << "Test 3: Vector Addition" << std::endl;
    // Create some vectors for testing
    Vector<double> v1(1.0, 2.0, 3.0);
    Vector<double> v2(2.0, 3.0, 4.0);

    // Test vector addition
    Vector<double> result1 = v1 + v2;
    assert(result1.getX() == 3.0 && result1.getY() == 5.0 && result1.getZ() == 7.0);

    std::cout << "Test 4: Vector Subtraction" << std::endl;
    // Test vector subtraction
    Vector<double> result2 = v1 - v2;
    assert(result2.getX() == -1.0 && result2.getY() == -1.0 && result2.getZ() == -1.0);

    std::cout << "Test 5: Vector Scalar Multiplication" << std::endl;
    // Test vector scalar multiplication
    Vector<double> result3 = v1 * 2.0;
    assert(result3.getX() == 2.0 && result3.getY() == 4.0 && result3.getZ() == 6.0);

    std::cout << "Test 6: Dot Product" << std::endl;
    // Test dot product
    double dotProduct = v1.dot(v2);
    assert(fabs(dotProduct - 20.0) < 1e-6);

    std::cout << "Test 7: Cross Product" << std::endl;
    // Test cross product
    Vector<double> crossProduct = v1.cross(v2);
    assert(crossProduct.getX() == -1.0 && crossProduct.getY() == 2.0 && crossProduct.getZ() == -1.0);

    std::cout << "Test 8: Magnitude" << std::endl;
    // Test magnitude
    double magnitude = v1.magnitude();
    assert(fabs(magnitude - 3.741657) < 1e-6);

    std::cout << "Test 9: Vector Normalization" << std::endl;
    // Test vector normalization
    Vector<double> normalized = v1.normalize();
    double mag = normalized.magnitude();
    assert(fabs(mag - 1.0) < 1e-6);

    std::cout << "Test 10: Projection" << std::endl;
    // Test projection
    Vector<double> projection = v1.projection(v2);
    assert(fabs(projection.getX() - 1.3793 < 1e-4) && fabs(projection.getY() - 2.0689 < 1e-4) && fabs(projection.getZ() - 2.7586 < 1e-4));

    std::cout << "Test 11: IsOrthogonal" << std::endl;
    // Test isOrthogonal
    bool orthogonal = v1.isOrthogonal(v2);
    assert(!orthogonal);

    std::cout << "Test 12: IsParallel" << std::endl;
    // Test isParallel
    bool parallel = v1.isParallel(v2);
    assert(!parallel);

    std::cout << "Test 13: Distance" << std::endl;
    // Test distance
    double distance = v1.distance(v2);
    assert(fabs(distance - 1.732051) < 1e-6);

    std::cout << "Test 14: Midpoint" << std::endl;
    // Test midpoint
    Vector<double> mid = v1.midpoint(v2);
    assert(mid.getX() == 1.5 && mid.getY() == 2.5 && mid.getZ() == 3.5);

    std::cout << "Test 15: IsZero" << std::endl;
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
