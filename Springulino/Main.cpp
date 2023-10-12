#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>


const double EPSILON = 1e-6;


struct Vec3 {
    double X, Y, Z;

    Vec3(double x, double y, double z) : X(x), Y(y), Z(z) {};

    Vec3() : X(0), Y(0), Z(0) {};

    double Magnitude() {
        return pow((X * X + Y * Y + Z * Z), 0.5);
    }

    Vec3& operator+=(const Vec3& other) {
        this->X += other.X;
        this->Y += other.Y;
        this->Z += other.Z;
        return *this;
    }

    Vec3& operator=(const Vec3& other) {
        if (this != &other) {
            this->X = other.X;
            this->Y = other.Y;
            this->Z = other.Z;
        }
        return *this;
    }

    Vec3& operator+=(double other) {
        this->X += other;
        this->Y += other;
        this->Z += other;
        return *this;
    }

    Vec3 operator-(const Vec3& other) const {
        return Vec3(X - other.X, Y - other.Y, Z - other.Z);
    }

    Vec3 operator*(double num) const {
        return Vec3(X * num, Y * num, Z * num);
    }

    Vec3 operator/(double num) const {
        return Vec3(X / num, Y / num, Z / num);
    }

    friend Vec3 operator/(double scalar, const Vec3& v) {
        if (v.X != 0 && v.Y != 0 && v.Z != 0) {
            return Vec3(scalar / v.X, scalar / v.Y, scalar / v.Z);
        }
        else {
            std::cerr << "Division by zero is not allowed." << std::endl;
            return v;
        }
    }

    Vec3& operator-(double num) {
        Vec3 newVec3(0, 0, 0);

        newVec3.X -= num;
        newVec3.Y -= num;
        newVec3.Z -= num;
        return newVec3;
    }

    Vec3 operator+(double num) {
        return Vec3(X+num, Y+num, Z+num);
    }

    Vec3 operator+(const Vec3& other) {
        return Vec3(this->X + other.X, this->Y + other.Y, this->Z + other.Z);
    }

    Vec3 operator-() const {
        return Vec3(-X, -Y, -Z);
    }
};

double dot(const Vec3& v1, const Vec3& v2) {
    return v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;
}

Vec3 cross(const Vec3& v1, const Vec3& v2) {
    return Vec3(
        v1.Y * v2.Z - v1.Z * v2.Y,
        v1.Z * v2.X - v1.X * v2.Z,
        v1.X * v2.Y - v1.Y * v2.X
    );
}

class Mass {
public:
    Mass(double mass, Vec3 initial_position, Vec3 initial_velocity = Vec3(0, 0, 0)) {
        this->mass = mass;
        this->pos = initial_position;
        this->vel = initial_velocity;
    }

    void applyForce(Vec3 force, double dt) {
        vel += (force / mass) * dt;
    }

    void Move(double dt) {
        pos += vel * dt;
    }

    Vec3 getPosition() {
        return pos;
    }

    double getMass() {
        return mass;
    }

    double mass;
    Vec3 pos;
    Vec3 vel;
    Vec3 accel;
    Vec3 force;
};

class Spring {
public:
    Spring(double spring_constant, double equilibrium_position) {
        this->k = spring_constant;
        this->l0 = equilibrium_position;
    }

    Vec3 calculateForce(const Vec3& dir, double distance) {
        Vec3 force;

        double F = -k * (distance - l0);

        double dirRatio = pow(abs(dir.X * dir.X + dir.Y * dir.Y) + dir.Z * dir.Z, 0.5) / F;

        force = dir / dirRatio;

        return force;
    }

    double k;
    double l0;
};

// Masses, Spring pair
struct MSPair {
    Mass* firstMass;
    Mass* secondMass;
    Spring* spring;

    double Distance() {
        Vec3 dvec = secondMass->pos - firstMass->pos;
        return pow(dvec.X * dvec.X + dvec.Y * dvec.Y + dvec.Z * dvec.Z, 0.5);
    }

    void Collide() {
        Vec3 v1, v2;

        v1 = (- secondMass->mass / (firstMass->vel - secondMass->vel) + firstMass->vel*firstMass->mass + secondMass->vel*secondMass->mass) / (firstMass->mass + secondMass->mass);

        v2 = v1 + (1 / (firstMass->vel - secondMass->vel));

        firstMass->vel = v1;
        secondMass->vel = v2;
    }

    MSPair(Mass* mass1, Mass* mass2, Spring* spring1) : firstMass(mass1), secondMass(mass2), spring(spring1) {};
};

struct PhysTriangle {
    MSPair* Vert1;
    MSPair* Vert2;
    MSPair* Vert3;

    PhysTriangle(MSPair* a, MSPair* b, MSPair* c) {
        Vert1 = a;
        Vert2 = b;
        Vert3 = c;
    }

    std::vector<Mass*> GetMasses() const {
        std::vector<Mass*> masses;

        // Function to add a Mass* to the vector if it's not already present
        auto addMassIfNotPresent = [&](Mass* mass) {
            if (mass && std::find(masses.begin(), masses.end(), mass) == masses.end()) {
                masses.push_back(mass);
            }
        };

        addMassIfNotPresent(Vert1->firstMass);
        addMassIfNotPresent(Vert1->secondMass);
        addMassIfNotPresent(Vert2->firstMass);
        addMassIfNotPresent(Vert2->secondMass);
        addMassIfNotPresent(Vert3->firstMass);
        addMassIfNotPresent(Vert3->secondMass);

        return masses;
    }

    bool LineSegmentIntersect(const Vec3& p1, const Vec3& q1, const Vec3& p2, const Vec3& q2) {
        Vec3 u = q1 - p1;
        Vec3 v = q2 - p2;
        Vec3 w = p1 - p2;

        Vec3 n = cross(u, v);

        float denom = dot(n, n);

        // Check if the line segments are parallel (denom is close to zero)
        if (denom < EPSILON) {
            return false;
        }

        float sI = dot(cross(v, w), n) / denom;
        float tI = dot(cross(u, w), n) / denom;

        // Check if the intersection point is within the line segments
        if (sI >= 0.0 && sI <= 1.0 && tI >= 0.0 && tI <= 1.0) {
            return true; // Line segments intersect
        }

        return false; // No intersection found
    }

    bool IsCollidingWith(const PhysTriangle& other) {
        Vec3 vertices1[] = { Vert1->firstMass->pos, Vert2->firstMass->pos, Vert3->firstMass->pos };
        Vec3 vertices2[] = { other.Vert1->firstMass->pos, other.Vert2->firstMass->pos, other.Vert3->firstMass->pos };

        // Check for edge-edge intersection between the triangles
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                int next_i = (i + 1) % 3;
                int next_j = (j + 1) % 3;
                if (LineSegmentIntersect(vertices1[i], vertices1[next_i], vertices2[j], vertices2[next_j])) {
                    return true; // Triangles intersect
                }
            }
        }

        return false; // No intersection found
    }
};

struct HitBox {
    Vec3 minPoint;
    Vec3 maxPoint;

    HitBox() : minPoint(Vec3(0, 0, 0)), maxPoint(Vec3(0, 0, 0)) {}

    HitBox(const Vec3& min, const Vec3& max) : minPoint(min), maxPoint(max) {}

    bool ContainsPoint(const Vec3& point) const {
        return point.X >= minPoint.X && point.X <= maxPoint.X &&
            point.Y >= minPoint.Y && point.Y <= maxPoint.Y &&
            point.Z >= minPoint.Z && point.Z <= maxPoint.Z;
    }

    bool Intersects(const HitBox& other) const {
        return !(maxPoint.X < other.minPoint.X || minPoint.X > other.maxPoint.X ||
            maxPoint.Y < other.minPoint.Y || minPoint.Y > other.maxPoint.Y ||
            maxPoint.Z < other.minPoint.Z || minPoint.Z > other.maxPoint.Z);
    }
};

struct PhysModel {
    std::vector<PhysTriangle> triangles;
    std::vector<Mass*> masses;

    PhysModel(const std::vector<PhysTriangle>& Triangles) : triangles(Triangles) {
        FindUniqueMasses();
    }

    void FindUniqueMasses() {
        for (const PhysTriangle& triangle : triangles) {
            AddUniqueMasses(triangle.GetMasses());
        }
    }

    void AddUniqueMasses(const std::vector<Mass*>& Masses) {
        for (Mass* mass : Masses) {
            if (std::find(masses.begin(), masses.end(), mass) == masses.end()) {
                masses.push_back(mass);
            }
        }
    }
};

struct Triangle {
    Vec3 vert1;
    Vec3 vert2;
    Vec3 vert3;

    Triangle(const Vec3& vertex1, const Vec3& vertex2, const Vec3& vertex3) : vert1(vertex1), vert2(vertex2), vert3(vertex3) {};
};

struct Mesh {
    Mesh(const std::vector<Triangle>& Triangles) : triangles(Triangles) {}

    std::vector<Triangle> triangles;
};

struct PhysProp {
    PhysProp(PhysModel Model, HitBox Hitbox, Mass* Position) {};

    PhysModel model;
    HitBox hitbox;
    Vec3 position;
    Vec3 velocity;
    Vec3 facing;
};

struct Prop {
    Mesh mesh;
    HitBox hitbox;
    Vec3 position;
    Vec3 velocity;
    Vec3 facing;
    Vec3 mass;
    Vec3 massCentre;
};

class Simulation {
public:
    Simulation(std::vector<MSPair*> mspairs, double gravity = 9.81) {
        this->gravity = gravity;
        this->pairs = mspairs;

        for (size_t i = 0; i < mspairs.size(); i++) {
            MSPair* pair = mspairs[i];
            size_t massesLength = masses.size();

            bool IsAIn = false;
            bool IsBIn = false;
            for (size_t j = 0; j < massesLength; j++) {
                if (masses[j] == pair->firstMass) {
                    IsAIn = true;
                }

                if (masses[j] == pair->secondMass) {
                    IsBIn = true;
                }
            }

            if (!IsAIn) {
                masses.push_back(pair->firstMass);
            }

            if (!IsBIn) {
                masses.push_back(pair->secondMass);
            }
        }
    };

    void SolveForces(double dt) {
        MSPair* pair;
        Vec3 ForceAB;

        for (size_t i = 0; i < pairs.size(); i++) {
            pair = pairs[i];

            ForceAB = pair->spring->calculateForce(pair->secondMass->pos - pair->firstMass->pos, pair->Distance());
            
            pair->firstMass->applyForce( ForceAB, dt);
            pair->secondMass->applyForce(-ForceAB, dt);
        }

        for (size_t j = 0; j < masses.size(); j++) {
            masses[j]->Move(dt);
        }
    };

    MSPair* operator[](int index) {
        return pairs[index];
    }

private:
    double gravity;
    std::vector<MSPair*> pairs;
    std::vector<Mass*> masses;
};

int main() {
    // Create masses, springs, and MSPairs
    Mass mass1(2.0, Vec3(0, 0, 0));
    Mass mass2(1.0, Vec3(6, 0, 0));
    Mass mass3(1.0, Vec3(3, 2, 0));
    Mass mass4(1.0, Vec3(9, 2, 0));
    Mass mass5(1.0, Vec3(1.5, 4, 0));
    Mass mass6(1.0, Vec3(7.5, 4, 0));

    Spring spring1(5.0, 2.0);
    Spring spring2(5.0, 2.0);
    Spring spring3(5.0, 2.0);
    Spring spring4(5.0, 2.0);
    Spring spring5(5.0, 2.0);
    Spring spring6(5.0, 2.0);

    MSPair pair1(&mass1, &mass2, &spring1);
    MSPair pair2(&mass1, &mass3, &spring2);
    MSPair pair3(&mass2, &mass3, &spring3);
    MSPair pair4(&mass4, &mass5, &spring4);
    MSPair pair5(&mass4, &mass6, &spring5);
    MSPair pair6(&mass5, &mass6, &spring6);

    PhysTriangle tri1(&pair1, &pair2, &pair3);
    PhysTriangle tri2(&pair4, &pair5, &pair6);

    Simulation sim({&pair1, &pair2, &pair3, &pair4, &pair5, &pair6});

    // Simulate the system for 10 time steps
    for (int i = 0; i < 1000; i++) {
        sim.SolveForces(0.1);
    }

    return 0;
}


////////////////////////////////////////////////////////////////

//#include <iostream>
//#include <vector>
//#include <cmath>
//#include <memory>
//
//struct Vec3 {
//    double X, Y, Z;
//
//    Vec3(double x, double y, double z) : X(x), Y(y), Z(z) {};
//
//    Vec3() : X(0), Y(0), Z(0) {};
//
//    double Magnitude() const {
//        return std::sqrt(X * X + Y * Y + Z * Z);
//    }
//
//    Vec3& operator+=(const Vec3& other) {
//        X += other.X;
//        Y += other.Y;
//        Z += other.Z;
//        return *this;
//    }
//
//    Vec3& operator=(const Vec3& other) {
//        if (this != &other) {
//            X = other.X;
//            Y = other.Y;
//            Z = other.Z;
//        }
//        return *this;
//    }
//
//    Vec3& operator+=(double other) {
//        X += other;
//        Y += other;
//        Z += other;
//        return *this;
//    }
//
//    Vec3 operator-(const Vec3& other) const {
//        return Vec3(X - other.X, Y - other.Y, Z - other.Z);
//    }
//
//    Vec3 operator*(double num) const {
//        return Vec3(X * num, Y * num, Z * num);
//    }
//
//    Vec3 operator/(double num) const {
//        return Vec3(X / num, Y / num, Z / num);
//    }
//
//    Vec3 operator-() const {
//        return Vec3(-X, -Y, -Z);
//    }
//};
//
//class Mass {
//public:
//    Mass(double mass, Vec3 initial_position, Vec3 initial_velocity = Vec3(0, 0, 0)) : mass(mass), pos(initial_position), vel(initial_velocity) {}
//
//    void applyForce(const Vec3& force, double dt) {
//        vel += (force / mass) * dt;
//    }
//
//    void Move(double dt) {
//        pos += vel * dt;
//    }
//
//    Vec3 getPosition() const {
//        return pos;
//    }
//
//    double getMass() const {
//        return mass;
//    }
//
//private:
//    double mass;
//    Vec3 pos;
//    Vec3 vel;
//};
//
//class Spring {
//public:
//    Spring(double spring_constant, double equilibrium_position) : k(spring_constant), l0(equilibrium_position) {}
//
//    Vec3 calculateForce(const Vec3& dir, double distance) const {
//        double F = -k * (distance - l0);
//        return dir * (1.0 / F);
//    }
//
//private:
//    double k;
//    double l0;
//};
//
//// Masses, Spring pair
//struct MSPair {
//    Mass* A;
//    Mass* B;
//    Spring* spring;
//
//    double Distance() const {
//        Vec3 dvec = B->getPosition() - A->getPosition();
//        return dvec.Magnitude();
//    }
//
//    MSPair(Mass* mass1, Mass* mass2, Spring* spring1) : A(mass1), B(mass2), spring(spring1) {}
//};
//
//class Simulation {
//public:
//    Simulation(std::vector<MSPair*>& mspairs, double gravity = 9.81) : gravity(gravity), pairs(mspairs) {
//        for (auto pair : pairs) {
//            bool IsAIn = false;
//            bool IsBIn = false;
//            for (auto mass : masses) {
//                if (mass == pair->A) {
//                    IsAIn = true;
//                }
//                if (mass == pair->B) {
//                    IsBIn = true;
//                }
//            }
//
//            if (!IsAIn) {
//                masses.push_back(pair->A);
//            }
//
//            if (!IsBIn) {
//                masses.push_back(pair->B);
//            }
//        }
//    }
//
//    void SolveForces(double dt) {
//        for (auto pair : pairs) {
//            Vec3 ForceAB = pair->spring->calculateForce(pair->B->getPosition() - pair->A->getPosition(), pair->Distance());
//            pair->A->applyForce(ForceAB, dt);
//            pair->B->applyForce(-ForceAB, dt);
//        }
//
//        for (auto mass : masses) {
//            mass->Move(dt);
//        }
//    }
//
//    MSPair* operator[](int index) {
//        return pairs[index];
//    }
//
//private:
//    double gravity;
//    std::vector<MSPair*> pairs;
//    std::vector<Mass*> masses;
//};
//
//int main() {
//    std::unique_ptr<Mass> mass1 = std::make_unique<Mass>(2.0, Vec3(0, 0, 0));
//    std::unique_ptr<Mass> mass2 = std::make_unique<Mass>(1.0, Vec3(6, 0, 0));
//    std::unique_ptr<Spring> spring = std::make_unique<Spring>(5.0, 2.0);
//
//    MSPair* pair = new MSPair(mass1.get(), mass2.get(), spring.get());
//
//    std::vector<MSPair*> pairs({ pair });
//
//    Simulation sim(pairs);
//
//    MSPair* pair1 = sim[0];
//    Mass* A = pair1->A;
//    Mass* B = pair1->B;
//
//    for (int i = 0; i < 10; i++) {
//        sim.SolveForces(0.1);
//        std::cout << "Centre of Mass: " << (A->getMass() * A->getPosition().X + B->getMass() * B->getPosition().X) / (A->getMass() + B->getMass()) << std::endl;
//    }
//
//    return 0;
//}
