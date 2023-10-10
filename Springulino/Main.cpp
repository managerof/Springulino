#include <iostream>
#include <vector>
#include <cmath>

struct Vec3 {
    double X, Y, Z;

    Vec3(double x, double y, double z) : X(x), Y(y), Z(z) {};

    Vec3() : X(0), Y(0), Z(0) {};

    double Magnitude() {
        return pow((X * X + Y * Y + Z * Z), 0.5);;
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

    Vec3& operator-(const Vec3& other) {
        Vec3 vec;
        
        vec.X = this->X - other.X;
        vec.Y = this->Y - other.Y;
        vec.Z = this->Z - other.Z;

        return vec;
    }

    Vec3 operator*(double num) const {
        return Vec3(X * num, Y * num, Z * num);
    }

    Vec3 operator/(double num) const {
        return Vec3(X / num, Y / num, Z / num);
    }

    Vec3& operator-(double num) {
        Vec3 newVec3(0, 0, 0);

        newVec3.X -= num;
        newVec3.Y -= num;
        newVec3.Z -= num;
        return newVec3;
    }

    Vec3 operator-() const {
        return Vec3(-X, -Y, -Z);
    }
};

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
    Mass* A;
    Mass* B;
    Spring* spring;

    double Distance() {
        Vec3 dvec = B->pos - A->pos;
        return pow(dvec.X * dvec.X + dvec.Y * dvec.Y + dvec.Z * dvec.Z, 0.5);
    }

    MSPair(Mass* mass1, Mass* mass2, Spring* spring1) : A(mass1), B(mass2), spring(spring1) {};
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
                if (masses[j] == pair->A) {
                    IsAIn = true;
                }

                if (masses[j] == pair->B) {
                    IsBIn = true;
                }
            }

            if (!IsAIn) {
                masses.push_back(pair->A);
            }

            if (!IsBIn) {
                masses.push_back(pair->B);
            }
        }
    };


    void SolveForces(double dt) {
        MSPair* pair;

        for (size_t i = 0; i < pairs.size(); i++) {
            pair = pairs[i];

            Vec3 ForceAB = pair->spring->calculateForce(pair->B->pos - pair->A->pos, pair->Distance());
            
            pair->A->applyForce( ForceAB, dt);
            pair->B->applyForce(-ForceAB, dt);
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
    Mass* mass1 = new Mass(2.0, Vec3(0, 0, 0));
    Mass* mass2 = new Mass(1.0, Vec3(6, 0, 0));
    Spring* spring = new Spring(5.0, 2.0);

    MSPair* pair = new MSPair(mass1, mass2, spring);

    std::vector<MSPair*> pairs({ pair });

    Simulation sim(pairs);

    MSPair* pair1 = sim[0];
    Mass* A = pair1->A;
    Mass* B = pair1->B;

    for (int i = 0; i < 10; i++) {
        sim.SolveForces(0.1);
        std::cout << "Centre of Mass: " << (A->mass * A->pos.X + B->mass * B->pos.X) / (A->mass + B->mass) << std::endl;
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
