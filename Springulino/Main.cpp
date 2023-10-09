#include <iostream>
#include <vector>

class Mass {
public:
    Mass(double mass, double initial_position, double initial_velocity) {
        this->mass = mass;
        this->pos = initial_position;
        this->vel = initial_velocity;
        this->accel = 0;
    }

    void applyForce(double force, double dt) {
        vel += (force / mass) * dt;
    }

    void Move(double dt) {
        pos += vel * dt;
    }

    double getPosition() {
        return pos;
    }

    double getMass() {
        return mass;
    }

    double mass;
    double pos;
    double vel;
    double accel;
};

class Spring {
public:
    Spring(double spring_constant, double equilibrium_position) {
        this->k = spring_constant;
        this->l0 = equilibrium_position;
    }

    double calculateForce(double position1, double position2) {
        return -k * (position1 - position2 - l0);
    }

    double k;
    double l0;
};

// Masses, Spring pair
struct MSPair {
    Mass* A;
    Mass* B;
    Spring* spring;

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
            for (size_t i = 0; i < massesLength; i++) {
                if (masses[0] == pair->A) {
                    IsAIn = true;
                }

                if (masses[0] == pair->B) {
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

            double ForceAB = pair->spring->calculateForce(pair->A->getPosition(), pair->B->getPosition());
            
            pair->A->applyForce( ForceAB, dt);
            pair->B->applyForce(-ForceAB, dt);
        }

        for (size_t i = 0; i < masses.size(); i++) {
            masses[i]->Move(dt);
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
    Mass* mass1 = new Mass(2.0, 0.0, 0.0);
    Mass* mass2 = new Mass(1.0, 6.0, 0.0);
    Spring* spring = new Spring(5.0, 2.0);

    MSPair* pair = new MSPair(mass1, mass2, spring);

    std::vector<MSPair*> pairs({ pair });

    Simulation sim(pairs);

    MSPair* pair1 = sim[0];
    Mass* A = pair1->A;
    Mass* B = pair1->B;

    for (int i = 0; i < 1000; i++) {
        sim.SolveForces(0.1);
        std::cout << "Centre of Mass: " << (A->mass * A->pos + B->mass * B->pos) / (A->mass + B->mass) << std::endl;
    }

    return 0;
}