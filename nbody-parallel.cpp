// Evan Stark - April 13th 2025 - ITCS 4145 001
// This program provides a parallel solution to
// the N-Body problem: a program that simulates the movement of
// celestial bodies/particles based on the variables of all other 
// bodies in the system.

// SOURCES USED:
// Starter code provided by Erik Saule.
// Additional help prevoided by TA Devashish Bhat.
// https://www.youtube.com/watch?v=29mF-kqNWpE (Video going over the basics of using OpenMP)
// https://hpc-tutorials.llnl.gov/openmp/ (Written tutorial of using OpenMP)
// https://www.youtube.com/watch?v=fiMRQSE-Ak8 (Video going over synchronization in OpenMP)

// All necessary libraries
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <vector>
#include <omp.h>

// Gravity.
double G = 6.674*std::pow(10,-11);

// Building a new simulation struct to hold a number of particles,
// initializing each one's mass, position, velocity, and force.
struct simulation 
{
    size_t nbpart;
    
    std::vector<double> mass;

    //position
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    //velocity
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;

    //force
    std::vector<double> fx;
    std::vector<double> fy;
    std::vector<double> fz;

    // Constructing a new simulation object.
    simulation(size_t nb):
        nbpart(nb), mass(nb),
        x(nb), y(nb), z(nb),
        vx(nb), vy(nb), vz(nb),
        fx(nb), fy(nb), fz(nb) 
    {}
};

// Assigning random values to each particle.
void random_init(simulation& s) 
{
    // Initializing a random number generator.
    std::random_device rd;  
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dismass(0.9, 1.);
    std::normal_distribution<double> dispos(0., 1.);
    std::normal_distribution<double> disvel(0., 1.);
    
    // Assigning new values to each particle via the RNG.
    // Parallelize the random particle generation.
    #pragma omp parallel for
    for (size_t i = 0; i<s.nbpart; ++i) {
        s.mass[i] = dismass(gen);

        s.x[i] = dispos(gen);
        s.y[i] = dispos(gen);
        s.z[i] = dispos(gen);
        s.z[i] = 0.;
        
        s.vx[i] = disvel(gen);
        s.vy[i] = disvel(gen);
        s.vz[i] = disvel(gen);
        s.vz[i] = 0.;
        s.vx[i] = s.y[i]*1.5;
        s.vy[i] = -s.x[i]*1.5;
    }

    return;
  
    //normalize velocity (using normalization found on some physicis blog)
    double meanmass = 0;
    double meanmassvx = 0;
    double meanmassvy = 0;
    double meanmassvz = 0;

    // Applying normalization.
    for (size_t i = 0; i<s.nbpart; ++i) {
        meanmass += s.mass[i];
        meanmassvx += s.mass[i] * s.vx[i];
        meanmassvy += s.mass[i] * s.vy[i];
        meanmassvz += s.mass[i] * s.vz[i];
    }
    for (size_t i = 0; i<s.nbpart; ++i) {
        s.vx[i] -= meanmassvx/meanmass;
        s.vy[i] -= meanmassvy/meanmass;
        s.vz[i] -= meanmassvz/meanmass;
    }
    
}

// Initializing the solar system simulation test case.
void init_solar(simulation& s) 
{
    enum Planets {SUN, MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE, MOON};
    s = simulation(10);

    // Masses in kg
    s.mass[SUN] = 1.9891 * std::pow(10, 30);
    s.mass[MERCURY] = 3.285 * std::pow(10, 23);
    s.mass[VENUS] = 4.867 * std::pow(10, 24);
    s.mass[EARTH] = 5.972 * std::pow(10, 24);
    s.mass[MARS] = 6.39 * std::pow(10, 23);
    s.mass[JUPITER] = 1.898 * std::pow(10, 27);
    s.mass[SATURN] = 5.683 * std::pow(10, 26);
    s.mass[URANUS] = 8.681 * std::pow(10, 25);
    s.mass[NEPTUNE] = 1.024 * std::pow(10, 26);
    s.mass[MOON] = 7.342 * std::pow(10, 22);

    // Positions (in meters) and velocities (in m/s)
    double AU = 1.496 * std::pow(10, 11); // Astronomical Unit

    // Assigning position and velocity values.
    s.x = {0, 0.39*AU, 0.72*AU, 1.0*AU, 1.52*AU, 5.20*AU, 9.58*AU, 19.22*AU, 30.05*AU, 1.0*AU + 3.844*std::pow(10, 8)};
    s.y = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    s.z = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    s.vx = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    s.vy = {0, 47870, 35020, 29780, 24130, 13070, 9680, 6800, 5430, 29780 + 1022};
    s.vz = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
}

// Updates the force that particle from applies on particle to
void update_force(simulation& s, size_t from, size_t to) {
    double softening = .1;
    double dist_sq = std::pow(s.x[from]-s.x[to],2)
        + std::pow(s.y[from]-s.y[to],2)
        + std::pow(s.z[from]-s.z[to],2);
    double F = G * s.mass[from]*s.mass[to]/(dist_sq+softening); // Force.

    // Direction and normalizing it.
    double dx = s.x[from]-s.x[to];
    double dy = s.y[from]-s.y[to];
    double dz = s.z[from]-s.z[to];
    double norm = std::sqrt(dx*dx+dy*dy+dz*dz);
    
    // Normalization in each direction.
    dx = dx/norm;
    dy = dy/norm;
    dz = dz/norm;

    // Apply forces
    s.fx[to] += dx*F;
    s.fy[to] += dy*F;
    s.fz[to] += dz*F;
}

// Resetting forces before each time step.
void reset_force(simulation& s) {
    for (size_t i=0; i<s.nbpart; ++i) {
        s.fx[i] = 0.;
        s.fy[i] = 0.;
        s.fz[i] = 0.;
    }
}

// Applying forces and updating the velocities in each direction.
void apply_force(simulation& s, size_t i, double dt) {
    s.vx[i] += s.fx[i]/s.mass[i]*dt;
    s.vy[i] += s.fy[i]/s.mass[i]*dt;
    s.vz[i] += s.fz[i]/s.mass[i]*dt;
}

// Update particle position using new velocity value.
void update_position(simulation& s, size_t i, double dt) {
    s.x[i] += s.vx[i]*dt;
    s.y[i] += s.vy[i]*dt;
    s.z[i] += s.vz[i]*dt;
}

// Dump the status of the state.
// Each step = one line.
void dump_state(simulation& s) {
    std::cout<<s.nbpart<<'\t';
    for (size_t i=0; i<s.nbpart; ++i) {
        std::cout<<s.mass[i]<<'\t';
        std::cout<<s.x[i]<<'\t'<<s.y[i]<<'\t'<<s.z[i]<<'\t';
        std::cout<<s.vx[i]<<'\t'<<s.vy[i]<<'\t'<<s.vz[i]<<'\t';
        std::cout<<s.fx[i]<<'\t'<<s.fy[i]<<'\t'<<s.fz[i]<<'\t';
    }
    std::cout<<'\n';
}

// Using this function to load particle information from an outside file.
void load_from_file(simulation& s, std::string filename) {
    std::ifstream in (filename);
    size_t nbpart;
    in>>nbpart;
    s = simulation(nbpart);
    for (size_t i=0; i<s.nbpart; ++i) {
        in>>s.mass[i];
        in >>  s.x[i] >>  s.y[i] >>  s.z[i];
        in >> s.vx[i] >> s.vy[i] >> s.vz[i];
        in >> s.fx[i] >> s.fy[i] >> s.fz[i];
    }
    
    // Throw exception if file cannot be read.
    if (!in.good())
        throw "kaboom";
}

// Updating force by comparing particle with every other particle. 
// Will be called by each thread after OpenMP parallel pragma.
void threaded_update_force(simulation& s, size_t beg, size_t end)
{
    for (size_t i = beg; i < end; i++) 
    {
        for (size_t j = 0; j < s.nbpart; j++) 
        {
            if (i != j) 
            {
                update_force(s, i, j);
            }
        }
    }
}

// Main driver function.
int main(int argc, char* argv[]) {
    // Check if the user has 5 arguments passed in.
    if (argc != 5) {
        std::cerr
        <<"usage: "<<argv[0]<<" <input> <dt> <nbstep> <printevery>"<<"\n"
        <<"input can be:"<<"\n"
        <<"a number (random initialization)"<<"\n"
        <<"planet (initialize with solar system)"<<"\n"
        <<"a filename (load from file in singleline tsv)"<<"\n";
        return -1;
    }
    
    // dt = how large the time step is.
    // nbstep = number of time steps to run.
    // printevery = how frequent to dump the current state.
    double dt = std::atof(argv[2]); //in seconds
    size_t nbstep = std::atol(argv[3]);
    size_t printevery = std::atol(argv[4]);
  
    // Initialize simulation with 3 particles.
    simulation s(3);

    //parse command line
    {
        size_t nbpart = std::atol(argv[1]); //return 0 if not a number
        if ( nbpart > 0) 
        {
            s = simulation(nbpart);
            random_init(s);
        } 
        
        else 
        {
            std::string inputparam = argv[1];
            
            if (inputparam == "planet") 
            {
                init_solar(s);
            } 
            
            else 
            {
                load_from_file(s, inputparam);
            }
        }    
    }

    // Calculating how many threads are needed.
    // Allocating ~100 particles per thread.
    int num_threads = s.nbpart / 100;
    
    // Need to call at least one thread!
    if (num_threads == 0) {
        num_threads = 1;
    }

    for (size_t step = 0; step< nbstep; step++) 
    {
        // Printing out dump state once every one in a while.
        if (step %printevery == 0)
            dump_state(s);

        // Resetting force before making calcs.
        reset_force(s);

        // Parallelizing force updates.
        #pragma omp parallel for num_threads(num_threads)
        for (size_t i = 0; i < num_threads; i++) 
        {
            size_t beg = (s.nbpart * i) / num_threads;  // region for thread to start at.
            size_t end;
            if (i = num_threads - 1)
            {
                end = s.nbpart;   // putting all leftover particles in the last thread.
            }

            else {
                end = (i+1 * s.nbpart) / num_threads;   // ending at a certain range otherwise.
            }

            threaded_update_force(s, beg, end);
        }

        // Apply forces and update each particle's position in parallel.
        #pragma omp parallel for num_threads(num_threads)
        for (size_t i=0; i<s.nbpart; ++i) 
        {
            apply_force(s, i, dt);
            update_position(s, i, dt);
        }
    }

    return 0;
}