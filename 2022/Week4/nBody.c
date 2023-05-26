#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Small epsilon
#define EPS 0.0000001

/**
 * @brief Macro for easy termination in case of errors. Usage is similar to printf:
 * FATAL("This prints the number 5 and exits with an error code: %d\n", 5);
 */
#define FATAL(fmt, ...)                                                                                              \
    do {                                                                                                             \
        fprintf(stderr, "Fatal error: %s in %s() on line %d:\n\t" fmt, __FILE__, __func__, __LINE__, ##__VA_ARGS__); \
        exit(EXIT_FAILURE);                                                                                          \
    } while (0)

/**
 * @brief 3D vector.
 */
typedef struct vec3 {
    float x;
    float y;
    float z;
} vec3;

/**
 * @brief A particle consisting of a position, velocity and mass.
 */
typedef struct Particle {
    vec3 pos;
    vec3 v;
    float mass;
} Particle;

/**
 * @brief Reads input particles from a given file. The first line of the file is the number of particles N. The next N
 * lines each contains 7 floats: the first 3 are the initial position of the particle, the next 3 are the current
 * velocity of the particle and the last float is the mass of the particle.
 *
 * @param inputFile Input file.
 * @param numParticles This pointer will be updated with the number of particles found in the file.
 * @return Particle* Array containing the particles from the file.
 */
Particle *readInput(FILE *inputFile, int *numParticles) {
    int n;
    fscanf(inputFile, "%d\n", &n);
    Particle *particles = malloc(n * sizeof(Particle));
    for (int i = 0; i < n; i++) {
        fscanf(inputFile, "%f %f %f ", &particles[i].pos.x, &particles[i].pos.y, &particles[i].pos.z);
        fscanf(inputFile, "%f %f %f ", &particles[i].v.x, &particles[i].v.y, &particles[i].v.z);
        fscanf(inputFile, "%f\n", &particles[i].mass);
    }
    *numParticles = n;
    return particles;
}

/**
 * @brief Saves the state of all particles to the provided file.
 *
 * @param file File to save the particle data to.
 * @param particles The particles to save.
 * @param numParticles Number of particles to save.
 */
void saveParticles(FILE *file, Particle *particles, int numParticles) {
    fprintf(file, "%d\n", numParticles);
    for (int i = 0; i < numParticles; i++) {
        fprintf(file, "%.1f %.1f %.1f ", particles[i].pos.x, particles[i].pos.y, particles[i].pos.z);
        fprintf(file, "%.1f %.1f %.1f ", particles[i].v.x, particles[i].v.y, particles[i].v.z);
        fprintf(file, "%.1f\n", particles[i].mass);
    }
}

/**
 * @brief Performs a serial N-Body simulation with the given particles for the provided number of time-steps.
 *
 * @param particles The particles to simulate.
 * @param n The total number of particles.
 * @param timeSteps The number of time-steps to simulate for.
 * @param numThreads Number of threads to use for the simulation. Currently unused.
 */
void startSimulation(Particle *particles, int n, int timeSteps, int numThreads) {
    // Accelerations
    vec3 *acc = malloc(n * sizeof(vec3));

    for (int t = 0; t < timeSteps; t++) {
        memset(acc, 0, n * sizeof(vec3));
        // Compute all-to-all forces and accelerations
        for (int q = 0; q < n; q++) {
            for (int j = 0; j < n; j++) {
                if (q == j) {
                    // Skip interaction with self
                    continue;
                }
                float rx = particles[j].pos.x - particles[q].pos.x;
                float ry = particles[j].pos.y - particles[q].pos.y;
                float rz = particles[j].pos.z - particles[q].pos.z;
                float dd = rx * rx + ry * ry + rz * rz + EPS;
                float d = 1 / sqrtf(dd * dd * dd);
                float s = particles[j].mass * d;
                acc[q].x += rx * s;
                acc[q].y += ry * s;
                acc[q].z += rz * s;
            }
        }
        // Update positions and velocities
        for (int i = 0; i < n; i++) {
            particles[i].pos.x += particles[i].v.x;
            particles[i].pos.y += particles[i].v.y;
            particles[i].pos.z += particles[i].v.z;
            particles[i].v.x += acc[i].x;
            particles[i].v.y += acc[i].y;
            particles[i].v.z += acc[i].z;
        }
    }
    free(acc);
}

/**
 * @brief N-Body simulation. Needs at least arguments: a number of time-steps to simulate and a number of
 * threads. Optionally, an input file can be given as final argument.
 * Usage: ./a.out <time-steps> <num-threads> [inputFile]
 * Example usage: ./a.out 5 2
 * This will read the configuration from stdin and perform the simulation for 5 time steps. In the serial
 * version, the number of threads argument is unused.
 * Basic Compilation: gcc nBody.c -lm
 * When evaluating performance: compile with additional -O3 flag
 *
 * @param argc Argument count.
 * @param argv Arguments. Two arguments are expected: a number of time-step and number of threads (in
 * that order). A third argument - an input file - is optional.
 * @return int Exit code.
 */
int main(int argc, char *argv[]) {
    if (argc != 3 && argc != 4) {
        FATAL(
            "Please provide the number of time steps to simulate and number of threads to use.\n"
            "Optionally, provide an input file to read from.\n");
    }
    int timeSteps = atoi(argv[1]);
    int numThreads = atoi(argv[2]);
    char *inputFilePath = argv[3];

    int numParticles;
    fprintf(stderr, "Reading input...\n\n");
    Particle *particles;

    if (inputFilePath == NULL) {
        particles = readInput(stdin, &numParticles);
    } else {
        FILE *inputFile = fopen(inputFilePath, "r");
        if (inputFile == NULL) {
            FATAL("Cannot open file %s\n", inputFilePath);
        }
        particles = readInput(inputFile, &numParticles);
        fclose(inputFile);
    }

    fprintf(stderr, "Starting n-body simulation using %d threads\n", numThreads);
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    startSimulation(particles, numParticles, timeSteps, numThreads);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double timeSpent = (end.tv_sec - start.tv_sec);
    timeSpent += (end.tv_nsec - start.tv_nsec) / 1e9;
    fprintf(stderr, "Simulation of %d particles complete -- Took: %f seconds.\n\n", numParticles, timeSpent);

    fprintf(stderr, "Saving particles...\n");
    saveParticles(stdout, particles, numParticles);

    free(particles);
    return EXIT_SUCCESS;
}