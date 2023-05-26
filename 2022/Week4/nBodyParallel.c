/* file: nBodyParallel.c */
/* author: Radu Rusu (email: rusu@student.rug.nl) */
/* date: 19 May 2022 */
/* version: 1.0 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

/* An example of barriers used for understanding the defining and initialization
https://gist.github.com/shelterz/4b13459668eec743f15be6c200aa91b2 */

pthread_barrier_t barrier;

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
 * @brief A parameters set consisting of start index, end index, number of particles, 
 * time steps and particles array.
 * 
 */
typedef struct Parameters {
    int start;
    int end;
    int n; 
    int timesteps;
    Particle *particles;
} Parameters;

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
 * @brief Performs a parallel N-Body simulation with the given particles for the provided number of time-steps.
 *
 * @param parameters A struct with: the particles to simulate, the total number of particles, the number of 
 * time-steps to simulate for and the bounds for each thread to use for the simulation. 
 */
void *calculateForceAccel(void *parameters) {
    /* Since the values of the struct are converted to void to be used in the function call, 
    we need to cast them back to Parameters type*/
    Parameters *p = (Parameters *)parameters;
    
    Particle *particles = p->particles;
    int start = p->start;
    int end = p->end;
    int numParticles = p->n;
    int timeSteps = p->timesteps;

    vec3 *acc = malloc(numParticles * sizeof(vec3));

    /* To avoid multiple thread creation for each timeStep, which takes a lot of time and performance, 
    each thread executes for the provided number of timeSteps */
    for (int t = 0; t < timeSteps; t++) {

        memset(acc, 0, numParticles * sizeof(vec3));

        /* The nested for-loops are parallelized by giving each thread a series of particles 
        (numParticles / numThreads) and it computes the interaction with all other ones exepting inself */
        for (int q = start; q < end; q++) {
            for (int j = 0; j < numParticles; j++) {
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
        
        /* The barrier is needed here to wait for all threads to finish calculating the 
        acceleration in the previous 2 loops, so that we get accurate results when 
        updating the position and velocity */
        pthread_barrier_wait(&barrier);

        /* Update positions and velocities */
        for (int i = start; i < end; i++) {
            particles[i].pos.x += particles[i].v.x;
            particles[i].pos.y += particles[i].v.y;
            particles[i].pos.z += particles[i].v.z;
            particles[i].v.x += acc[i].x;
            particles[i].v.y += acc[i].y;
            particles[i].v.z += acc[i].z;
        }

        /* The barrier is needed here to wait for all threads to finish updating the position 
        and velocity of the particles before freeing the array and returning */
        pthread_barrier_wait(&barrier);
    }

    free(acc);
    return NULL;
}

/**
 * @brief Starts the threads and organizes the necessary parameters for their execution.
 *
 * @param particles The particles to simulate.
 * @param n The total number of particles.
 * @param timeSteps The number of time-steps to simulate for.
 * @param numThreads Number of threads to use for the simulation. Currently unused.
 */
void startParallelSimulation(Particle *particles, int n, int timeSteps, int numThreads) {
    pthread_t threadID[numThreads];
    /* Each thread should have its parameters and they should be sent by value to the function that is 
    executed in parallel, otherwise multiple threads access the same memory reference and the results are incorrect. */
    Parameters parameters[numThreads];

    for(int i = 0; i < numThreads; i++) {

        /* These lines compute the array bounds for each thread, making sure that the case when numParticles 
        can not be divided without a remainder to the numThreads is handled*/
        int remainder = n % numThreads;
        int low = i * ((n - remainder) / numThreads);
        int high = (i + 1) * ((n - remainder) / numThreads);

        if (i == numThreads - 1)
            high += remainder;

        parameters[i].start = low;
        parameters[i].end = high;
        parameters[i].particles = particles;
        parameters[i].n = n;
        parameters[i].timesteps = timeSteps;

        pthread_create(&threadID[i], NULL, calculateForceAccel, &parameters[i]);
    }

    /* For-loop to ensure that all threads, which were created, have executed their part 
    and terminated with a successful exit code*/
    for(int i = 0; i < numThreads; i++){
        pthread_join(threadID[i], NULL);
    }
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

    /* After obtaining the number of threads, the barrier is initialized with the number and 
    pointer to the global variable */
    pthread_barrier_init(&barrier, NULL, numThreads);

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

    /* Starting the simulation */
    startParallelSimulation(particles, numParticles, timeSteps, numThreads);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double timeSpent = (end.tv_sec - start.tv_sec);
    timeSpent += (end.tv_nsec - start.tv_nsec) / 1e9;
    fprintf(stderr, "Simulation of %d particles complete -- Took: %f seconds.\n\n", numParticles, timeSpent);

    fprintf(stderr, "Saving particles...\n");
    saveParticles(stdout, particles, numParticles);

    free(particles);
    pthread_barrier_destroy(&barrier);
    return EXIT_SUCCESS;
}