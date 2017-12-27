#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

#define PI 3.14159265358979323846

typedef struct _body 
{
    float id;
    float x, y; //position
    float ax, ay; //acceleration
    float vx, vy; //velocity
    float mass;
} body;

static const float DELTA_TIME = 0.2;
static const int NUMBER_OF_BODIES = 100;
static const float SIMULATION_SOFTENING_LENTH = 100.0f;
static const float SIMULATION_SOFTENING_LENTH_SQUARED = 10000.0f;
static const float DURATION = 10.0f;

body *GenerateDebugData (int NUMBER_OF_BODIES);
void CalculateNewtonGravityAcceleration (body *a, body *b, float *ax, float *ay);
static void integratebodies();

int main(int argc, char **argv) {

    body *bodies = GenerateDebugData(NUMBER_OF_BODIES);

    printf("%d/n", NUMBER_OF_BODIES);
    printf("%f/n", DURATION);
    printf("%f/n", DELTA_TIME);
    for(int i = 0; i < NUMBER_OF_BODIES; ++i)
    {
        printf("Body #%f\n", bodies[i].id);
		printf("x = %f, y = %f\n", bodies[i].x, bodies[i].y);
		printf("ax = %f, ay = %f\n", bodies[i].ax, bodies[i].ay);
		printf("vx = %f, vy = %f\n", bodies[i].vx, bodies[i].vy);
		printf("mass = %f\n\n", bodies[i].mass);
    }

    MPI_Init(&argc, &argv);

	MPI_Datatype MPI_BODY;
	MPI_Type_contiguous(5, MPI_FLOAT, &MPI_BODY);
	MPI_Type_commit(&MPI_BODY);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *sendcounts = malloc(sizeof(int) * world_size);;    // array describing how many elements to send to each process
    int *displs = malloc(sizeof(body) * world_size);;        // array describing the displacements where each segment begins
    int rem = NUMBER_OF_BODIES % world_size;
    int sum = 0; 
    for (int i = 0; i < world_size; i++) {
        sendcounts[i] = NUMBER_OF_BODIES / world_size;
        if (rem > 0) {
            sendcounts[i]++;
            rem--;
        }
        displs[i] = sum;
        sum += sendcounts[i];
    }
    
    for (float time_active = 0; time_active < DURATION; time_active += DELTA_TIME) {
        body *local_bodies =  (body *) malloc(sizeof(*local_bodies) * sendcounts[rank]);

        MPI_Scatterv(bodies, sendcounts, displs, MPI_BODY, local_bodies, sendcounts[rank], MPI_BODY, 0, MPI_COMM_WORLD);

        for(size_t i = 0 ; i < sendcounts[rank]; ++i ) 
        {
            float total_ax = 0.0f, total_ay = 0.0f;
            for(size_t j = 0 ; j < NUMBER_OF_BODIES; ++j ) 
            {
                if( local_bodies[i].id == bodies[j].id ) { continue ; }
                float ax, ay;
                CalculateNewtonGravityAcceleration(&local_bodies[i], &local_bodies[j], &ax , &ay);
                total_ax += ax;
                total_ay += ay;
            }
            local_bodies[i].ax = total_ax;
            local_bodies[i].ay = total_ay;
        }

        MPI_Gatherv(&local_bodies, sendcounts[rank], MPI_BODY, bodies, sendcounts, displs, MPI_BODY, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
	    free(local_bodies);

        if(rank == 0)
        {
            integratebodies(bodies);
        }
    }
    MPI_Finalize();

    return 0;
}

static void integratebodies (body *bodies)
{
	for (size_t i = 0; i < NUMBER_OF_BODIES; ++i) {
		bodies[i].vx += bodies[i].ax * DELTA_TIME ;
		bodies[i].vy += bodies[i].ay * DELTA_TIME ;
		bodies[i].x += bodies[i].vx * DELTA_TIME ;
		bodies[i].y += bodies[i].vy * DELTA_TIME ;
        printf("%f  %f\n", bodies[i].ax, bodies[i].ay);
	}
}

body *GenerateDebugData (int NUMBER_OF_BODIES) {

    srand(time(NULL));

    const float accelerationScale = 100.0f;
    const float galacticPlaneY = 0.0f;

    body *result = (body *) malloc(sizeof(*result) * NUMBER_OF_BODIES);

    for (int i = 0; i < NUMBER_OF_BODIES; ++i)
    {
        float rnd_val = (float)rand() / (float)(RAND_MAX);
        float angle = ((float) i / NUMBER_OF_BODIES) * 2.0f * PI + ((rnd_val - 0.5f) * 0.5f);

        body *planet = (body *) malloc(sizeof(*planet));
        planet->id = i;
		planet->mass = (rand() + 0.5f) * 100000.0f;
		planet->ax = 0;
		planet->ay = 0;
		planet->vx = cos(angle) * accelerationScale * rand();
		planet->vy = sin(angle) * accelerationScale * rand();

		planet->x = (float)rand() / (float)(RAND_MAX);
		planet->y = galacticPlaneY;
		result[i] = *planet;
    }

    	return result;
}

void CalculateNewtonGravityAcceleration (body *a, body *b, float *ax, float *ay) {
    float galacticPlaneX = b->x - a->x;
	float galacticPlaneY = b->y - a->y;

	float distanceSquared = (galacticPlaneX * galacticPlaneX) + (galacticPlaneY * galacticPlaneY) + SIMULATION_SOFTENING_LENTH_SQUARED;

	float distanceSquaredCubed = distanceSquared * distanceSquared * distanceSquared;

	float inverse = 1.0f / sqrtf(distanceSquaredCubed);

	float scale = b->mass * inverse;

	*ax = galacticPlaneX * scale;
	*ay = galacticPlaneY * scale;
}
