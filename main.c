#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <mpi.h>

typedef struct _body 
{
    float x , y ; // p o s i t i o n
    float ax , ay ; // a c c e l e r a t i o n
    float vx , vy ; // v e l o c i t y
    float mass ; 
} body ;

static const float DELTA_TIME = 0.2;
static const int NUMBER_OF_BODIES = 50;

void integrate(body *body, float DELTA_TIME );
void GenerateDebugData (int planetCount);
void SimulateWithBruteforce ();
Vector2 CalculateNewtonGravityAcceleration (body *a, body *b, float *ax, float *ay);
body bodies[NUMBER_OF_BODIES];

int main(int argc, char **argv) {
    
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for(size_t i = 0 ; i < NUMBER_OF_BODIES; ++i ) 
    {
        float total_ax = 0.0f, total_ay = 0.0f;
        for(size_t j = 0 ; j < N; ++j ) 
        {
            if( i == j ) { continue ; }
            float ax, ay;
            calculate_newton_gravity_acceleration(&bodies[i], &bodies[j], &ax , &ay);
            total_ax += ax;
            total_ay += ay;
        }
        bodies[i].ax = total_ax;
        bodies[i].ay = total_ay;
    }

    for(size_t i = 0 ; i < N; ++i )
    {
        integrate(&bodies[i],DT);
    }
    
    MPI_Finalize();

    return 0;
}

void SimulateWithBruteforce () {
    if (!IntegrateMovement) return;

    for (size_t i = 0; i < NUMBER_OF_BODIES; ++i) 
    {
        body planet = bodies[i];
        if (!planet.IsAlive) continue;

        Vector2 acceleration = Vector2.zero;
        foreach (PlanetController anotherPlanet in planets) 
        {
            if (planet == anotherPlanet || !anotherPlanet.IsAlive) continue;
            acceleration += CalculateNewtonGravityAcceleration (planet, anotherPlanet);
        }

        planet.Acceleration = acceleration;
    }
}

Vector2 CalculateNewtonGravityAcceleration (body *a, body *b, float *ax, float *ay) {
    ++interactions;
    Vector2 acceleration = Vector2.zero;

    Vector2 galacticPlaneR = secondBody.Position - firstBody.Position;

    float distanceSquared = galacticPlaneR.sqrMagnitude + simulationSofteningLengthSquared;
    float distanceSquaredCubed = distanceSquared * distanceSquared * distanceSquared;
    float inverse = 1.0f / Mathf.Sqrt (distanceSquaredCubed);
    float scale = secondBody.Mass * inverse;
    acceleration += galacticPlaneR * scale;

    return acceleration;
}


void integrate(body *body) {
    body−>vx += body−>ax * DELTA_TIME ;
    body−>vy += body−>ay * DELTA_TIME ;
    body−>x += body−>vx * DELTA_TIME ;
    body−>y += body−>vy * DELTA_TIME ;
}
void GenerateDebugData (int planetCount) {
    const float accelerationScale = 100.0f;
    const float galacticPlaneY = 0.0f;
    
    for (int i = 0; i < planetCount; ++i) {
        float angle =
            ((float) i / planetCount) * 2.0f * Mathf.PI +
            ((Random.value - 0.5f) * 0.5f);

        body planet = Instantiate (
                    PlanetTemplate,
                    new Vector3 (
                        Random.value,
                        galacticPlaneY,
                        Random.value
                    ),
                    Quaternion.identity
                );

            planet.GetComponent<Rigidbody> ().velocity =
                new Vector3 (
                    Mathf.Cos (angle),
                    0.0f,
                    Mathf.Sin (angle)
                ) * accelerationScale * Random.value;

            float initialMass =
                planet.Mass;
            planet.Mass =
                initialMass * Random.value + initialMass * 0.5f;

            float scale =
                (planet.Mass / (initialMass * 1.5f)) + 0.1f;
            planet.transform.localScale =
                new Vector3 (scale, scale, scale);

            SimulatorInstance.AddPlanet(planet);
        }
    }

