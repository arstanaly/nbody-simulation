#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <mpi.h>

//Use it only to generate planets using method GenerateDebugData

typedef struct _body {
float x , y ; // p o s i t i o n
float ax , ay ; // a c c e l e r a t i o n
float vx , vy ; // v e l o c i t y
float mass ; 
} body ;
static const delta_time;
static const int NUMBER_OF_BODIES = 50;

body *bodies = (body *) malloc(sizeof(*bodies) * NUMBER_OF_BODIES);

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


int main(int argc, char **argv) {
    
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    
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

Vector2 CalculateNewtonGravityAcceleration (body *a, body *b, float *ax, float *ay) 
{
        ++interactions;
        if (ShowInteractions)
            DrawDebugLines (firstBody, secondBody);

        Vector2 acceleration =
            Vector2.zero;

        Vector2 galacticPlaneR =
            secondBody.Position - firstBody.Position;

        float distanceSquared =
            galacticPlaneR.sqrMagnitude + simulationSofteningLengthSquared;
        float distanceSquaredCubed =
            distanceSquared * distanceSquared * distanceSquared;
        float inverse =
            1.0f / Mathf.Sqrt (distanceSquaredCubed);
        float scale =
            secondBody.Mass * inverse;

        acceleration +=
            galacticPlaneR * scale;

        return acceleration;
    }

void integrate(body *body, float delta_time )
{
    body−>vx += body−>ax * delta_time ;
    body−>vy += body−>ay * delta_time ;
    body−>x += body−>vx * delta_time ;
    body−>y += body−>vy * delta_time ;
}