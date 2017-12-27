//Use it only to generate planets using method GenerateDebugData

typedef struct _body {
float x , y ; // p o s i t i o n
float ax , ay ; // a c c e l e r a t i o n
float vx , vy ; // v e l o c i t y
float mass ; 
} body ;
static const delta_time;
static const int number_of_planets = 10;

body bodies[number_of_planets];

    private float galacticPlaneY =
        0.0f;
    private void AddPlanet (Vector3 position) {
        PlanetController planet =
            Instantiate(
                PlanetTemplate,
                position,
                Quaternion.identity
            );

        SimulatorInstance.AddPlanet (planet);
    }

    private void GenerateDebugData () {
        const int planetCount =
            180;
        const float accelerationScale =
            100.0f;

        for (int i = 0; i < planetCount; ++i) {
            float angle =
                ((float) i / planetCount) * 2.0f * Mathf.PI +
                    ((Random.value - 0.5f) * 0.5f);

            PlanetController planet =
                Instantiate (
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

