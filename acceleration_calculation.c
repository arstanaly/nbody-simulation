private void SimulateWithBruteforce () {
        if (!IntegrateMovement)
            return;

        foreach (PlanetController planet in planets) {
            if (!planet.IsAlive)
                continue;

            Vector2 acceleration = Vector2.zero;
            foreach (PlanetController anotherPlanet in planets) {
                if (planet == anotherPlanet || !anotherPlanet.IsAlive)
                    continue;

                acceleration +=
                    CalculateNewtonGravityAcceleration (
                        planet, anotherPlanet
                    );
            }

            planet.Acceleration =
                acceleration;
        }
    }

private Vector2 CalculateNewtonGravityAcceleration (
                        IBody firstBody,
                        IBody secondBody
                    ) {
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