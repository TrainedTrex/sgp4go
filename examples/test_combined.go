package main

import (
	"fmt"
	"time"

	"github.com/TrainedTrex/sgp4go"
)

func main() {
	// Example TLE for GEO satellite
	line1 := "1 28884U 05041A   25258.23425557 -.00000019  00000-0  00000-0 0  9996"
	line2 := "2 28884   2.8166  80.2739 0017303 131.7673 259.8922  0.98940294 72726"

	// Parse the TLE
	tle, err := sgp4go.ParseTLE(line1, line2)
	if err != nil {
		fmt.Printf("Error parsing TLE: %v\n", err)
		return
	}

	fmt.Println("=== Testing Combined SGP4/SDP4 Implementation ===")
	fmt.Printf("Satellite Number: %s\n", tle.SatelliteNumber)
	fmt.Printf("Mean Motion: %f revs/day\n", tle.MeanMotion)
	fmt.Printf("Orbital Period: %f minutes\n", tle.GetOrbitalPeriod())

	// Test the new combined implementation
	fmt.Println("\n--- Testing Combined Implementation ---")
	pos, vel, err := sgp4go.PropagateSatelliteCombined(tle, tle.JulianEpoch)
	if err != nil {
		fmt.Printf("Error with combined implementation: %v\n", err)
	} else {
		fmt.Printf("Position at epoch (Earth radii): X=%f, Y=%f, Z=%f\n", pos.X, pos.Y, pos.Z)
		fmt.Printf("Velocity at epoch (Earth radii/min): X=%f, Y=%f, Z=%f\n", vel.X, vel.Y, vel.Z)
	}

	// Test with kilometers
	posKm, velKm, err := sgp4go.PropagateSatelliteCombined(tle, tle.JulianEpoch, "km")
	if err != nil {
		fmt.Printf("Error with combined implementation (km): %v\n", err)
	} else {
		fmt.Printf("Position at epoch (km): X=%f, Y=%f, Z=%f\n", posKm.X, posKm.Y, posKm.Z)
		fmt.Printf("Velocity at epoch (km/s): X=%f, Y=%f, Z=%f\n", velKm.X, velKm.Y, velKm.Z)
	}

	// Compare with original implementation
	fmt.Println("\n--- Comparing with Original Implementation ---")
	posOrig, velOrig, err := sgp4go.PropagateSatellite(tle, tle.JulianEpoch)
	if err != nil {
		fmt.Printf("Error with original implementation: %v\n", err)
	} else {
		fmt.Printf("Original Position (Earth radii): X=%f, Y=%f, Z=%f\n", posOrig.X, posOrig.Y, posOrig.Z)
		fmt.Printf("Original Velocity (Earth radii/min): X=%f, Y=%f, Z=%f\n", velOrig.X, velOrig.Y, velOrig.Z)
	}

	// Test propagation to a future time
	fmt.Println("\n--- Testing Future Propagation ---")
	futureTime := tle.JulianEpoch + 1.0 // 1 day later
	posFuture, velFuture, err := sgp4go.PropagateSatelliteCombined(tle, futureTime, "km")
	if err != nil {
		fmt.Printf("Error with future propagation: %v\n", err)
	} else {
		fmt.Printf("Position 1 day later (km): X=%f, Y=%f, Z=%f\n", posFuture.X, posFuture.Y, posFuture.Z)
		fmt.Printf("Velocity 1 day later (km/s): X=%f, Y=%f, Z=%f\n", velFuture.X, velFuture.Y, velFuture.Z)
	}
}