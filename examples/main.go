package main

import (
	"fmt"
	"time"

	"github.com/TrainedTrex/sgp4go"
)

func main() {
	// Example TLE for ISS
	line1 := "1 25544U 98067A   25224.47423135  .00010254  00000+0  18430-3 0  9994"
	line2 := "2 25544  51.6349  28.0780 0001172 181.8871 178.2114 15.50477111523893"

	// Parse the TLE
	tle, err := sgp4go.ParseTLE(line1, line2)
	if err != nil {
		fmt.Printf("Error parsing TLE: %v\n", err)
		return
	}

	fmt.Println("tle.Epoch", tle.Epoch)
	fmt.Println("tle.JulianEpoch", tle.JulianEpoch)
	fmt.Println("tle.Epoch", sgp4go.JulianDateToTime(tle.JulianEpoch))
	fmt.Println("tle.Epoch", sgp4go.JulianToTimeAstronomical(tle.JulianEpoch))

	// Demonstrate accessing TLE fields
	fmt.Printf("Satellite Number: %s\n", tle.SatelliteNumber)
	fmt.Printf("Epoch: %f\n", tle.Epoch)
	fmt.Printf("Julian Epoch: %f\n", tle.JulianEpoch)
	fmt.Printf("Inclination: %f degrees\n", tle.Inclination)
	fmt.Printf("Eccentricity: %f\n", tle.Eccentricity)
	fmt.Printf("Mean Motion: %f revs/day\n", tle.MeanMotion)
	fmt.Printf("Orbital Period: %f minutes\n", tle.GetOrbitalPeriod())

	// Convert to satellite data for propagation
	sat, err := sgp4go.ConvertSatelliteData(tle)
	if err != nil {
		fmt.Printf("Error converting satellite data: %v\n", err)
		return
	}

	// Satellite data fields
	fmt.Printf("\nSatellite Data:\n")
	fmt.Printf("Is Deep Space: %t\n", sat.IsDeepSpace())
	fmt.Printf("Inclination (degrees): %f\n", sat.GetInclination())
	fmt.Printf("Eccentricity: %f\n", sat.GetEccentricity())
	fmt.Printf("Mean Motion: %f revs/day\n", sat.GetMeanMotion())
	fmt.Printf("Orbital Period: %f minutes\n", sat.GetOrbitalPeriod())

	fmt.Printf("\nPropagation Examples:\n")

	pos1, vel1, err := sgp4go.PropagateSatellite(tle, tle.JulianEpoch)
	if err != nil {
		fmt.Printf("Error propagating satellite: %v\n", err)
		return
	}
	fmt.Printf("Position at epoch (Earth radii): X=%f, Y=%f, Z=%f\n", pos1.X, pos1.Y, pos1.Z)
	fmt.Printf("Velocity at epoch (Earth radii/min): X=%f, Y=%f, Z=%f\n", vel1.X, vel1.Y, vel1.Z)
	// Convert to kilometers
	sgp4go.ConvertPositionToKilometers(&pos1)
	sgp4go.ConvertVelocityToKilometersPerSecond(&vel1)
	fmt.Printf("Position at epoch (kilometers): X=%f, Y=%f, Z=%f\n", pos1.X, pos1.Y, pos1.Z)
	fmt.Printf("Velocity at epoch (km/s): X=%f, Y=%f, Z=%f\n", vel1.X, vel1.Y, vel1.Z)

	jnow := sgp4go.JulianToTimeAstronomical(tle.JulianEpoch)
	fmt.Println("jnow", tle.JulianEpoch)
	fmt.Println("jnow", jnow)
	jnowJul := sgp4go.TimeToJulianDate(jnow.Add(time.Minute * 92))
	fmt.Println("jnowJul", jnowJul)
	pos2, vel2, err := sgp4go.PropagateSatellite(tle, jnowJul)
	if err != nil {
		fmt.Printf("Error propagating satellite: %v\n", err)
		return
	}
	fmt.Printf("Position now (Earth radii): X=%f, Y=%f, Z=%f\n", pos2.X, pos2.Y, pos2.Z)
	fmt.Printf("Velocity now (Earth radii/min): X=%f, Y=%f, Z=%f\n", vel2.X, vel2.Y, vel2.Z)

	// Convert to kilometers using convenience functions - FIXED: use 92min time instead of now
	time92min := jnow.Add(time.Minute * 92)
	pos2Km, vel2Km, err := sgp4go.PropagateSatelliteFromTimeInKilometers(tle, time92min)
	if err != nil {
		fmt.Printf("Error propagating satellite: %v\n", err)
		return
	}
	fmt.Printf("Position now (kilometers): X=%f, Y=%f, Z=%f\n", pos2Km.X, pos2Km.Y, pos2Km.Z)
	fmt.Printf("Velocity now (km/s): X=%f, Y=%f, Z=%f\n", vel2Km.X, vel2Km.Y, vel2Km.Z)

	fmt.Printf("\nAll Accessible Fields:\n")
	fmt.Printf("TLE.SatelliteNumber: %s\n", tle.SatelliteNumber)
	fmt.Printf("TLE.Epoch: %f\n", tle.Epoch)
	fmt.Printf("TLE.JulianEpoch: %f\n", tle.JulianEpoch)
	fmt.Printf("TLE.MeanMotionDot: %f\n", tle.MeanMotionDot)
	fmt.Printf("TLE.MeanMotionDDot: %f\n", tle.MeanMotionDDot)
	fmt.Printf("TLE.BStar: %f\n", tle.BStar)
	fmt.Printf("TLE.BStarExp: %d\n", tle.BStarExp)
	fmt.Printf("TLE.ElementSet: %s\n", tle.ElementSet)
	fmt.Printf("TLE.Inclination: %f\n", tle.Inclination)
	fmt.Printf("TLE.RightAscension: %f\n", tle.RightAscension)
	fmt.Printf("TLE.Eccentricity: %f\n", tle.Eccentricity)
	fmt.Printf("TLE.ArgumentPerigee: %f\n", tle.ArgumentPerigee)
	fmt.Printf("TLE.MeanAnomaly: %f\n", tle.MeanAnomaly)
	fmt.Printf("TLE.MeanMotion: %f\n", tle.MeanMotion)

	fmt.Printf("\nSatelliteData fields:\n")
	fmt.Printf("Sat.XMO: %f\n", sat.XMO)
	fmt.Printf("Sat.XNodeO: %f\n", sat.XNodeO)
	fmt.Printf("Sat.OmegaO: %f\n", sat.OmegaO)
	fmt.Printf("Sat.EO: %f\n", sat.EO)
	fmt.Printf("Sat.XIncl: %f\n", sat.XIncl)
	fmt.Printf("Sat.XNO: %f\n", sat.XNO)
	fmt.Printf("Sat.JulianEpoch: %f\n", sat.JulianEpoch)
	fmt.Printf("Sat.XKE: %f\n", sat.XKE)
	fmt.Printf("Sat.BStar: %f\n", sat.BStar)
	fmt.Printf("Sat.Eccentricity: %f\n", sat.EO)
	fmt.Printf("Sat.IFlag: %d\n", sat.IFlag)
	fmt.Printf("Sat.IDeep: %d\n", sat.IDeep)
}
