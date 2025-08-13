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

	pos1, _, err := sgp4go.PropagateSatellite(tle, tle.JulianEpoch)
	if err != nil {
		fmt.Printf("Error propagating satellite: %v\n", err)
		return
	}
	fmt.Printf("Position at epoch: X=%f, Y=%f, Z=%f km\n", pos1.X, pos1.Y, pos1.Z)

	now := time.Now()
	pos2, _, err := sgp4go.PropagateSatelliteFromTime(tle, now)
	if err != nil {
		fmt.Printf("Error propagating satellite: %v\n", err)
		return
	}
	fmt.Printf("Position now: X=%f, Y=%f, Z=%f km\n", pos2.X, pos2.Y, pos2.Z)

	minutesSinceEpoch := 1440.0 // 1 day after epoch
	pos3, _, err := sgp4go.PropagateSatelliteFromMinutes(tle, minutesSinceEpoch)
	if err != nil {
		fmt.Printf("Error propagating satellite: %v\n", err)
		return
	}
	fmt.Printf("Position 1 day after epoch: X=%f, Y=%f, Z=%f km\n", pos3.X, pos3.Y, pos3.Z)

	fmt.Printf("\nAll Accessible Fields:\n")
	fmt.Printf("TLE.SatelliteNumber: %s\n", tle.SatelliteNumber)
	fmt.Printf("TLE.Epoch: %f\n", tle.Epoch)
	fmt.Printf("TLE.JulianEpoch: %f (equivalent to jdsatepoch)\n", tle.JulianEpoch)
	fmt.Printf("TLE.MeanMotionDot: %f\n", tle.MeanMotionDot)
	fmt.Printf("TLE.MeanMotionDDot: %f\n", tle.MeanMotionDDot)
	fmt.Printf("TLE.BStar: %f\n", tle.BStar)
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
	fmt.Printf("Sat.IFlag: %d\n", sat.IFlag)
	fmt.Printf("Sat.IDeep: %d\n", sat.IDeep)
}
