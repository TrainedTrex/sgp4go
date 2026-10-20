package main

import (
	"fmt"
	"time"

	"github.com/TrainedTrex/sgp4go"
)

// OMM JSON or GP JSON (plain numbers) parses the same way.
const ommJSON = `[{
	"OBJECT_NAME": "STARLINK-38128",
	"OBJECT_ID": "2026-160A",
	"EPOCH": "2026-10-19T11:54:15.765120",
	"MEAN_MOTION": "15.30428694",
	"ECCENTRICITY": "0.00024627",
	"INCLINATION": "69.9982",
	"RA_OF_ASC_NODE": "117.4636",
	"ARG_OF_PERICENTER": "271.8593",
	"MEAN_ANOMALY": "88.2301",
	"NORAD_CAT_ID": "100001",
	"ELEMENT_SET_NO": "999",
	"REV_AT_EPOCH": "1644",
	"BSTAR": "0.00059883618000",
	"MEAN_MOTION_DOT": "0.00017582",
	"MEAN_MOTION_DDOT": "0.0000000000000"
}]`

func main() {
	oms, err := sgp4go.ParseOMM([]byte(ommJSON))
	if err != nil {
		fmt.Printf("Error parsing OMM: %v\n", err)
		return
	}
	om := oms[0]
	fmt.Printf("Object: %s (%s), NORAD %d\n", om.ObjectName, om.ObjectID, om.NoradCatID)
	fmt.Printf("Epoch: %s\n", om.Epoch.Format(time.RFC3339Nano))

	// Convert to a TLE and propagate as usual
	tle := om.TLE()
	jd := sgp4go.TimeToJulianDate(om.Epoch.Add(90 * time.Minute))
	pos, vel, err := sgp4go.PropagateSatellite(tle, jd, "km")
	if err != nil {
		fmt.Printf("Error propagating satellite: %v\n", err)
		return
	}
	fmt.Printf("Position at epoch+90min (km): X=%f, Y=%f, Z=%f\n", pos.X, pos.Y, pos.Z)
	fmt.Printf("Velocity at epoch+90min (km/s): X=%f, Y=%f, Z=%f\n", vel.X, vel.Y, vel.Z)

	// The same object as an Alpha-5 TLE: "A0001" decodes to catalog number 100001
	alpha, err := sgp4go.ParseTLE(
		"1 A0001U 26160A   26292.49601580  .00017582  00000-0  59884-3 0  9990",
		"2 A0001  69.9982 117.4636 0002463 271.8593  88.2301 15.30428694 16443")
	if err != nil {
		fmt.Printf("Error parsing TLE: %v\n", err)
		return
	}
	posTLE, _, err := sgp4go.PropagateSatellite(alpha, jd, "km")
	if err != nil {
		fmt.Printf("Error propagating satellite: %v\n", err)
		return
	}
	fmt.Printf("\nAlpha-5 TLE %q -> catalog number %d\n", alpha.SatelliteNumber, alpha.CatalogNumber)
	fmt.Printf("TLE position at epoch+90min (km): X=%f, Y=%f, Z=%f\n", posTLE.X, posTLE.Y, posTLE.Z)
}
