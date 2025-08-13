package main

import (
	"math"
	"testing"
	"time"

	"github.com/TrainedTrex/sgp4go"
)

func TestParseTLE(t *testing.T) {
	// Test TLE for a known satellite
	line1 := "1 25544U 98067A   25224.47423135  .00010254  00000+0  18430-3 0  9994"
	line2 := "2 25544  51.6349  28.0780 0001172 181.8871 178.2114 15.50477111523893"

	tle, err := sgp4go.ParseTLE(line1, line2)
	if err != nil {
		t.Fatalf("ParseTLE failed: %v", err)
	}

	// Check parsed values
	if tle.SatelliteNumber != "25544" {
		t.Errorf("Expected satellite number 25544, got %s", tle.SatelliteNumber)
	}

	// The epoch should be 25224.47423135 (from the TLE)
	if math.Abs(tle.Epoch-25224.47423135) > 0.001 {
		t.Errorf("Expected epoch 25224.47423135, got %f", tle.Epoch)
	}

	// The inclination should be 51.6349 (from the TLE)
	if math.Abs(tle.Inclination-51.6349) > 0.001 {
		t.Errorf("Expected inclination 51.6349, got %f", tle.Inclination)
	}

	// The eccentricity should be 0.0001172 (from the TLE)
	if math.Abs(tle.Eccentricity-0.0001172) > 0.0000001 {
		t.Errorf("Expected eccentricity 0.0001172, got %f", tle.Eccentricity)
	}

	// The mean motion should be 15.504771 (from the TLE)
	if math.Abs(tle.MeanMotion-15.504771) > 0.001 {
		t.Errorf("Expected mean motion 15.504771, got %f", tle.MeanMotion)
	}
}

func TestConvertSatelliteData(t *testing.T) {
	// Create a test TLE
	tle := &sgp4go.TLE{
		SatelliteNumber: "25544",
		Epoch:           25224.47423135,
		Inclination:     51.6349,
		Eccentricity:    0.0001172,
		MeanMotion:      15.504771,
		RightAscension:  28.0780,
		ArgumentPerigee: 181.8871,
		MeanAnomaly:     178.2114,
		MeanMotionDot:   0.0,
		MeanMotionDDot:  0.00010254,
		BStar:           0.18430e-3,
		ElementSet:      "999",
	}

	satData, err := sgp4go.ConvertSatelliteData(tle)
	if err != nil {
		t.Fatalf("ConvertSatelliteData failed: %v", err)
	}

	// Check that conversion was successful
	// IFlag should be 0 for SGP4 (near-earth) satellites
	if satData.IFlag != 0 {
		t.Errorf("Expected IFlag 0, got %d", satData.IFlag)
	}

	if satData.XKE <= 0 {
		t.Errorf("Expected positive XKE, got %f", satData.XKE)
	}

	// Check that angles were converted to radians
	if math.Abs(satData.XIncl-sgp4go.Radians(51.6349)) > 0.001 {
		t.Errorf("Expected inclination in radians, got %f", satData.XIncl)
	}
}

func TestMathematicalFunctions(t *testing.T) {
	// Test Power function
	result := sgp4go.Power(2.0, 3.0)
	expected := 8.0
	if math.Abs(result-expected) > 0.001 {
		t.Errorf("Power(2,3) expected %f, got %f", expected, result)
	}

	// Test Cube function
	result = sgp4go.Cube(3.0)
	expected = 27.0
	if math.Abs(result-expected) > 0.001 {
		t.Errorf("Cube(3) expected %f, got %f", expected, result)
	}

	// Test Fmod2p function
	result = sgp4go.Fmod2p(3 * sgp4go.PI)
	expected = sgp4go.PI
	if math.Abs(result-expected) > 0.001 {
		t.Errorf("Fmod2p(3π) expected %f, got %f", expected, result)
	}

	// Test AcTan function
	result = sgp4go.AcTan(1.0, 1.0)
	expected = sgp4go.PI / 4.0
	if math.Abs(result-expected) > 0.001 {
		t.Errorf("AcTan(1,1) expected %f, got %f", expected, result)
	}
}

func TestVectorOperations(t *testing.T) {
	v1 := sgp4go.Vector{X: 1.0, Y: 2.0, Z: 3.0}
	v2 := sgp4go.Vector{X: 4.0, Y: 5.0, Z: 6.0}

	// Test vector addition
	result := sgp4go.VecAdd(v1, v2)
	expected := sgp4go.Vector{X: 5.0, Y: 7.0, Z: 9.0}
	if result.X != expected.X || result.Y != expected.Y || result.Z != expected.Z {
		t.Errorf("VecAdd expected %v, got %v", expected, result)
	}

	// Test vector subtraction
	result = sgp4go.VecSub(v2, v1)
	expected = sgp4go.Vector{X: 3.0, Y: 3.0, Z: 3.0}
	if result.X != expected.X || result.Y != expected.Y || result.Z != expected.Z {
		t.Errorf("VecSub expected %v, got %v", expected, result)
	}

	// Test dot product
	resultDot := sgp4go.Dot(v1, v2)
	expectedDot := 1.0*4.0 + 2.0*5.0 + 3.0*6.0
	if math.Abs(resultDot-expectedDot) > 0.001 {
		t.Errorf("Dot expected %f, got %f", expectedDot, resultDot)
	}

	// Test magnitude calculation
	sgp4go.Magnitude(&v1)
	expectedMag := math.Sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0)
	if math.Abs(v1.Magnitude-expectedMag) > 0.001 {
		t.Errorf("Magnitude expected %f, got %f", expectedMag, v1.Magnitude)
	}
}

func TestTimeFunctions(t *testing.T) {
	// Test Julian date conversion
	epoch := 25224.47423135 // 2025-08-12
	jd := sgp4go.JulianDateOfEpoch(epoch)

	// The result should be a reasonable Julian date (around 2460100 for 2025)
	if jd < 2400000 || jd > 2500000 {
		t.Errorf("Julian date %f seems unreasonable for epoch %f", jd, epoch)
	}

	// Test time conversion with a fixed time to avoid precision issues
	fixedTime := time.Date(2025, 8, 12, 0, 0, 0, 0, time.UTC) // Use midnight to avoid time-of-day issues
	jd2 := sgp4go.TimeToJulianDate(fixedTime)

	// Convert back and check
	time2 := sgp4go.JulianDateToTime(jd2)

	// Allow for some precision loss in the conversion
	diff := math.Abs(float64(fixedTime.Unix()) - float64(time2.Unix()))
	if diff > 1 { // Allow 1 second difference
		t.Errorf("Time conversion error too large: %f seconds", diff)
	}

	// Test that the date components match
	if fixedTime.Year() != time2.Year() ||
		fixedTime.Month() != time2.Month() ||
		fixedTime.Day() != time2.Day() {
		t.Errorf("Date conversion failed: expected %v, got %v", fixedTime.Format("2000-01-02"), time2.Format("2000-01-02"))
	}
}
