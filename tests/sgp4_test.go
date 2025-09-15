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

func TestSGP4Validation(t *testing.T) {
	// Test TLE for ISS
	line1 := "1 25544U 98067A   25224.47423135  .00010254  00000+0  18430-3 0  9994"
	line2 := "2 25544  51.6349  28.0780 0001172 181.8871 178.2114 15.50477111523893"

	tle, err := sgp4go.ParseTLE(line1, line2)
	if err != nil {
		t.Fatalf("ParseTLE failed: %v", err)
	}

	// Expected values from the verified true solution
	expectedEpoch := sgp4go.Vector{
		X: 5995.284156,
		Y: 3198.229573,
		Z: 0.004114,
	}
	expectedEpochVel := sgp4go.Vector{
		X: -2.243097,
		Y: 4.190838,
		Z: 6.009074,
	}

	expected92min := sgp4go.Vector{
		X: 6112.175919,
		Y: 2954.952808,
		Z: -293.036839,
	}
	expected92minVel := sgp4go.Vector{
		X: -1.843258,
		Y: 4.393272,
		Z: 5.999953,
	}

	// Test 1: Position at epoch
	posEpoch, velEpoch, err := sgp4go.PropagateSatellite(tle, tle.JulianEpoch, "km")
	t.Logf("Epoch: %v", sgp4go.JulianToTimeAstronomical(tle.JulianEpoch))

	if err != nil {
		t.Fatalf("PropagateSatellite at epoch failed: %v", err)
	}

	// Always show the values for comparison
	t.Logf("=== EPOCH COMPARISON ===")
	t.Logf("Expected Position: X=%f, Y=%f, Z=%f", expectedEpoch.X, expectedEpoch.Y, expectedEpoch.Z)
	t.Logf("Got Position:      X=%f, Y=%f, Z=%f", posEpoch.X, posEpoch.Y, posEpoch.Z)
	t.Logf("Position Diff:     X=%f, Y=%f, Z=%f", posEpoch.X-expectedEpoch.X, posEpoch.Y-expectedEpoch.Y, posEpoch.Z-expectedEpoch.Z)

	t.Logf("Expected Velocity: X=%f, Y=%f, Z=%f", expectedEpochVel.X, expectedEpochVel.Y, expectedEpochVel.Z)
	t.Logf("Got Velocity:      X=%f, Y=%f, Z=%f", velEpoch.X, velEpoch.Y, velEpoch.Z)
	t.Logf("Velocity Diff:     X=%f, Y=%f, Z=%f", velEpoch.X-expectedEpochVel.X, velEpoch.Y-expectedEpochVel.Y, velEpoch.Z-expectedEpochVel.Z)

	// Check epoch position with very strict tolerance
	posError := math.Sqrt(math.Pow(posEpoch.X-expectedEpoch.X, 2) +
		math.Pow(posEpoch.Y-expectedEpoch.Y, 2) +
		math.Pow(posEpoch.Z-expectedEpoch.Z, 2))

	if posError > 0.01 { // Allow 10m error (more realistic for SGP4)
		t.Errorf("Epoch position error too large: %f km", posError)
	}

	// Check epoch velocity with very strict tolerance
	velError := math.Sqrt(math.Pow(velEpoch.X-expectedEpochVel.X, 2) +
		math.Pow(velEpoch.Y-expectedEpochVel.Y, 2) +
		math.Pow(velEpoch.Z-expectedEpochVel.Z, 2))

	if velError > 0.001 { // Allow 1 m/s error (more realistic for SGP4)
		t.Errorf("Epoch velocity error too large: %f km/s", velError)
	}

	// Test 2: Position after 92 minutes
	// Calculate the exact time for 92 minutes after epoch
	epochTime := sgp4go.JulianToTimeAstronomical(tle.JulianEpoch)
	time92min := epochTime.Add(time.Minute * 92)
	jd92min := sgp4go.TimeToJulianDate(time92min)

	pos92min, vel92min, err := sgp4go.PropagateSatellite(tle, jd92min, "km")
	if err != nil {
		t.Fatalf("PropagateSatellite at 92min failed: %v", err)
	}
	t.Logf("92min: %v", sgp4go.JulianToTimeAstronomical(jd92min))

	// Always show the values for comparison
	t.Logf("=== 92 MINUTES COMPARISON ===")
	t.Logf("Expected Position: X=%f, Y=%f, Z=%f", expected92min.X, expected92min.Y, expected92min.Z)
	t.Logf("Got Position:      X=%f, Y=%f, Z=%f", pos92min.X, pos92min.Y, pos92min.Z)
	t.Logf("Position Diff:     X=%f, Y=%f, Z=%f", pos92min.X-expected92min.X, pos92min.Y-expected92min.Y, pos92min.Z-expected92min.Z)

	t.Logf("Expected Velocity: X=%f, Y=%f, Z=%f", expected92minVel.X, expected92minVel.Y, expected92minVel.Z)
	t.Logf("Got Velocity:      X=%f, Y=%f, Z=%f", vel92min.X, vel92min.Y, vel92min.Z)
	t.Logf("Velocity Diff:     X=%f, Y=%f, Z=%f", vel92min.X-expected92minVel.X, vel92min.Y-expected92minVel.Y, vel92min.Z-expected92minVel.Z)

	// Check 92min position with very strict tolerance
	posError92 := math.Sqrt(math.Pow(pos92min.X-expected92min.X, 2) +
		math.Pow(pos92min.Y-expected92min.Y, 2) +
		math.Pow(pos92min.Z-expected92min.Z, 2))

	if posError92 > 0.01 { // Allow 10m error (more realistic for SGP4)
		t.Errorf("92min position error too large: %f km", posError92)
	}

	// Check 92min velocity with very strict tolerance
	velError92 := math.Sqrt(math.Pow(vel92min.X-expected92minVel.X, 2) +
		math.Pow(vel92min.Y-expected92minVel.Y, 2) +
		math.Pow(vel92min.Z-expected92minVel.Z, 2))

	if velError92 > 0.001 { // Allow 1 m/s error (more realistic for SGP4)
		t.Errorf("92min velocity error too large: %f km/s", velError92)
	}

	// Debug output
	t.Logf("Epoch time: %v", epochTime)
	t.Logf("92min time: %v", time92min)
	t.Logf("Epoch JD: %f", tle.JulianEpoch)
	t.Logf("92min JD: %f", jd92min)
	t.Logf("Time difference: %f minutes", (jd92min-tle.JulianEpoch)*1440.0)
}

func TestSGP4ValidationGEO(t *testing.T) {
	// Test TLE for GEO
	line1 := "1 28884U 05041A   25258.23425557 -.00000019  00000-0  00000-0 0  9996"
	line2 := "2 28884   2.8166  80.2739 0017303 131.7673 259.8922  0.98940294 72726"

	tle, err := sgp4go.ParseTLE(line1, line2)
	if err != nil {
		t.Fatalf("ParseTLE failed: %v", err)
	}

	// Expected values from the verified true solution
	expectedEpoch := sgp4go.Vector{
		X: -15729.706099,
		Y: 39525.614082,
		Z: 1079.884828,
	}
	expectedEpochVel := sgp4go.Vector{
		X: -2.837587,
		Y: -1.138411,
		Z: 0.128763,
	}

	expected92min := sgp4go.Vector{
		X: -29760.145461,
		Y: 30329.305621,
		Z: 1688.168220,
	}
	expected92minVel := sgp4go.Vector{
		X: -2.178714,
		Y: -2.150033,
		Z: 0.088681,
	}

	// Test 1: Position at epoch
	posEpoch, velEpoch, err := sgp4go.PropagateSatellite(tle, tle.JulianEpoch, "km")
	if err != nil {
		t.Fatalf("PropagateSatellite at epoch failed: %v", err)
	}

	// Always show the values for comparison
	t.Logf("=== EPOCH COMPARISON ===")
	t.Logf("Expected Position: X=%f, Y=%f, Z=%f", expectedEpoch.X, expectedEpoch.Y, expectedEpoch.Z)
	t.Logf("Got Position:      X=%f, Y=%f, Z=%f", posEpoch.X, posEpoch.Y, posEpoch.Z)
	t.Logf("Position Diff:     X=%f, Y=%f, Z=%f", posEpoch.X-expectedEpoch.X, posEpoch.Y-expectedEpoch.Y, posEpoch.Z-expectedEpoch.Z)

	t.Logf("Expected Velocity: X=%f, Y=%f, Z=%f", expectedEpochVel.X, expectedEpochVel.Y, expectedEpochVel.Z)
	t.Logf("Got Velocity:      X=%f, Y=%f, Z=%f", velEpoch.X, velEpoch.Y, velEpoch.Z)
	t.Logf("Velocity Diff:     X=%f, Y=%f, Z=%f", velEpoch.X-expectedEpochVel.X, velEpoch.Y-expectedEpochVel.Y, velEpoch.Z-expectedEpochVel.Z)

	// Check epoch position with very strict tolerance
	posError := math.Sqrt(math.Pow(posEpoch.X-expectedEpoch.X, 2) +
		math.Pow(posEpoch.Y-expectedEpoch.Y, 2) +
		math.Pow(posEpoch.Z-expectedEpoch.Z, 2))

	if posError > 0.01 {
		t.Errorf("Epoch position error too large: %f km", posError)
	}

	// Check epoch velocity with very strict tolerance
	velError := math.Sqrt(math.Pow(velEpoch.X-expectedEpochVel.X, 2) +
		math.Pow(velEpoch.Y-expectedEpochVel.Y, 2) +
		math.Pow(velEpoch.Z-expectedEpochVel.Z, 2))

	if velError > 0.001 {
		t.Errorf("Epoch velocity error too large: %f km/s", velError)
	}

	// Test 2: Position after 92 minutes
	// Calculate the exact time for 92 minutes after epoch
	epochTime := sgp4go.JulianToTimeAstronomical(tle.JulianEpoch)
	time92min := epochTime.Add(time.Minute * 92)
	jd92min := sgp4go.TimeToJulianDate(time92min)

	pos92min, vel92min, err := sgp4go.PropagateSatellite(tle, jd92min, "km")
	if err != nil {
		t.Fatalf("PropagateSatellite at 92min failed: %v", err)
	}

	// Always show the values for comparison
	t.Logf("=== 92 MINUTES COMPARISON ===")
	t.Logf("Expected Position: X=%f, Y=%f, Z=%f", expected92min.X, expected92min.Y, expected92min.Z)
	t.Logf("Got Position:      X=%f, Y=%f, Z=%f", pos92min.X, pos92min.Y, pos92min.Z)
	t.Logf("Position Diff:     X=%f, Y=%f, Z=%f", pos92min.X-expected92min.X, pos92min.Y-expected92min.Y, pos92min.Z-expected92min.Z)

	t.Logf("Expected Velocity: X=%f, Y=%f, Z=%f", expected92minVel.X, expected92minVel.Y, expected92minVel.Z)
	t.Logf("Got Velocity:      X=%f, Y=%f, Z=%f", vel92min.X, vel92min.Y, vel92min.Z)
	t.Logf("Velocity Diff:     X=%f, Y=%f, Z=%f", vel92min.X-expected92minVel.X, vel92min.Y-expected92minVel.Y, vel92min.Z-expected92minVel.Z)

	// Check 92min position with very strict tolerance
	posError92 := math.Sqrt(math.Pow(pos92min.X-expected92min.X, 2) +
		math.Pow(pos92min.Y-expected92min.Y, 2) +
		math.Pow(pos92min.Z-expected92min.Z, 2))

	if posError92 > 0.01 {
		t.Errorf("92min position error too large: %f km", posError92)
	}

	// Check 92min velocity with very strict tolerance
	velError92 := math.Sqrt(math.Pow(vel92min.X-expected92minVel.X, 2) +
		math.Pow(vel92min.Y-expected92minVel.Y, 2) +
		math.Pow(vel92min.Z-expected92minVel.Z, 2))

	if velError92 > 0.001 {
		t.Errorf("92min velocity error too large: %f km/s", velError92)
	}

	// Debug output
	t.Logf("Epoch time: %v", epochTime)
	t.Logf("92min time: %v", time92min)
	t.Logf("Epoch JD: %f", tle.JulianEpoch)
	t.Logf("92min JD: %f", jd92min)
	t.Logf("Time difference: %f minutes", (jd92min-tle.JulianEpoch)*1440.0)
}
