package sgp4go

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

// ParseTLE parses a two-line element set and returns a TLE struct
func ParseTLE(line1, line2 string) (*TLE, error) {
	// Pad lines to 69 characters if needed
	if len(line1) < 69 {
		line1 = line1 + strings.Repeat(" ", 69-len(line1))
	}
	if len(line2) < 69 {
		line2 = line2 + strings.Repeat(" ", 69-len(line2))
	}

	tle := &TLE{}

	// Parse line 1
	var err error
	tle.SatelliteNumber = strings.TrimSpace(line1[2:7])

	tle.Epoch, err = ParseFloat(line1[18:32])
	if err != nil {
		return nil, fmt.Errorf("error parsing epoch: %v", err)
	}

	tle.MeanMotionDot, err = ParseFloat(line1[33:43])
	if err != nil {
		return nil, fmt.Errorf("error parsing mean motion dot: %v", err)
	}

	// Parse mean motion ddot (with exponent)
	// Format: 6 digits followed by exponent in next 2 positions
	ddotStr := line1[44:50]
	ddotExpStr := strings.TrimSpace(line1[50:52])
	if len(strings.TrimSpace(ddotStr)) >= 5 {
		ddotVal, err := ParseFloat(ddotStr)
		if err != nil {
			return nil, fmt.Errorf("error parsing mean motion ddot: %v", err)
		}
		exp, err := strconv.Atoi(ddotExpStr)
		if err != nil {
			return nil, fmt.Errorf("error parsing mean motion ddot exponent: %v", err)
		}
		tle.MeanMotionDDot = ddotVal * 1e-5
		tle.MeanMotionDDotExp = exp
	}

	// Parse B*
	// Format: 6 digits followed by exponent in next 2 positions
	bstarStr := line1[53:59]
	bstarExpStr := strings.TrimSpace(line1[59:61])
	if len(strings.TrimSpace(bstarStr)) >= 5 {
		bstarVal, err := ParseFloat(bstarStr)
		if err != nil {
			return nil, fmt.Errorf("error parsing B*: %v", err)
		}
		exp, err := strconv.Atoi(bstarExpStr)
		if err != nil {
			return nil, fmt.Errorf("error parsing B* exponent: %v", err)
		}
		tle.BStar = bstarVal * 1e-5
		tle.BStarExp = exp
	}

	elsetNum, err := strconv.Atoi(strings.TrimSpace(line1[64:68]))
	if err != nil {
		return nil, fmt.Errorf("error parsing element set number: %v", err)
	}
	tle.ElementSet = fmt.Sprintf("%03d", elsetNum)

	// Parse line 2
	tle.Inclination, err = ParseFloat(line2[8:16])
	if err != nil {
		return nil, fmt.Errorf("error parsing inclination: %v", err)
	}

	tle.RightAscension, err = ParseFloat(line2[17:25])
	if err != nil {
		return nil, fmt.Errorf("error parsing right ascension: %v", err)
	}

	// Parse eccentricity (with implied decimal point)
	eccStr := strings.TrimSpace(line2[26:33])
	if len(eccStr) >= 7 {
		eccVal, err := ParseFloat("0." + eccStr)
		if err != nil {
			return nil, fmt.Errorf("error parsing eccentricity: %v", err)
		}
		tle.Eccentricity = eccVal
	}

	tle.ArgumentPerigee, err = ParseFloat(line2[34:42])
	if err != nil {
		return nil, fmt.Errorf("error parsing argument of perigee: %v", err)
	}

	tle.MeanAnomaly, err = ParseFloat(line2[43:51])
	if err != nil {
		return nil, fmt.Errorf("error parsing mean anomaly: %v", err)
	}

	tle.MeanMotion, err = ParseFloat(line2[52:63])
	if err != nil {
		return nil, fmt.Errorf("error parsing mean motion: %v", err)
	}

	// Calculate Julian epoch
	tle.JulianEpoch = JulianDateOfEpoch(tle.Epoch)

	return tle, nil
}

// ParseFloat parses a float from a string, handling various formats
func ParseFloat(s string) (float64, error) {
	s = strings.TrimSpace(s)
	if s == "" {
		return 0.0, nil
	}

	// Handle scientific notation
	if strings.Contains(s, "E") || strings.Contains(s, "e") {
		return strconv.ParseFloat(s, 64)
	}

	// Handle regular decimal numbers
	return strconv.ParseFloat(s, 64)
}

// ConvertSatelliteData converts TLE data to the internal format used by SGP4
func ConvertSatelliteData(tle *TLE) (*SatelliteData, error) {
	if tle == nil {
		return nil, fmt.Errorf("TLE cannot be nil")
	}

	sat := &SatelliteData{
		TLE: *tle,
	}

	// Convert orbital elements to proper units
	sat.XNDD6O = tle.MeanMotionDDot * Power(10.0, float64(tle.MeanMotionDDotExp))
	sat.BStar = tle.BStar * Power(10.0, float64(tle.BStarExp)) / AE
	sat.XNodeO = Radians(tle.RightAscension)
	sat.OmegaO = Radians(tle.ArgumentPerigee)
	sat.XMO = Radians(tle.MeanAnomaly)
	sat.XIncl = Radians(tle.Inclination)
	sat.XNO = tle.MeanMotion * TWOPI / XMNPDA
	sat.XNDT2O = tle.MeanMotionDot * TWOPI / Square(XMNPDA)
	sat.XNDD6O = sat.XNDD6O * TWOPI / Cube(XMNPDA)
	sat.EO = tle.Eccentricity
	sat.JulianEpoch = tle.JulianEpoch

	// Calculate derived constants
	sat.XKE = CalculateXKE()

	// Determine whether Deep-Space Model is needed
	a1 := Power(sat.XKE/sat.XNO, TOTHRD)
	temp := 1.5 * CK2 * (3*Power(math.Cos(sat.XIncl), 2) - 1) / Power(1-sat.EO*sat.EO, 1.5)
	del1 := temp / (a1 * a1)
	ao := a1 * (1 - del1*(0.5*TOTHRD+del1*(1+134.0/81.0*del1)))
	delo := temp / (ao * ao)
	xnodp := sat.XNO / (1 + delo)

	if (TWOPI / xnodp) >= 225.0 {
		sat.IDeep = 1
	} else {
		sat.IDeep = 0
	}

	sat.IFlag = 0

	return sat, nil
}

// CalculateXKE calculates the Earth's gravitational parameter in Earth radii^3/min^2
func CalculateXKE() float64 {
	return math.Sqrt(3600.0 * GE / Cube(XKMPER))
}

// ConvertSatState converts position and velocity from Earth radii to kilometers
func ConvertSatState(pos, vel *Vector) {
	// Convert position from Earth radii to kilometers
	pos.X *= XKMPER
	pos.Y *= XKMPER
	pos.Z *= XKMPER

	// Convert velocity from Earth radii/minute to kilometers/second
	vel.X *= XKMPER / 60.0
	vel.Y *= XKMPER / 60.0
	vel.Z *= XKMPER / 60.0

	// Recalculate magnitudes
	Magnitude(pos)
	Magnitude(vel)
}
