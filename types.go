package sgp4go

import "time"

// Vector represents a 3D vector with magnitude
type Vector struct {
	X, Y, Z, Magnitude float64
}

// State represents position and velocity vectors
type State struct {
	Position Vector
	Velocity Vector
}

// TLE represents a Two-Line Element set
type TLE struct {
	SatelliteNumber   string
	Epoch             float64
	JulianEpoch       float64
	GregorianEpoch    time.Time
	MeanMotionDot     float64
	MeanMotionDDot    float64
	MeanMotionDDotExp int
	BStar             float64
	BStarExp          int
	ElementSet        string
	Inclination       float64
	RightAscension    float64
	Eccentricity      float64
	ArgumentPerigee   float64
	MeanAnomaly       float64
	MeanMotion        float64
}

// GetJulianEpoch returns the Julian epoch (equivalent to jdsatepoch in other packages)
func (t *TLE) GetJulianEpoch() float64 {
	return t.JulianEpoch
}

// GetEpoch returns the TLE epoch in the original format
func (t *TLE) GetEpoch() float64 {
	return t.Epoch
}

// GetSatelliteNumber returns the satellite number
func (t *TLE) GetSatelliteNumber() string {
	return t.SatelliteNumber
}

// GetElementSet returns the element set number
func (t *TLE) GetElementSet() string {
	return t.ElementSet
}

// GetOrbitalPeriod returns the orbital period in minutes
func (t *TLE) GetOrbitalPeriod() float64 {
	if t.MeanMotion > 0 {
		return 86400.0 / t.MeanMotion // Mean Motion in revolutions per day
	}
	return 0
}

// SatelliteData represents the internal satellite data used by SGP4
type SatelliteData struct {
	// Original TLE data
	TLE TLE

	// Converted orbital elements
	XMO, XNodeO, OmegaO, EO, XIncl, XNO, XNDT2O, XNDD6O, BStar float64
	JulianEpoch, XKE                                           float64

	// Model flags
	IFlag, IDeep int

	// Derived constants for SGP4
	EQSQ, SINIQ, COSIQ, RTEQSQ, AO, COSQ2, SINOMO, COSOMO float64
	BSQ, XLLDOT, OMGDT, XNODOT, XNODP                     float64

	// Deep space variables
	XLL, OMGASM, XNODES, EM, XINC, XN, T float64
	QOMS2T                               float64
}

// GetJulianEpoch returns the Julian epoch
func (s *SatelliteData) GetJulianEpoch() float64 {
	return s.JulianEpoch
}

// GetOrbitalPeriod returns the orbital period in minutes
func (s *SatelliteData) GetOrbitalPeriod() float64 {
	if s.XNO > 0 {
		return TWOPI / s.XNO
	}
	return 0
}

// IsDeepSpace returns true if this satellite requires the deep space model (SDP4)
func (s *SatelliteData) IsDeepSpace() bool {
	return s.IDeep == 1
}

// GetInclination returns the inclination in degrees
func (s *SatelliteData) GetInclination() float64 {
	return Degrees(s.XIncl)
}

// GetEccentricity returns the eccentricity
func (s *SatelliteData) GetEccentricity() float64 {
	return s.EO
}

// GetMeanMotion returns the mean motion in revolutions per day
func (s *SatelliteData) GetMeanMotion() float64 {
	return s.XNO * XMNPDA / TWOPI
}

// Satellite represents a satellite with its data and state
type Satellite struct {
	Name     string
	Data     SatelliteData
	Selected bool
}

// SatelliteDatabase represents a collection of satellites
type SatelliteDatabase struct {
	Satellites []Satellite
	Epoch      time.Time
}

// PropagationResult represents the result of orbital propagation
type PropagationResult struct {
	Time     time.Time
	Position Vector
	Velocity Vector
	Error    error
}
