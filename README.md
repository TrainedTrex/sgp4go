# SGP4Go

A Go implementation of the NORAD SGP4/SDP4 orbital propagation models, translated from the original code by David Vallado & Dr. TS Kelso.

> Beta version

## Features

- **Complete SGP4/SDP4 Implementation**: Full translation of the NORAD orbital models
- **TLE Parsing**: Parse standard Two-Line Element sets
- **Simple Propagation API**: Single function with optional units parameter
- **Flexible Units**: Optional units parameter - Earth radii (default) or kilometers
- **Testing**: Unit tests and validation against known results
- **Go Idioms**: Go error handling, type safety, and modern Go coding patterns

## Installation

```bash
go get github.com/TrainedTrex/sgp4go
```

## Quick Start

```go
package main

import (
    "fmt"
    "time"
    "github.com/TrainedTrex/sgp4go"
)

func main() {
    // Parse a TLE
    line1 := "1 25544U 98067A   25224.47423135  .00010254  00000+0  18430-3 0  9994"
    line2 := "2 25544  51.6349  28.0780 0001172 181.8871 178.2114 15.50477111523893"
    
    tle, err := sgp4.ParseTLE(line1, line2)
    if err != nil {
        panic(err)
    }

    // Access TLE fields
    fmt.Printf("Satellite: %s\n", tle.SatelliteNumber)
    fmt.Printf("Epoch: %f\n", tle.Epoch)
    fmt.Printf("Julian Epoch: %f\n", tle.JulianEpoch)
    fmt.Printf("Inclination: %f degrees\n", tle.Inclination)
    fmt.Printf("Orbital Period: %f minutes\n", tle.GetOrbitalPeriod())

    // Propagate using different methods
    now := time.Now()
    
    // Convert time to Julian date for propagation
    jd := sgp4.TimeToJulianDate(now)
    
    // Method 1: Earth radii (default)
    pos, vel, err := sgp4.PropagateSatellite(tle, jd)
    if err != nil {
        panic(err)
    }
    fmt.Printf("Position (Earth radii): X=%f, Y=%f, Z=%f\n", pos.X, pos.Y, pos.Z)
    
    // Method 2: Kilometers (with "km" parameter)
    posKm, velKm, err := sgp4.PropagateSatellite(tle, jd, "km")
    if err != nil {
        panic(err)
    }
    fmt.Printf("Position (km): X=%f, Y=%f, Z=%f\n", posKm.X, posKm.Y, posKm.Z)
    fmt.Printf("Velocity (km/s): X=%f, Y=%f, Z=%f\n", velKm.X, velKm.Y, velKm.Z)
}
```

## API Reference

### TLE Struct

```go
type TLE struct {
    SatelliteNumber   string  // Satellite catalog number
    Epoch             float64 // TLE epoch (year + day of year)
    JulianEpoch       float64 // Julian date of epoch
    MeanMotionDot     float64 // First derivative of mean motion
    MeanMotionDDot    float64 // Second derivative of mean motion
    MeanMotionDDotExp int     // Exponent for MeanMotionDDot
    BStar             float64 // B* drag term
    BStarExp          int     // Exponent for BStar
    ElementSet        string  // Element set number
    Inclination       float64 // Inclination (degrees)
    RightAscension    float64 // Right ascension of ascending node (degrees)
    Eccentricity      float64 // Eccentricity
    ArgumentPerigee   float64 // Argument of perigee (degrees)
    MeanAnomaly       float64 // Mean anomaly (degrees)
    MeanMotion        float64 // Mean motion (revolutions per day)
}
```

**Convenience Methods:**
- `GetJulianEpoch() float64` - Returns Julian epoch
- `GetEpoch() float64` - Returns TLE epoch
- `GetSatelliteNumber() string` - Returns satellite number
- `GetElementSet() string` - Returns element set number
- `GetOrbitalPeriod() float64` - Returns orbital period in minutes

### SatelliteData Struct

All fields are publicly accessible:

```go
type SatelliteData struct {
    TLE TLE // Original TLE data
    
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
```

**Convenience Methods:**
- `GetJulianEpoch() float64` - Returns Julian epoch
- `GetOrbitalPeriod() float64` - Returns orbital period in minutes
- `IsDeepSpace() bool` - Returns true if deep space model is needed
- `GetInclination() float64` - Returns inclination in degrees
- `GetEccentricity() float64` - Returns eccentricity
- `GetMeanMotion() float64` - Returns mean motion in revs/day

### Main Functions

#### TLE Parsing
```go
func ParseTLE(line1, line2 string) (*TLE, error)
func ParseFloat(s string) (float64, error)
```

#### Satellite Data Conversion
```go
func ConvertSatelliteData(tle *TLE) (*SatelliteData, error)
```

#### Propagation (High-Level)
```go
func PropagateSatellite(tle *TLE, time float64, units ...string) (Vector, Vector, error)
```

**Note:** The `units` parameter is optional. Pass `"km"` to get results in kilometers, or omit for Earth radii (default). For time.Time objects, use `TimeToJulianDate()` to convert first.

#### Propagation (Low-Level)
```go
func SGP(time float64, sat *SatelliteData) (Vector, Vector, error)
func SGP4(tsince float64, sat *SatelliteData) (Vector, Vector, error)
func SDP4(tsince float64, sat *SatelliteData) (Vector, Vector, error)
```

#### Time Conversions
```go
func JulianDateOfEpoch(epoch float64) float64
func JulianDate(year, month, day int) float64
func JulianDateToTime(jd float64) time.Time
func TimeToJulianDate(t time.Time) float64
func DaysSinceEpoch(julianEpoch, julianDate float64) float64
func MinutesSinceEpoch(julianEpoch, julianDate float64) float64
```

#### Mathematical Functions
```go
func Sign(arg float64) int
func Square(arg float64) float64
func Cube(arg float64) float64
func Power(arg, pwr float64) float64
func Radians(arg float64) float64
func Degrees(arg float64) float64
func Tan(arg float64) float64
func ArcSin(arg float64) float64
func ArcCos(arg float64) float64
func Modulus(arg1, arg2 float64) float64
func Fmod2p(arg float64) float64
func AcTan(sinx, cosx float64) float64
```

#### Vector Operations
```go
func Magnitude(v *Vector)
func VecAdd(v1, v2 Vector) Vector
func VecSub(v1, v2 Vector) Vector
func ScalarMultiply(k float64, v Vector) Vector
func Dot(v1, v2 Vector) float64
func Angle(v1, v2 Vector) float64
func Cross(v1, v2 Vector) Vector
func Normalize(v *Vector)
```

#### Unit Conversions
```go
func ConvertPositionToKilometers(pos *Vector)
func ConvertVelocityToKilometersPerSecond(vel *Vector)
func ConvertPositionAndVelocityToKilometers(pos, vel *Vector)
```

## Examples

See the `examples/` directory for complete examples:

- `main.go` - Basic usage example

## Testing

Run the test suite:

```bash
go test ./tests/...
```

## Implementation Notes

This is a direct translation from the original SGP4/SDP4 implementation by David Vallado and Dr. TS Kelso. The code maintains the same algorithms and numerical operations as the original while taking advantage of Go's type safety and error handling.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Original SGP4/SDP4 implementation by David Vallado & Dr. TS Kelso
- Based on NORAD SGP4 orbital models
- Translation maintains fidelity to the original algorithms
