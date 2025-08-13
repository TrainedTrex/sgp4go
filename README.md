# SGP4Go

A Go implementation of the NORAD SGP4/SDP4 orbital propagation models, translated from the original Pascal code by Dr. TS Kelso.

## Features

- **Complete SGP4/SDP4 Implementation**: Full translation of the NORAD orbital models
- **TLE Parsing**: Parse standard Two-Line Element sets
- **Multiple Propagation Methods**: Support for Julian dates, Go time.Time, and minutes since epoch
- **User-Friendly API**: All fields publicly accessible with convenience methods
- **Comprehensive Testing**: Unit tests and validation against known results
- **Go Idioms**: Proper error handling, type safety, and modern Go patterns

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

    // Access TLE fields (all publicly accessible!)
    fmt.Printf("Satellite: %s\n", tle.SatelliteNumber)
    fmt.Printf("Epoch: %f\n", tle.Epoch)
    fmt.Printf("Julian Epoch: %f\n", tle.JulianEpoch) // Equivalent to jdsatepoch
    fmt.Printf("Inclination: %f degrees\n", tle.Inclination)
    fmt.Printf("Orbital Period: %f minutes\n", tle.GetOrbitalPeriod())

    // Propagate using different methods
    now := time.Now()
    
    // Method 1: Using Go time.Time
    pos, vel, err := sgp4.PropagateSatelliteFromTime(tle, now)
    if err != nil {
        panic(err)
    }
    fmt.Printf("Position: X=%f, Y=%f, Z=%f km\n", pos.X, pos.Y, pos.Z)
}
```

## API Reference

### TLE Struct

All fields are publicly accessible:

```go
type TLE struct {
    SatelliteNumber   string  // Satellite catalog number
    Epoch             float64 // TLE epoch (year + day of year)
    JulianEpoch       float64 // Julian date of epoch (equivalent to jdsatepoch)
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
func CalculateXKE() float64
```

#### Propagation (High-Level)
```go
func PropagateSatellite(tle *TLE, time float64) (Vector, Vector, error)
func PropagateSatelliteFromTime(tle *TLE, t time.Time) (Vector, Vector, error)
func PropagateSatelliteFromMinutes(tle *TLE, minutesSinceEpoch float64) (Vector, Vector, error)
```

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

## Examples

See the `examples/` directory for complete examples:

- `main.go` - Basic usage example
- `user_friendly_example.go` - Demonstrates all accessible fields and methods

## Testing

Run the test suite:

```bash
go test ./tests/...
```

## Implementation Notes

This is a direct translation from the original Pascal SGP4/SDP4 implementation by Dr. TS Kelso. The code maintains the same algorithms and numerical precision as the original while taking advantage of Go's type safety and error handling.

### Key Differences from Original

1. **No Global Variables**: All state is encapsulated in structs
2. **Public Fields**: All important fields are publicly accessible (capitalized)
3. **Error Handling**: Proper Go-style error handling instead of returning 0.0
4. **Type Safety**: Compile-time type checking
5. **Convenience Methods**: Added getter methods for commonly needed values

### User-Friendly Features

Unlike some other SGP4 implementations, this package ensures that:

- **All TLE fields are publicly accessible** (no hidden fields like `jdsatepoch`)
- **All satellite data fields are publicly accessible**
- **Multiple propagation methods** for different use cases
- **Comprehensive convenience methods** for common operations
- **Clear documentation** of all accessible fields and methods

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Original SGP4/SDP4 implementation by Dr. TS Kelso
- Based on NORAD SGP4 orbital models
- Translation maintains fidelity to the original algorithms
