package sgp4go

import (
	"math"
)

// Sign returns the sign of a number (-1, 0, or 1)
func Sign(arg float64) int {
	if arg > 0 {
		return 1
	} else if arg < 0 {
		return -1
	}
	return 0
}

// Cube returns the cube of a number
func Cube(arg float64) float64 {
	return arg * arg * arg
}

// Square returns the square of a number
func Square(arg float64) float64 {
	return arg * arg
}

// Power returns arg raised to the power pwr
func Power(arg, pwr float64) float64 {
	if arg > 0 {
		return math.Exp(pwr * math.Log(arg))
	}
	// In the original Pascal code, this would print an error
	// For Go, we'll return NaN to indicate invalid operation
	return math.NaN()
}

// Radians converts degrees to radians
func Radians(arg float64) float64 {
	return arg * PI / 180.0
}

// Degrees converts radians to degrees
func Degrees(arg float64) float64 {
	return arg * 180.0 / PI
}

// Tan returns the tangent of an angle
func Tan(arg float64) float64 {
	return math.Sin(arg) / math.Cos(arg)
}

// ArcSin returns the arcsine of a number
func ArcSin(arg float64) float64 {
	if math.Abs(arg) >= 1.0 {
		return float64(Sign(arg)) * PI / 2.0
	}
	return math.Atan(arg / math.Sqrt(1.0-arg*arg))
}

// ArcCos returns the arccosine of a number
func ArcCos(arg float64) float64 {
	return PI/2.0 - ArcSin(arg)
}

// Modulus returns the remainder of arg1 divided by arg2
func Modulus(arg1, arg2 float64) float64 {
	return arg1 - arg2*math.Floor(arg1/arg2)
}

// Fmod2p returns the modulus of arg with respect to 2π
func Fmod2p(arg float64) float64 {
	return Modulus(arg, TWOPI)
}

// AcTan returns the arctangent of sinx/cosx, handling quadrant correctly
func AcTan(sinx, cosx float64) float64 {
	if cosx == 0.0 {
		if sinx > 0.0 {
			return PI / 2.0
		} else {
			return -PI / 2.0
		}
	}

	angle := math.Atan(sinx / cosx)
	if cosx < 0.0 {
		angle += PI
	} else if sinx < 0.0 {
		angle += TWOPI
	}

	return angle
}

// Magnitude calculates the magnitude of a vector
func Magnitude(v *Vector) {
	v.Magnitude = math.Sqrt(v.X*v.X + v.Y*v.Y + v.Z*v.Z)
}

// VecAdd adds two vectors and stores result in v3
func VecAdd(v1, v2 Vector) Vector {
	return Vector{
		X: v1.X + v2.X,
		Y: v1.Y + v2.Y,
		Z: v1.Z + v2.Z,
	}
}

// VecSub subtracts v2 from v1 and stores result in v3
func VecSub(v1, v2 Vector) Vector {
	return Vector{
		X: v1.X - v2.X,
		Y: v1.Y - v2.Y,
		Z: v1.Z - v2.Z,
	}
}

// ScalarMultiply multiplies a vector by a scalar
func ScalarMultiply(k float64, v Vector) Vector {
	return Vector{
		X: k * v.X,
		Y: k * v.Y,
		Z: k * v.Z,
	}
}

// Dot returns the dot product of two vectors
func Dot(v1, v2 Vector) float64 {
	return v1.X*v2.X + v1.Y*v2.Y + v1.Z*v2.Z
}

// Angle returns the angle between two vectors in radians
func Angle(v1, v2 Vector) float64 {
	dot := Dot(v1, v2)
	mag1 := math.Sqrt(v1.X*v1.X + v1.Y*v1.Y + v1.Z*v1.Z)
	mag2 := math.Sqrt(v2.X*v2.X + v2.Y*v2.Y + v2.Z*v2.Z)

	if mag1 == 0.0 || mag2 == 0.0 {
		return 0.0
	}

	cosAngle := dot / (mag1 * mag2)
	if cosAngle > 1.0 {
		cosAngle = 1.0
	} else if cosAngle < -1.0 {
		cosAngle = -1.0
	}

	return math.Acos(cosAngle)
}

// Cross returns the cross product of two vectors
func Cross(v1, v2 Vector) Vector {
	return Vector{
		X: v1.Y*v2.Z - v1.Z*v2.Y,
		Y: v1.Z*v2.X - v1.X*v2.Z,
		Z: v1.X*v2.Y - v1.Y*v2.X,
	}
}

// Normalize normalizes a vector to unit length
func Normalize(v *Vector) {
	mag := math.Sqrt(v.X*v.X + v.Y*v.Y + v.Z*v.Z)
	if mag > 0.0 {
		v.X /= mag
		v.Y /= mag
		v.Z /= mag
		v.Magnitude = 1.0
	}
}

// Min returns the minimum of two integers
func Min(arg1, arg2 int) int {
	if arg1 < arg2 {
		return arg1
	}
	return arg2
}

// Max returns the maximum of two integers
func Max(arg1, arg2 int) int {
	if arg1 > arg2 {
		return arg1
	}
	return arg2
}

// RMin returns the minimum of two floats
func RMin(arg1, arg2 float64) float64 {
	if arg1 < arg2 {
		return arg1
	}
	return arg2
}

// RMax returns the maximum of two floats
func RMax(arg1, arg2 float64) float64 {
	if arg1 > arg2 {
		return arg1
	}
	return arg2
}

// ConvertPositionToKilometers converts position from Earth radii to kilometers
func ConvertPositionToKilometers(pos *Vector) {
	pos.X *= XKMPER
	pos.Y *= XKMPER
	pos.Z *= XKMPER
	Magnitude(pos)
}

// ConvertVelocityToKilometersPerSecond converts velocity from Earth radii/minute to kilometers/second
func ConvertVelocityToKilometersPerSecond(vel *Vector) {
	vel.X *= XKMPER / 60.0
	vel.Y *= XKMPER / 60.0
	vel.Z *= XKMPER / 60.0
	Magnitude(vel)
}

// ConvertPositionAndVelocityToKilometers converts both position and velocity to kilometers
func ConvertPositionAndVelocityToKilometers(pos, vel *Vector) {
	ConvertPositionToKilometers(pos)
	ConvertVelocityToKilometersPerSecond(vel)
}
