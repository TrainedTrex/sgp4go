package sgp4go

// Physical constants from WGS-72 model
const (
	// Earth constants
	AE     = 1.0          // Earth equatorial radius in Earth radii
	TOTHRD = 2.0 / 3.0    // 2/3
	XKMPER = 6378.135     // Earth equatorial radius - kilometers (WGS '72)
	F      = 1.0 / 298.26 // Earth flattening (WGS '72)
	GE     = 398600.8     // Earth gravitational constant (WGS '72)

	// Harmonic coefficients
	J2 = 1.0826158e-3 // J2 harmonic (WGS '72)
	J3 = -2.53881e-6  // J3 harmonic (WGS '72)
	J4 = -1.65597e-6  // J4 harmonic (WGS '72)

	// Derived constants
	CK2 = J2 / 2.0        // J2/2
	CK4 = -3.0 * J4 / 8.0 // -3*J4/8
	XJ3 = J3              // J3

	// Atmospheric constants
	QO = AE + 120.0/XKMPER // Atmospheric boundary
	S  = AE + 78.0/XKMPER  // Atmospheric boundary

	// Numerical constants
	E6A = 1e-6 // Small number for numerical stability

	// Deep space model flags
	DPINIT = 1 // Deep-space initialization code
	DPSEC  = 2 // Deep-space secular code
	DPPER  = 3 // Deep-space periodic code

	// Time constants
	TWOPI  = 2.0 * PI // 2π
	XMNPDA = 1440.0   // Minutes per day
)

// Mathematical constants
const (
	PI = 3.14159265358979323846 // π
)

// Maximum number of satellites that can be stored
const MAX_SATS = 250
