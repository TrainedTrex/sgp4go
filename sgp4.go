package sgp4go

import (
	"fmt"
	"math"
)

// SGP4 propagates a satellite's position and velocity using the SGP4 model
// tsince is the time since epoch in minutes
// sat is the satellite data
// Returns position and velocity vectors in Earth radii and Earth radii/minute
func SGP4(tsince float64, sat *SatelliteData) (Vector, Vector, error) {
	if sat == nil {
		return Vector{}, Vector{}, fmt.Errorf("satellite data cannot be nil")
	}

	// Initialize derived constants if needed
	if sat.IFlag == 0 {
		defineDerivedConstants(sat)
	}

	// Call the appropriate model based on deep space flag
	if sat.IDeep == 1 {
		return SDP4(tsince, sat)
	} else {
		return SGP4NearEarth(tsince, sat)
	}
}

// SGP4NearEarth implements the near-earth SGP4 model
func SGP4NearEarth(tsince float64, sat *SatelliteData) (Vector, Vector, error) {
	// Local variables for SGP4 calculations
	var (
		a1, a3ovk2, ao, aodp, aycof, betao, betao2, c1, c1sq, c2, c3, c4, c5, coef, coef1           float64
		cosio, d2, d3, d4, del1, delmo, delo, eeta, eosq, eta, etasq                                float64
		omgcof, omgdot, perige, pinvsq, psisq, qoms24, s4, sinio, sinmo, t2cof, t3cof, t4cof, t5cof float64
		temp, temp1, temp2, temp3, theta2, theta4, tsi, x1m5th, x1mth2, x3thm1, x7thm1              float64
		xhdot1, xlcof, xmcof, xmdot, xnodcf, xnodot                                                 float64
		isimp                                                                                       int
	)

	// Recover original mean motion (xnodp) and semimajor axis (aodp) from input elements
	if sat.IFlag == 0 {
		a1 = Power(sat.XKE/sat.XNO, TOTHRD)
		cosio = math.Cos(sat.XIncl)
		theta2 = cosio * cosio
		x3thm1 = 3*theta2 - 1
		eosq = sat.EO * sat.EO
		betao2 = 1 - eosq
		betao = math.Sqrt(betao2)
		del1 = 1.5 * CK2 * x3thm1 / (a1 * a1 * betao * betao2)
		ao = a1 * (1 - del1*(0.5*TOTHRD+del1*(1+134.0/81.0*del1)))
		delo = 1.5 * CK2 * x3thm1 / (ao * ao * betao * betao2)
		sat.XNODP = sat.XNO / (1 + delo)
		aodp = ao / (1 - delo)

		// For perigee less than 220 kilometers, the isimp flag is set
		isimp = 0
		if (aodp*(1-sat.EO)-AE)*XKMPER < 220.0 {
			isimp = 1
		}

		// For perigee below 156 km, the values of s and qoms2t are altered
		s4 = S
		qoms24 = sat.QOMS2T
		perige = (aodp*(1-sat.EO) - AE) * XKMPER
		if perige >= 156.0 {
			// Continue with normal processing
		} else {
			s4 = perige - 78.0
			if perige > 98.0 {
				// Continue with s4 = perige - 78
			} else {
				s4 = 20.0
			}
			qoms24 = Power((120.0-s4)*AE/XKMPER, 4)
			s4 = s4/XKMPER + AE
		}

		// Calculate other derived constants
		pinvsq = 1.0 / (aodp * aodp * betao2 * betao2)
		tsi = 1.0 / (aodp - s4)
		eta = aodp * sat.EO * tsi
		etasq = eta * eta
		eeta = sat.EO * eta
		psisq = math.Abs(1.0 - etasq)
		coef = qoms24 * Power(tsi, 4)
		coef1 = coef / Power(psisq, 3.5)
		c2 = coef1 * sat.XNODP * (aodp*(1+1.5*etasq+eeta*(4+etasq)) + 0.75*CK2*tsi/psisq*x3thm1*(8+3*etasq*(8+etasq)))
		c1 = sat.BStar * c2
		sinio = math.Sin(sat.XIncl)
		a3ovk2 = -XJ3 / CK2 * Power(AE, 3)
		c3 = coef * tsi * a3ovk2 * sat.XNODP * AE * sinio / sat.EO
		x1mth2 = 1 - theta2
		c4 = 2 * sat.XNODP * coef1 * aodp * betao2 * (eta*(2+0.5*etasq) + sat.EO*(0.5+2*etasq) - 2*CK2*tsi/(aodp*psisq)*(-3*x3thm1*(1-2*eeta+etasq*(1.5-0.5*eeta))+0.75*x1mth2*(2*etasq-eeta*(1+etasq))*math.Cos(2*sat.OmegaO)))
		c5 = 2 * coef1 * aodp * betao2 * (1 + 2.75*(etasq+eeta) + eeta*etasq)
		theta4 = theta2 * theta2
		temp1 = 3 * CK2 * pinvsq * sat.XNODP
		temp2 = temp1 * CK2 * pinvsq
		temp3 = 1.25 * CK4 * pinvsq * pinvsq * sat.XNODP
		xmdot = sat.XNODP + 0.5*temp1*betao*x3thm1 + 0.0625*temp2*betao*(13-78*theta2+137*theta4)
		x1m5th = 1 - 5*theta2
		omgdot = -0.5*temp1*x1m5th + 0.0625*temp2*(7-114*theta2+395*theta4) + temp3*(3-36*theta2+49*theta4)
		xhdot1 = -temp1 * cosio
		xnodot = xhdot1 + (0.5*temp2*(4-19*theta2)+2*temp3*(3-7*theta2))*cosio
		omgcof = sat.BStar * c3 * math.Cos(sat.OmegaO)
		xmcof = -TOTHRD * coef * sat.BStar * AE / eeta
		xnodcf = 3.5 * betao2 * xhdot1 * c1
		t2cof = 1.5 * c1
		xlcof = 0.125 * a3ovk2 * sinio * (3 + 5*cosio) / (1 + cosio)
		aycof = 0.25 * a3ovk2 * sinio
		delmo = Power(1+eta*math.Cos(sat.XMO), 3)
		sinmo = math.Sin(sat.XMO)
		x7thm1 = 7*theta2 - 1

		if isimp == 1 {
			// Simplified model - skip complex coefficient calculations
			// This corresponds to the goto 90 in the original code
		} else {
			c1sq = c1 * c1
			d2 = 4 * aodp * tsi * c1sq
			temp = d2 * tsi * c1 / 3
			d3 = (17*aodp + s4) * temp
			d4 = 0.5 * temp * aodp * tsi * (221*aodp + 31*s4) * c1
			t3cof = d2 + 2*c1sq
			t4cof = 0.25 * (3*d3 + c1*(12*d2+10*c1sq))
			t5cof = 0.2 * (3*d4 + 12*c1*d3 + 6*d2*d2 + 15*c1sq*(2*d2+c1sq))
		}

		sat.IFlag = 0
	}

	// Update for secular gravity and atmospheric drag
	xmdf := sat.XMO + xmdot*tsince
	omgadf := sat.OmegaO + omgdot*tsince
	xnoddf := sat.XNodeO + xnodot*tsince
	omega := omgadf
	xmp := xmdf
	tsq := tsince * tsince
	xnode := xnoddf + xnodcf*tsq
	tempa := 1 - c1*tsince
	tempe := sat.BStar * c4 * tsince
	templ := t2cof * tsq

	if isimp == 1 {
		// Simplified model - skip complex secular updates
		// This corresponds to the goto 110 in the original code
		// For low perigee satellites, we use linear variation in sqrt a
		// and quadratic variation in mean anomaly, dropping c3, delta omega,
		// and delta m terms
	} else {
		delomg := omgcof * tsince
		delm := xmcof * (Power(1+eta*math.Cos(xmdf), 3) - delmo)
		temp := delomg + delm
		xmp = xmdf + temp
		omega = omgadf - temp
		tcube := tsq * tsince
		tfour := tsince * tcube
		tempa = tempa - d2*tsq - d3*tcube - d4*tfour
		tempe = tempe + sat.BStar*c5*(math.Sin(xmp)-sinmo)
		templ = templ + t3cof*tcube + tfour*(t4cof+tsince*t5cof)
	}

	a := aodp * Power(tempa, 2)
	e := sat.EO - tempe
	xl := xmp + omega + xnode + sat.XNODP*templ
	beta := math.Sqrt(1 - e*e)
	xn := sat.XKE / Power(a, 1.5)

	// Long period periodics
	axn := e * math.Cos(omega)
	longTemp := 1.0 / (a * beta * beta)
	xll := longTemp * xlcof * axn
	aynl := longTemp * aycof
	xlt := xl + xll
	ayn := e*math.Sin(omega) + aynl

	// Solve Kepler's Equation
	capu := Fmod2p(xlt - xnode)
	keplerTemp2 := capu
	var epw, keplerTemp3, keplerTemp4, keplerTemp5, keplerTemp6, sinepw, cosepw float64
	for i := 0; i < 10; i++ {
		sinepw = math.Sin(keplerTemp2)
		cosepw = math.Cos(keplerTemp2)
		keplerTemp3 = axn * sinepw
		keplerTemp4 = ayn * cosepw
		keplerTemp5 = axn * cosepw
		keplerTemp6 = ayn * sinepw
		epw = (capu-keplerTemp4+keplerTemp3-keplerTemp2)/(1-keplerTemp5-keplerTemp6) + keplerTemp2
		if math.Abs(epw-keplerTemp2) <= E6A {
			break
		}
		keplerTemp2 = epw
	}

	// Short period preliminary quantities
	ecose := keplerTemp5 + keplerTemp6
	esine := keplerTemp3 - keplerTemp4
	elsq := axn*axn + ayn*ayn
	shortTemp := 1 - elsq
	pl := a * shortTemp
	r := a * (1 - ecose)
	shortTemp1 := 1.0 / r
	rdot := sat.XKE * math.Sqrt(a) * esine * shortTemp1
	rfdot := sat.XKE * math.Sqrt(pl) * shortTemp1
	shortTemp2 := a * shortTemp1
	betal := math.Sqrt(shortTemp)
	shortTemp3 := 1.0 / (1 + betal)
	cosu := shortTemp2 * (cosepw - axn + ayn*esine*shortTemp3)
	sinu := shortTemp2 * (sinepw - ayn - axn*esine*shortTemp3)
	u := AcTan(sinu, cosu)
	sin2u := 2 * sinu * cosu
	cos2u := 2*cosu*cosu - 1
	shortTemp = 1.0 / pl
	shortTemp1 = CK2 * shortTemp
	shortTemp2 = shortTemp1 * shortTemp

	// Update for short periodics
	rk := r*(1-1.5*shortTemp2*betal*x3thm1) + 0.5*shortTemp1*x1mth2*cos2u
	uk := u - 0.25*shortTemp2*x7thm1*sin2u
	xnodek := xnode + 1.5*shortTemp2*cosio*sin2u
	xinck := sat.XIncl + 1.5*shortTemp2*cosio*sinio*cos2u
	rdotk := rdot - xn*shortTemp1*x1mth2*sin2u
	rfdotk := rfdot + xn*shortTemp1*(x1mth2*cos2u+1.5*x3thm1)

	// Orientation vectors
	sinuk := math.Sin(uk)
	cosuk := math.Cos(uk)
	sinik := math.Sin(xinck)
	cosik := math.Cos(xinck)
	sinnok := math.Sin(xnodek)
	cosnok := math.Cos(xnodek)
	xmx := -sinnok * cosik
	xmy := cosnok * cosik
	ux := xmx*sinuk + cosnok*cosuk
	uy := xmy*sinuk + sinnok*cosuk
	uz := sinik * sinuk
	vx := xmx*cosuk - cosnok*sinuk
	vy := xmy*cosuk - sinnok*sinuk
	vz := sinik * cosuk

	// Position and velocity
	pos := Vector{
		X: rk * ux,
		Y: rk * uy,
		Z: rk * uz,
	}
	vel := Vector{
		X: rdotk*ux + rfdotk*vx,
		Y: rdotk*uy + rfdotk*vy,
		Z: rdotk*uz + rfdotk*vz,
	}

	Magnitude(&pos)
	Magnitude(&vel)

	return pos, vel, nil
}

// SDP4 implements the deep-space SDP4 model
func SDP4(tsince float64, sat *SatelliteData) (Vector, Vector, error) {
	// This is a placeholder for the deep-space model
	// The full SDP4 implementation is quite complex and would require
	// significant additional code. For now, keep it close to earth, this returns an error
	return Vector{}, Vector{}, fmt.Errorf("SDP4 deep-space model not yet implemented")
}

func defineDerivedConstants(sat *SatelliteData) {
	sat.XKE = math.Sqrt(3600.0 * GE / Cube(XKMPER))
	sat.QOMS2T = Power(QO-S, 4)
}

// SGP is the main interface function that automatically selects SGP4 or SDP4
func SGP(time float64, sat *SatelliteData) (Vector, Vector, error) {
	if sat == nil {
		return Vector{}, Vector{}, fmt.Errorf("satellite data cannot be nil")
	}

	// Calculate time since epoch in minutes
	tsince := (time - sat.JulianEpoch) * 1440.0

	// Call the appropriate model
	return SGP4(tsince, sat)
}

// PropagateSatellite is a high-level convenience function that takes a TLE and time
// and returns the satellite's position and velocity at that time
// units can be "km" for kilometers or any other value for Earth radii (default)
func PropagateSatellite(tle *TLE, time float64, units ...string) (Vector, Vector, error) {
	if tle == nil {
		return Vector{}, Vector{}, fmt.Errorf("TLE cannot be nil")
	}

	sat, err := ConvertSatelliteData(tle)
	if err != nil {
		return Vector{}, Vector{}, fmt.Errorf("failed to convert TLE: %v", err)
	}

	pos, vel, err := SGP(time, sat)
	if err != nil {
		return Vector{}, Vector{}, err
	}

	if len(units) > 0 && units[0] == "km" {
		ConvertPositionAndVelocityToKilometers(&pos, &vel)
	}

	return pos, vel, nil
}
