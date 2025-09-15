package sgp4go

import (
	"fmt"
	"math"
)

// ElsetRec represents the satellite record structure from the C code
// This is an exact replica of the C ElsetRec struct from SGP4.h
type ElsetRec struct {
	// Basic satellite information
	WhichConst, EpochYr, EpochTynumrev, Error int
	SatID                                     [6]byte
	OperationMode, Init, Method               byte

	// Orbital elements
	A, Altp, Alta, EpochDays, Jdsatepoch, JdsatepochF float64
	Nddot, Ndot, Bstar, Rcse, Inclo, Nodeo, Ecco      float64
	Argpo, Mo, NoKozai                                float64

	// Additional TLE fields
	Classification byte
	Intldesg       [12]byte
	Ephtype        int
	Elnum, Revnum  int64

	// Unkozai'd variable
	NoUnkozai float64

	// Singly averaged variables
	Am, Em, Im, Om, Om2, Mm, Nm, T float64

	// Constant parameters
	Tumin, Mu, Radiusearthkm, Xke, J2, J3, J4, J3oj2 float64

	// Additional elements
	DiaMm              int64
	PeriodSec, RcsM2   float64
	Active, NotOrbital byte

	// Temporary variables
	Ep, Inclp, Nodep, Argpp, Mp float64

	// Near earth variables
	Isimp                                      int
	Aycof, Con41, Cc1, Cc4, Cc5, D2, D3, D4    float64
	Delmo, Eta, Argpdot, Omgcof, Sinmao        float64
	T2cof, T3cof, T4cof, T5cof, X1mth2, X7thm1 float64
	Mdot, Nodedot, Xlcof, Xmcof, Nodecf        float64

	// Deep space variables
	Irez                                                                 int
	D2201, D2211, D3210, D3222, D4410, D4422, D5220, D5232, D5421, D5433 float64
	Dedt, Del1, Del2, Del3, Didt, Dmdt, Dnodt, Domdt                     float64
	E3, Ee2, Peo, Pgho, Pho, Pinco, Plo, Se2, Se3, Sgh2, Sgh3, Sgh4      float64
	Sh2, Sh3, Si2, Si3, Sl2, Sl3, Sl4, Gsto, Xfact                       float64
	Xgh2, Xgh3, Xgh4, Xh2, Xh3, Xi2, Xi3, Xl2, Xl3, Xl4                  float64
	Xlamo, Zmol, Zmos, Atime, Xli, Xni                                   float64
	Snodm, Cnodm, Sinim, Cosim, Sinomm, Cosomm                           float64
	Day, Emsq, Gam, Rtemsq                                               float64
	S1, S2, S3, S4, S5, S6, S7                                           float64
	Ss1, Ss2, Ss3, Ss4, Ss5, Ss6, Ss7                                    float64
	Sz1, Sz2, Sz3, Sz11, Sz12, Sz13, Sz21, Sz22, Sz23, Sz31, Sz32, Sz33  float64
	Z1, Z2, Z3, Z11, Z12, Z13, Z21, Z22, Z23, Z31, Z32, Z33              float64
	Argpm, Inclm, Nodem, Dndt, Eccsq                                     float64

	// For initl
	Ainv, Ao, Con42, Cosio, Cosio2, Omeosq, Posq, Rp, Rteosq, Sinio float64
}

// Constants from the C code
const (
	WGS72OLD = 1
	WGS72    = 2
	WGS84    = 3
	PI_C     = 3.14159265358979323846
	TWOPI_C  = 2.0 * PI_C
	DEG2RAD  = PI_C / 180.0
)

// GetGravConst sets gravitational constants based on whichconst
func GetGravConst(whichconst int, rec *ElsetRec) {
	rec.WhichConst = whichconst
	switch whichconst {
	case WGS72OLD:
		rec.Mu = 398600.79964
		rec.Radiusearthkm = 6378.135
		rec.Xke = 0.0743669161
		rec.Tumin = 1.0 / rec.Xke
		rec.J2 = 0.001082616
		rec.J3 = -0.00000253881
		rec.J4 = -0.00000165597
		rec.J3oj2 = rec.J3 / rec.J2
	case WGS72:
		rec.Mu = 398600.8
		rec.Radiusearthkm = 6378.135
		rec.Xke = 60.0 / math.Sqrt(rec.Radiusearthkm*rec.Radiusearthkm*rec.Radiusearthkm/rec.Mu)
		rec.Tumin = 1.0 / rec.Xke
		rec.J2 = 0.001082616
		rec.J3 = -0.00000253881
		rec.J4 = -0.00000165597
		rec.J3oj2 = rec.J3 / rec.J2
	default: // WGS84
		rec.Mu = 398600.5
		rec.Radiusearthkm = 6378.137
		rec.Xke = 60.0 / math.Sqrt(rec.Radiusearthkm*rec.Radiusearthkm*rec.Radiusearthkm/rec.Mu)
		rec.Tumin = 1.0 / rec.Xke
		rec.J2 = 0.00108262998905
		rec.J3 = -0.00000253215306
		rec.J4 = -0.00000161098761
		rec.J3oj2 = rec.J3 / rec.J2
	}
}

// InitL initializes the SGP4 propagator - exact replica of initl from C code
func InitL(epoch float64, rec *ElsetRec) {
	// Local variables
	var ak, d1, del, adel, po, x2o3 float64
	var ds70, ts70, tfrac, c1, thgr70, fk5r, c1p2p float64

	x2o3 = 2.0 / 3.0

	// Calculate auxiliary epoch quantities
	rec.Eccsq = rec.Ecco * rec.Ecco
	rec.Omeosq = 1.0 - rec.Eccsq
	rec.Rteosq = math.Sqrt(rec.Omeosq)
	rec.Cosio = math.Cos(rec.Inclo)
	rec.Cosio2 = rec.Cosio * rec.Cosio

	// Un-kozai the mean motion
	ak = math.Pow(rec.Xke/rec.NoKozai, x2o3)
	d1 = 0.75 * rec.J2 * (3.0*rec.Cosio2 - 1.0) / (rec.Rteosq * rec.Omeosq)
	del = d1 / (ak * ak)
	adel = ak * (1.0 - del*del - del*(1.0/3.0+134.0*del*del/81.0))
	del = d1 / (adel * adel)
	rec.NoUnkozai = rec.NoKozai / (1.0 + del)

	rec.Ao = math.Pow(rec.Xke/(rec.NoUnkozai), x2o3)
	rec.Sinio = math.Sin(rec.Inclo)
	po = rec.Ao * rec.Omeosq
	rec.Con42 = 1.0 - 5.0*rec.Cosio2
	rec.Con41 = -rec.Con42 - rec.Cosio2 - rec.Cosio2
	rec.Ainv = 1.0 / rec.Ao
	rec.Posq = po * po
	rec.Rp = rec.Ao * (1.0 - rec.Ecco)
	rec.Method = 'n'

	// Modern approach to finding sidereal time
	ts70 = epoch - 7305.0
	ds70 = math.Floor(ts70 + 1.0e-8)
	tfrac = ts70 - ds70
	c1 = 1.72027916940703639e-2
	thgr70 = 1.7321343856509374
	fk5r = 5.07551419432269442e-15
	c1p2p = c1 + TWOPI_C
	_ = math.Mod(thgr70+c1*ds70+c1p2p*tfrac+ts70*ts70*fk5r, TWOPI_C)
	rec.Gsto = Gstime(epoch + 2433281.5)
}

// Gstime finds the Greenwich sidereal time
func Gstime(jdut1 float64) float64 {
	var temp, tut1 float64

	tut1 = (jdut1 - 2451545.0) / 36525.0
	temp = -6.2e-6*tut1*tut1*tut1 + 0.093104*tut1*tut1 +
		(876600.0*3600+8640184.812866)*tut1 + 67310.54841
	temp = math.Mod(temp*DEG2RAD/240.0, TWOPI_C)

	if temp < 0.0 {
		temp += TWOPI_C
	}

	return temp
}

// SGP4Init initializes the SGP4 propagator - exact replica of sgp4init from C code
func SGP4Init(opsmode byte, satrec *ElsetRec) bool {
	// Local variables
	var cc1sq, cc2, cc3, coef, coef1, cosio4 float64
	var eeta, etasq, perige, pinvsq, psisq, qzms24 float64
	var sfour, tc, temp, temp1, temp2, temp3, tsi, xpidot float64
	var xhdot1, qzms2t, ss, x2o3 float64
	var r, v [3]float64
	var delmotemp, qzms2ttemp, qzms24temp float64

	epoch := (satrec.Jdsatepoch + satrec.JdsatepochF) - 2433281.5

	// Initialization
	const temp4 = 1.5e-12

	// Set all near earth variables to zero
	satrec.Isimp = 0
	satrec.Method = 'n'
	satrec.Aycof = 0.0
	satrec.Con41 = 0.0
	satrec.Cc1 = 0.0
	satrec.Cc4 = 0.0
	satrec.Cc5 = 0.0
	satrec.D2 = 0.0
	satrec.D3 = 0.0
	satrec.D4 = 0.0
	satrec.Delmo = 0.0
	satrec.Eta = 0.0
	satrec.Argpdot = 0.0
	satrec.Omgcof = 0.0
	satrec.Sinmao = 0.0
	satrec.T = 0.0
	satrec.T2cof = 0.0
	satrec.T3cof = 0.0
	satrec.T4cof = 0.0
	satrec.T5cof = 0.0
	satrec.X1mth2 = 0.0
	satrec.X7thm1 = 0.0
	satrec.Mdot = 0.0
	satrec.Nodedot = 0.0
	satrec.Xlcof = 0.0
	satrec.Xmcof = 0.0
	satrec.Nodecf = 0.0

	// Set all deep space variables to zero
	satrec.Irez = 0
	satrec.D2201 = 0.0
	satrec.D2211 = 0.0
	satrec.D3210 = 0.0
	satrec.D3222 = 0.0
	satrec.D4410 = 0.0
	satrec.D4422 = 0.0
	satrec.D5220 = 0.0
	satrec.D5232 = 0.0
	satrec.D5421 = 0.0
	satrec.D5433 = 0.0
	satrec.Dedt = 0.0
	satrec.Del1 = 0.0
	satrec.Del2 = 0.0
	satrec.Del3 = 0.0
	satrec.Didt = 0.0
	satrec.Dmdt = 0.0
	satrec.Dnodt = 0.0
	satrec.Domdt = 0.0
	satrec.E3 = 0.0
	satrec.Ee2 = 0.0
	satrec.Peo = 0.0
	satrec.Pgho = 0.0
	satrec.Pho = 0.0
	satrec.Pinco = 0.0
	satrec.Plo = 0.0
	satrec.Se2 = 0.0
	satrec.Se3 = 0.0
	satrec.Sgh2 = 0.0
	satrec.Sgh3 = 0.0
	satrec.Sgh4 = 0.0
	satrec.Sh2 = 0.0
	satrec.Sh3 = 0.0
	satrec.Si2 = 0.0
	satrec.Si3 = 0.0
	satrec.Sl2 = 0.0
	satrec.Sl3 = 0.0
	satrec.Sl4 = 0.0
	satrec.Gsto = 0.0
	satrec.Xfact = 0.0
	satrec.Xgh2 = 0.0
	satrec.Xgh3 = 0.0
	satrec.Xgh4 = 0.0
	satrec.Xh2 = 0.0
	satrec.Xh3 = 0.0
	satrec.Xi2 = 0.0
	satrec.Xi3 = 0.0
	satrec.Xl2 = 0.0
	satrec.Xl3 = 0.0
	satrec.Xl4 = 0.0
	satrec.Xlamo = 0.0
	satrec.Zmol = 0.0
	satrec.Zmos = 0.0
	satrec.Atime = 0.0
	satrec.Xli = 0.0
	satrec.Xni = 0.0

	// Earth constants
	GetGravConst(satrec.WhichConst, satrec)

	satrec.Error = 0
	satrec.OperationMode = opsmode

	// Single averaged mean elements
	satrec.Am = 0.0
	satrec.Em = 0.0
	satrec.Im = 0.0
	satrec.Om = 0.0
	satrec.Om2 = 0.0
	satrec.Mm = 0.0
	satrec.Nm = 0.0

	// Earth constants
	ss = 78.0/satrec.Radiusearthkm + 1.0
	qzms2ttemp = (120.0 - 78.0) / satrec.Radiusearthkm
	qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp
	x2o3 = 2.0 / 3.0

	satrec.Init = 'y'
	satrec.T = 0.0

	InitL(epoch, satrec)

	satrec.A = math.Pow(satrec.NoUnkozai*satrec.Tumin, -2.0/3.0)
	satrec.Altp = satrec.A*(1.0+satrec.Ecco) - 1.0
	satrec.Alta = satrec.A*(1.0-satrec.Ecco) - 1.0
	satrec.Error = 0

	if (satrec.Omeosq >= 0.0) || (satrec.NoUnkozai >= 0.0) {
		satrec.Isimp = 0
		if satrec.Rp < (220.0/satrec.Radiusearthkm + 1.0) {
			satrec.Isimp = 1
		}
		sfour = ss
		qzms24 = qzms2t
		perige = (satrec.Rp - 1.0) * satrec.Radiusearthkm

		// For perigees below 156 km, s and qoms2t are altered
		if perige < 156.0 {
			sfour = perige - 78.0
			if perige < 98.0 {
				sfour = 20.0
			}
			qzms24temp = (120.0 - sfour) / satrec.Radiusearthkm
			qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp
			sfour = sfour/satrec.Radiusearthkm + 1.0
		}
		pinvsq = 1.0 / satrec.Posq

		tsi = 1.0 / (satrec.Ao - sfour)
		satrec.Eta = satrec.Ao * satrec.Ecco * tsi
		etasq = satrec.Eta * satrec.Eta
		eeta = satrec.Ecco * satrec.Eta
		psisq = math.Abs(1.0 - etasq)
		coef = qzms24 * math.Pow(tsi, 4.0)
		coef1 = coef / math.Pow(psisq, 3.5)
		cc2 = coef1 * satrec.NoUnkozai * (satrec.Ao*(1.0+1.5*etasq+eeta*
			(4.0+etasq)) + 0.375*satrec.J2*tsi/psisq*satrec.Con41*
			(8.0+3.0*etasq*(8.0+etasq)))
		satrec.Cc1 = satrec.Bstar * cc2
		cc3 = 0.0
		if satrec.Ecco > 1.0e-4 {
			cc3 = -2.0 * coef * tsi * satrec.J3oj2 * satrec.NoUnkozai * satrec.Sinio / satrec.Ecco
		}
		satrec.X1mth2 = 1.0 - satrec.Cosio2
		satrec.Cc4 = 2.0 * satrec.NoUnkozai * coef1 * satrec.Ao * satrec.Omeosq *
			(satrec.Eta*(2.0+0.5*etasq) + satrec.Ecco*
				(0.5+2.0*etasq) - satrec.J2*tsi/(satrec.Ao*psisq)*
				(-3.0*satrec.Con41*(1.0-2.0*eeta+etasq*
					(1.5-0.5*eeta))+0.75*satrec.X1mth2*
					(2.0*etasq-eeta*(1.0+etasq))*math.Cos(2.0*satrec.Argpo)))
		satrec.Cc5 = 2.0 * coef1 * satrec.Ao * satrec.Omeosq * (1.0 + 2.75*
			(etasq+eeta) + eeta*etasq)
		cosio4 = satrec.Cosio2 * satrec.Cosio2
		temp1 = 1.5 * satrec.J2 * pinvsq * satrec.NoUnkozai
		temp2 = 0.5 * temp1 * satrec.J2 * pinvsq
		temp3 = -0.46875 * satrec.J4 * pinvsq * pinvsq * satrec.NoUnkozai
		satrec.Mdot = satrec.NoUnkozai + 0.5*temp1*satrec.Rteosq*satrec.Con41 + 0.0625*
			temp2*satrec.Rteosq*(13.0-78.0*satrec.Cosio2+137.0*cosio4)
		satrec.Argpdot = -0.5*temp1*satrec.Con42 + 0.0625*temp2*
			(7.0-114.0*satrec.Cosio2+395.0*cosio4) +
			temp3*(3.0-36.0*satrec.Cosio2+49.0*cosio4)
		xhdot1 = -temp1 * satrec.Cosio
		satrec.Nodedot = xhdot1 + (0.5*temp2*(4.0-19.0*satrec.Cosio2)+
			2.0*temp3*(3.0-7.0*satrec.Cosio2))*satrec.Cosio
		xpidot = satrec.Argpdot + satrec.Nodedot
		satrec.Omgcof = satrec.Bstar * cc3 * math.Cos(satrec.Argpo)
		satrec.Xmcof = 0.0
		if satrec.Ecco > 1.0e-4 {
			satrec.Xmcof = -x2o3 * coef * satrec.Bstar / eeta
		}
		satrec.Nodecf = 3.5 * satrec.Omeosq * xhdot1 * satrec.Cc1
		satrec.T2cof = 1.5 * satrec.Cc1
		// sgp4fix for divide by zero with xinco = 180 deg
		if math.Abs(satrec.Cosio+1.0) > 1.5e-12 {
			satrec.Xlcof = -0.25 * satrec.J3oj2 * satrec.Sinio * (3.0 + 5.0*satrec.Cosio) / (1.0 + satrec.Cosio)
		} else {
			satrec.Xlcof = -0.25 * satrec.J3oj2 * satrec.Sinio * (3.0 + 5.0*satrec.Cosio) / temp4
		}
		satrec.Aycof = -0.5 * satrec.J3oj2 * satrec.Sinio
		// sgp4fix use multiply for speed instead of pow
		delmotemp = 1.0 + satrec.Eta*math.Cos(satrec.Mo)
		satrec.Delmo = delmotemp * delmotemp * delmotemp
		satrec.Sinmao = math.Sin(satrec.Mo)
		satrec.X7thm1 = 7.0*satrec.Cosio2 - 1.0

		// Deep space initialization
		if (2 * PI_C / satrec.NoUnkozai) >= 225.0 {
			satrec.Method = 'd'
			satrec.Isimp = 1
			tc = 0.0
			satrec.Inclm = satrec.Inclo

			Dscom(epoch, satrec.Ecco, satrec.Argpo, tc, satrec.Inclo, satrec.Nodeo, satrec.NoUnkozai, satrec)

			satrec.Ep = satrec.Ecco
			satrec.Inclp = satrec.Inclo
			satrec.Nodep = satrec.Nodeo
			satrec.Argpp = satrec.Argpo
			satrec.Mp = satrec.Mo

			Dpper(satrec.E3, satrec.Ee2, satrec.Peo, satrec.Pgho,
				satrec.Pho, satrec.Pinco, satrec.Plo, satrec.Se2,
				satrec.Se3, satrec.Sgh2, satrec.Sgh3, satrec.Sgh4,
				satrec.Sh2, satrec.Sh3, satrec.Si2, satrec.Si3,
				satrec.Sl2, satrec.Sl3, satrec.Sl4, satrec.T,
				satrec.Xgh2, satrec.Xgh3, satrec.Xgh4, satrec.Xh2,
				satrec.Xh3, satrec.Xi2, satrec.Xi3, satrec.Xl2,
				satrec.Xl3, satrec.Xl4, satrec.Zmol, satrec.Zmos, satrec.Init, satrec,
				satrec.OperationMode)

			satrec.Ecco = satrec.Ep
			satrec.Inclo = satrec.Inclp
			satrec.Nodeo = satrec.Nodep
			satrec.Argpo = satrec.Argpp
			satrec.Mo = satrec.Mp

			satrec.Argpm = 0.0
			satrec.Nodem = 0.0
			satrec.Mm = 0.0

			Dsinit(tc, xpidot, satrec)
		}

		// Set variables if not deep space
		if satrec.Isimp != 1 {
			cc1sq = satrec.Cc1 * satrec.Cc1
			satrec.D2 = 4.0 * satrec.Ao * tsi * cc1sq
			temp = satrec.D2 * tsi * satrec.Cc1 / 3.0
			satrec.D3 = (17.0*satrec.Ao + sfour) * temp
			satrec.D4 = 0.5 * temp * satrec.Ao * tsi * (221.0*satrec.Ao + 31.0*sfour) * satrec.Cc1
			satrec.T3cof = satrec.D2 + 2.0*cc1sq
			satrec.T4cof = 0.25 * (3.0*satrec.D3 + satrec.Cc1*
				(12.0*satrec.D2+10.0*cc1sq))
			satrec.T5cof = 0.2 * (3.0*satrec.D4 +
				12.0*satrec.Cc1*satrec.D3 +
				6.0*satrec.D2*satrec.D2 +
				15.0*cc1sq*(2.0*satrec.D2+cc1sq))
		}
	}

	// Finally propagate to zero epoch to initialize all others
	SGP4(satrec, 0.0, r[:], v[:])

	satrec.Init = 'n'

	return true
}

// Dscom provides deep space common items - exact replica from C code
func Dscom(epoch, ep, argpp, tc, inclp, nodep, np float64, rec *ElsetRec) {
	// Constants
	const zes = 0.01675
	const zel = 0.05490
	const c1ss = 2.9864797e-6
	const c1l = 4.7968065e-7
	const zsinis = 0.39785416
	const zcosis = 0.91744867
	const zcosgs = 0.1945905
	const zsings = -0.98088458

	// Local variables
	var lsflg int
	var a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, betasq, cc, ctem, stem float64
	var x1, x2, x3, x4, x5, x6, x7, x8, xnodce, xnoi float64
	var zcosg, zcosgl, zcosh, zcoshl, zcosi, zcosil, zsing, zsingl, zsinh, zsinhl, zsini, zsinil, zx, zy float64

	rec.Nm = np
	rec.Em = ep
	rec.Snodm = math.Sin(nodep)
	rec.Cnodm = math.Cos(nodep)
	rec.Sinomm = math.Sin(argpp)
	rec.Cosomm = math.Cos(argpp)
	rec.Sinim = math.Sin(inclp)
	rec.Cosim = math.Cos(inclp)
	rec.Emsq = rec.Em * rec.Em
	betasq = 1.0 - rec.Emsq
	rec.Rtemsq = math.Sqrt(betasq)

	// Initialize lunar solar terms
	rec.Peo = 0.0
	rec.Pinco = 0.0
	rec.Plo = 0.0
	rec.Pgho = 0.0
	rec.Pho = 0.0
	rec.Day = epoch + 18261.5 + tc/1440.0
	xnodce = math.Mod(4.5236020-9.2422029e-4*rec.Day, TWOPI_C)
	stem = math.Sin(xnodce)
	ctem = math.Cos(xnodce)
	zcosil = 0.91375164 - 0.03568096*ctem
	zsinil = math.Sqrt(1.0 - zcosil*zcosil)
	zsinhl = 0.089683511 * stem / zsinil
	zcoshl = math.Sqrt(1.0 - zsinhl*zsinhl)
	rec.Gam = 5.8351514 + 0.0019443680*rec.Day
	zx = 0.39785416 * stem / zsinil
	zy = zcoshl*ctem + 0.91744867*zsinhl*stem
	zx = math.Atan2(zx, zy)
	zx = rec.Gam + zx - xnodce
	zcosgl = math.Cos(zx)
	zsingl = math.Sin(zx)

	// Do solar terms
	zcosg = zcosgs
	zsing = zsings
	zcosi = zcosis
	zsini = zsinis
	zcosh = rec.Cnodm
	zsinh = rec.Snodm
	cc = c1ss
	xnoi = 1.0 / rec.Nm

	for lsflg = 1; lsflg <= 2; lsflg++ {
		a1 = zcosg*zcosh + zsing*zcosi*zsinh
		a3 = -zsing*zcosh + zcosg*zcosi*zsinh
		a7 = -zcosg*zsinh + zsing*zcosi*zcosh
		a8 = zsing * zsini
		a9 = zsing*zsinh + zcosg*zcosi*zcosh
		a10 = zcosg * zsini
		a2 = rec.Cosim*a7 + rec.Sinim*a8
		a4 = rec.Cosim*a9 + rec.Sinim*a10
		a5 = -rec.Sinim*a7 + rec.Cosim*a8
		a6 = -rec.Sinim*a9 + rec.Cosim*a10

		x1 = a1*rec.Cosomm + a2*rec.Sinomm
		x2 = a3*rec.Cosomm + a4*rec.Sinomm
		x3 = -a1*rec.Sinomm + a2*rec.Cosomm
		x4 = -a3*rec.Sinomm + a4*rec.Cosomm
		x5 = a5 * rec.Sinomm
		x6 = a6 * rec.Sinomm
		x7 = a5 * rec.Cosomm
		x8 = a6 * rec.Cosomm

		rec.Z31 = 12.0*x1*x1 - 3.0*x3*x3
		rec.Z32 = 24.0*x1*x2 - 6.0*x3*x4
		rec.Z33 = 12.0*x2*x2 - 3.0*x4*x4
		rec.Z1 = 3.0*(a1*a1+a2*a2) + rec.Z31*rec.Emsq
		rec.Z2 = 6.0*(a1*a3+a2*a4) + rec.Z32*rec.Emsq
		rec.Z3 = 3.0*(a3*a3+a4*a4) + rec.Z33*rec.Emsq
		rec.Z11 = -6.0*a1*a5 + rec.Emsq*(-24.0*x1*x7-6.0*x3*x5)
		rec.Z12 = -6.0*(a1*a6+a3*a5) + rec.Emsq*
			(-24.0*(x2*x7+x1*x8)-6.0*(x3*x6+x4*x5))
		rec.Z13 = -6.0*a3*a6 + rec.Emsq*(-24.0*x2*x8-6.0*x4*x6)
		rec.Z21 = 6.0*a2*a5 + rec.Emsq*(24.0*x1*x5-6.0*x3*x7)
		rec.Z22 = 6.0*(a4*a5+a2*a6) + rec.Emsq*
			(24.0*(x2*x5+x1*x6)-6.0*(x4*x7+x3*x8))
		rec.Z23 = 6.0*a4*a6 + rec.Emsq*(24.0*x2*x6-6.0*x4*x8)
		rec.Z1 = rec.Z1 + rec.Z1 + betasq*rec.Z31
		rec.Z2 = rec.Z2 + rec.Z2 + betasq*rec.Z32
		rec.Z3 = rec.Z3 + rec.Z3 + betasq*rec.Z33
		rec.S3 = cc * xnoi
		rec.S2 = -0.5 * rec.S3 / rec.Rtemsq
		rec.S4 = rec.S3 * rec.Rtemsq
		rec.S1 = -15.0 * rec.Em * rec.S4
		rec.S5 = x1*x3 + x2*x4
		rec.S6 = x2*x3 + x1*x4
		rec.S7 = x2*x4 - x1*x3

		// Do lunar terms
		if lsflg == 1 {
			rec.Ss1 = rec.S1
			rec.Ss2 = rec.S2
			rec.Ss3 = rec.S3
			rec.Ss4 = rec.S4
			rec.Ss5 = rec.S5
			rec.Ss6 = rec.S6
			rec.Ss7 = rec.S7
			rec.Sz1 = rec.Z1
			rec.Sz2 = rec.Z2
			rec.Sz3 = rec.Z3
			rec.Sz11 = rec.Z11
			rec.Sz12 = rec.Z12
			rec.Sz13 = rec.Z13
			rec.Sz21 = rec.Z21
			rec.Sz22 = rec.Z22
			rec.Sz23 = rec.Z23
			rec.Sz31 = rec.Z31
			rec.Sz32 = rec.Z32
			rec.Sz33 = rec.Z33
			zcosg = zcosgl
			zsing = zsingl
			zcosi = zcosil
			zsini = zsinil
			zcosh = zcoshl*rec.Cnodm + zsinhl*rec.Snodm
			zsinh = rec.Snodm*zcoshl - rec.Cnodm*zsinhl
			cc = c1l
		}
	}

	rec.Zmol = math.Mod(4.7199672+0.22997150*rec.Day-rec.Gam, TWOPI_C)
	rec.Zmos = math.Mod(6.2565837+0.017201977*rec.Day, TWOPI_C)

	// Do solar terms
	rec.Se2 = 2.0 * rec.Ss1 * rec.Ss6
	rec.Se3 = 2.0 * rec.Ss1 * rec.Ss7
	rec.Si2 = 2.0 * rec.Ss2 * rec.Sz12
	rec.Si3 = 2.0 * rec.Ss2 * (rec.Sz13 - rec.Sz11)
	rec.Sl2 = -2.0 * rec.Ss3 * rec.Sz2
	rec.Sl3 = -2.0 * rec.Ss3 * (rec.Sz3 - rec.Sz1)
	rec.Sl4 = -2.0 * rec.Ss3 * (-21.0 - 9.0*rec.Emsq) * zes
	rec.Sgh2 = 2.0 * rec.Ss4 * rec.Sz32
	rec.Sgh3 = 2.0 * rec.Ss4 * (rec.Sz33 - rec.Sz31)
	rec.Sgh4 = -18.0 * rec.Ss4 * zes
	rec.Sh2 = -2.0 * rec.Ss2 * rec.Sz22
	rec.Sh3 = -2.0 * rec.Ss2 * (rec.Sz23 - rec.Sz21)

	// Do lunar terms
	rec.Ee2 = 2.0 * rec.S1 * rec.S6
	rec.E3 = 2.0 * rec.S1 * rec.S7
	rec.Xi2 = 2.0 * rec.S2 * rec.Z12
	rec.Xi3 = 2.0 * rec.S2 * (rec.Z13 - rec.Z11)
	rec.Xl2 = -2.0 * rec.S3 * rec.Z2
	rec.Xl3 = -2.0 * rec.S3 * (rec.Z3 - rec.Z1)
	rec.Xl4 = -2.0 * rec.S3 * (-21.0 - 9.0*rec.Emsq) * zel
	rec.Xgh2 = 2.0 * rec.S4 * rec.Z32
	rec.Xgh3 = 2.0 * rec.S4 * (rec.Z33 - rec.Z31)
	rec.Xgh4 = -18.0 * rec.S4 * zel
	rec.Xh2 = -2.0 * rec.S2 * rec.Z22
	rec.Xh3 = -2.0 * rec.S2 * (rec.Z23 - rec.Z21)
}

// Dpper provides deep space long period periodic contributions - exact replica from C code
func Dpper(e3, ee2, peo, pgho, pho, pinco, plo, se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4, t, xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4, zmol, zmos float64, init byte, rec *ElsetRec, opsmode byte) {
	// Local variables
	var alfdp, betdp, cosip, cosop, dalf, dbet, dls, f2, f3, pe, pgh, ph, pinc, pl float64
	var sel, ses, sghl, sghs, shll, shs, sil, sinip, sinop, sinzf, sis, sll, sls, xls, xnoh, zf, zm, zel, zes, znl, zns float64

	// Constants
	zns = 1.19459e-5
	zes = 0.01675
	znl = 1.5835218e-4
	zel = 0.05490

	// Calculate time varying periodics
	zm = zmos + zns*t
	// Be sure that the initial call has time set to zero
	if init == 'y' {
		zm = zmos
	}
	zf = zm + 2.0*zes*math.Sin(zm)
	sinzf = math.Sin(zf)
	f2 = 0.5*sinzf*sinzf - 0.25
	f3 = -0.5 * sinzf * math.Cos(zf)
	ses = se2*f2 + se3*f3
	sis = si2*f2 + si3*f3
	sls = sl2*f2 + sl3*f3 + sl4*sinzf
	sghs = sgh2*f2 + sgh3*f3 + sgh4*sinzf
	shs = sh2*f2 + sh3*f3
	zm = zmol + znl*t
	if init == 'y' {
		zm = zmol
	}
	zf = zm + 2.0*zel*math.Sin(zm)
	sinzf = math.Sin(zf)
	f2 = 0.5*sinzf*sinzf - 0.25
	f3 = -0.5 * sinzf * math.Cos(zf)
	sel = ee2*f2 + e3*f3
	sil = xi2*f2 + xi3*f3
	sll = xl2*f2 + xl3*f3 + xl4*sinzf
	sghl = xgh2*f2 + xgh3*f3 + xgh4*sinzf
	shll = xh2*f2 + xh3*f3
	pe = ses + sel
	pinc = sis + sil
	pl = sls + sll
	pgh = sghs + sghl
	ph = shs + shll

	if init == 'n' {
		pe = pe - peo
		pinc = pinc - pinco
		pl = pl - plo
		pgh = pgh - pgho
		ph = ph - pho
		rec.Inclp = rec.Inclp + pinc
		rec.Ep = rec.Ep + pe
		sinip = math.Sin(rec.Inclp)
		cosip = math.Cos(rec.Inclp)

		// Apply periodics directly
		// sgp4fix for lyddane choice
		// use next line for gsfc version and perturbed inclination
		if rec.Inclp >= 0.2 {
			ph = ph / sinip
			pgh = pgh - cosip*ph
			rec.Argpp = rec.Argpp + pgh
			rec.Nodep = rec.Nodep + ph
			rec.Mp = rec.Mp + pl
		} else {
			// Apply periodics with lyddane modification
			sinop = math.Sin(rec.Nodep)
			cosop = math.Cos(rec.Nodep)
			alfdp = sinip * sinop
			betdp = sinip * cosop
			dalf = ph*cosop + pinc*cosip*sinop
			dbet = -ph*sinop + pinc*cosip*cosop
			alfdp = alfdp + dalf
			betdp = betdp + dbet
			rec.Nodep = math.Mod(rec.Nodep, TWOPI_C)
			// sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if (rec.Nodep < 0.0) && (opsmode == 'a') {
				rec.Nodep = rec.Nodep + TWOPI_C
			}
			xls = rec.Mp + rec.Argpp + cosip*rec.Nodep
			dls = pl + pgh - pinc*rec.Nodep*sinip
			xls = xls + dls
			xls = math.Mod(xls, TWOPI_C)
			xnoh = rec.Nodep
			rec.Nodep = math.Atan2(alfdp, betdp)
			// sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if (rec.Nodep < 0.0) && (opsmode == 'a') {
				rec.Nodep = rec.Nodep + TWOPI_C
			}
			if math.Abs(xnoh-rec.Nodep) > PI_C {
				if rec.Nodep < xnoh {
					rec.Nodep = rec.Nodep + TWOPI_C
				} else {
					rec.Nodep = rec.Nodep - TWOPI_C
				}
			}
			rec.Mp = rec.Mp + pl
			rec.Argpp = xls - rec.Mp - cosip*rec.Nodep
		}
	}
}

// Dsinit provides deep space contributions to mean motion dot - exact replica from C code
func Dsinit(tc, xpidot float64, rec *ElsetRec) {
	// Local variables
	var ainv2, aonv, cosisq, eoc, f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543 float64
	var g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533 float64
	var ses, sgs, sghl, sghs, shs, shll, sis, sini2, sls, temp, temp1, theta, xno2 float64
	var q22, q31, q33, root22, root44, root54, rptim, root32, root52, x2o3, znl, emo, zns, emsqo float64

	q22 = 1.7891679e-6
	q31 = 2.1460748e-6
	q33 = 2.2123015e-7
	root22 = 1.7891679e-6
	root44 = 7.3636953e-9
	root54 = 2.1765803e-9
	rptim = 4.37526908801129966e-3 // this equates to 7.29211514668855e-5 rad/sec
	root32 = 3.7393792e-7
	root52 = 1.1428639e-7
	x2o3 = 2.0 / 3.0
	znl = 1.5835218e-4
	zns = 1.19459e-5

	// Deep space initialization
	rec.Irez = 0
	if (rec.Nm < 0.0052359877) && (rec.Nm > 0.0034906585) {
		rec.Irez = 1
	}
	if (rec.Nm >= 8.26e-3) && (rec.Nm <= 9.24e-3) && (rec.Em >= 0.5) {
		rec.Irez = 2
	}

	// Do solar terms
	ses = rec.Ss1 * zns * rec.Ss5
	sis = rec.Ss2 * zns * (rec.Sz11 + rec.Sz13)
	sls = -zns * rec.Ss3 * (rec.Sz1 + rec.Sz3 - 14.0 - 6.0*rec.Emsq)
	sghs = rec.Ss4 * zns * (rec.Sz31 + rec.Sz33 - 6.0)
	shs = -zns * rec.Ss2 * (rec.Sz21 + rec.Sz23)
	// sgp4fix for 180 deg incl
	if (rec.Inclm < 5.2359877e-2) || (rec.Inclm > PI_C-5.2359877e-2) {
		shs = 0.0
	}
	if rec.Sinim != 0.0 {
		shs = shs / rec.Sinim
	}
	sgs = sghs - rec.Cosim*shs

	// Do lunar terms
	rec.Dedt = ses + rec.S1*znl*rec.S5
	rec.Didt = sis + rec.S2*znl*(rec.Z11+rec.Z13)
	rec.Dmdt = sls - znl*rec.S3*(rec.Z1+rec.Z3-14.0-6.0*rec.Emsq)
	sghl = rec.S4 * znl * (rec.Z31 + rec.Z33 - 6.0)
	shll = -znl * rec.S2 * (rec.Z21 + rec.Z23)
	// sgp4fix for 180 deg incl
	if (rec.Inclm < 5.2359877e-2) || (rec.Inclm > PI_C-5.2359877e-2) {
		shll = 0.0
	}
	rec.Domdt = sgs + sghl
	rec.Dnodt = shs
	if rec.Sinim != 0.0 {
		rec.Domdt = rec.Domdt - rec.Cosim/rec.Sinim*shll
		rec.Dnodt = rec.Dnodt + shll/rec.Sinim
	}

	// Calculate deep space resonance effects
	rec.Dndt = 0.0
	theta = math.Mod(rec.Gsto+tc*rptim, TWOPI_C)
	rec.Em = rec.Em + rec.Dedt*rec.T
	rec.Inclm = rec.Inclm + rec.Didt*rec.T
	rec.Argpm = rec.Argpm + rec.Domdt*rec.T
	rec.Nodem = rec.Nodem + rec.Dnodt*rec.T
	rec.Mm = rec.Mm + rec.Dmdt*rec.T

	// Initialize the resonance terms
	if rec.Irez != 0 {
		aonv = math.Pow(rec.Nm/rec.Xke, x2o3)

		// Geopotential resonance for 12 hour orbits
		if rec.Irez == 2 {
			cosisq = rec.Cosim * rec.Cosim
			emo = rec.Em
			rec.Em = rec.Ecco
			emsqo = rec.Emsq
			rec.Emsq = rec.Eccsq
			eoc = rec.Em * rec.Emsq
			g201 = -0.306 - (rec.Em-0.64)*0.440

			if rec.Em <= 0.65 {
				g211 = 3.616 - 13.2470*rec.Em + 16.2900*rec.Emsq
				g310 = -19.302 + 117.3900*rec.Em - 228.4190*rec.Emsq + 156.5910*eoc
				g322 = -18.9068 + 109.7927*rec.Em - 214.6334*rec.Emsq + 146.5816*eoc
				g410 = -41.122 + 242.6940*rec.Em - 471.0940*rec.Emsq + 313.9530*eoc
				g422 = -146.407 + 841.8800*rec.Em - 1629.014*rec.Emsq + 1083.4350*eoc
				g520 = -532.114 + 3017.977*rec.Em - 5740.032*rec.Emsq + 3708.2760*eoc
			} else {
				g211 = -72.099 + 331.819*rec.Em - 508.738*rec.Emsq + 266.724*eoc
				g310 = -346.844 + 1582.851*rec.Em - 2415.925*rec.Emsq + 1246.113*eoc
				g322 = -342.585 + 1554.908*rec.Em - 2366.899*rec.Emsq + 1215.972*eoc
				g410 = -1052.797 + 4758.686*rec.Em - 7193.992*rec.Emsq + 3651.957*eoc
				g422 = -3581.690 + 16178.110*rec.Em - 24462.770*rec.Emsq + 12422.520*eoc
				if rec.Em > 0.715 {
					g520 = -5149.66 + 29936.92*rec.Em - 54087.36*rec.Emsq + 31324.56*eoc
				} else {
					g520 = 1464.74 - 4664.75*rec.Em + 3763.64*rec.Emsq
				}
			}
			if rec.Em < 0.7 {
				g533 = -919.22770 + 4988.6100*rec.Em - 9064.7700*rec.Emsq + 5542.21*eoc
				g521 = -822.71072 + 4568.6173*rec.Em - 8491.4146*rec.Emsq + 5337.524*eoc
				g532 = -853.66600 + 4690.2500*rec.Em - 8624.7700*rec.Emsq + 5341.4*eoc
			} else {
				g533 = -37995.780 + 161616.52*rec.Em - 229838.20*rec.Emsq + 109377.94*eoc
				g521 = -51752.104 + 218913.95*rec.Em - 309468.16*rec.Emsq + 146349.42*eoc
				g532 = -40023.880 + 170470.89*rec.Em - 242699.48*rec.Emsq + 115605.82*eoc
			}

			sini2 = rec.Sinim * rec.Sinim
			f220 = 0.75 * (1.0 + 2.0*rec.Cosim + cosisq)
			f221 = 1.5 * sini2
			f321 = 1.875 * rec.Sinim * (1.0 - 2.0*rec.Cosim - 3.0*cosisq)
			f322 = -1.875 * rec.Sinim * (1.0 + 2.0*rec.Cosim - 3.0*cosisq)
			f441 = 35.0 * sini2 * f220
			f442 = 39.3750 * sini2 * sini2
			f522 = 9.84375 * rec.Sinim * (sini2*(1.0-2.0*rec.Cosim-5.0*cosisq) +
				0.33333333*(-2.0+4.0*rec.Cosim+6.0*cosisq))
			f523 = rec.Sinim * (4.92187512*sini2*(-2.0-4.0*rec.Cosim+
				10.0*cosisq) + 6.56250012*(1.0+2.0*rec.Cosim-3.0*cosisq))
			f542 = 29.53125 * rec.Sinim * (2.0 - 8.0*rec.Cosim + cosisq*
				(-12.0+8.0*rec.Cosim+10.0*cosisq))
			f543 = 29.53125 * rec.Sinim * (-2.0 - 8.0*rec.Cosim + cosisq*
				(12.0+8.0*rec.Cosim-10.0*cosisq))
			xno2 = rec.Nm * rec.Nm
			ainv2 = aonv * aonv
			temp1 = 3.0 * xno2 * ainv2
			temp = temp1 * root22
			rec.D2201 = temp * f220 * g201
			rec.D2211 = temp * f221 * g211
			temp1 = temp1 * aonv
			temp = temp1 * root32
			rec.D3210 = temp * f321 * g310
			rec.D3222 = temp * f322 * g322
			temp1 = temp1 * aonv
			temp = 2.0 * temp1 * root44
			rec.D4410 = temp * f441 * g410
			rec.D4422 = temp * f442 * g422
			temp1 = temp1 * aonv
			temp = temp1 * root52
			rec.D5220 = temp * f522 * g520
			rec.D5232 = temp * f523 * g532
			temp = 2.0 * temp1 * root54
			rec.D5421 = temp * f542 * g521
			rec.D5433 = temp * f543 * g533
			rec.Xlamo = math.Mod(rec.Mo+rec.Nodeo+rec.Nodeo-theta-theta, TWOPI_C)
			rec.Xfact = rec.Mdot + rec.Dmdt + 2.0*(rec.Nodedot+rec.Dnodt-rptim) - rec.NoUnkozai
			rec.Em = emo
			rec.Emsq = emsqo
		}

		// Synchronous resonance terms
		if rec.Irez == 1 {
			g200 = 1.0 + rec.Emsq*(-2.5+0.8125*rec.Emsq)
			g310 = 1.0 + 2.0*rec.Emsq
			g300 = 1.0 + rec.Emsq*(-6.0+6.60937*rec.Emsq)
			f220 = 0.75 * (1.0 + rec.Cosim) * (1.0 + rec.Cosim)
			f311 = 0.9375*rec.Sinim*rec.Sinim*(1.0+3.0*rec.Cosim) - 0.75*(1.0+rec.Cosim)
			f330 = 1.0 + rec.Cosim
			f330 = 1.875 * f330 * f330 * f330
			rec.Del1 = 3.0 * rec.Nm * rec.Nm * aonv * aonv
			rec.Del2 = 2.0 * rec.Del1 * f220 * g200 * q22
			rec.Del3 = 3.0 * rec.Del1 * f330 * g300 * q33 * aonv
			rec.Del1 = rec.Del1 * f311 * g310 * q31 * aonv
			rec.Xlamo = math.Mod(rec.Mo+rec.Nodeo+rec.Argpo-theta, TWOPI_C)
			rec.Xfact = rec.Mdot + xpidot - rptim + rec.Dmdt + rec.Domdt + rec.Dnodt - rec.NoUnkozai
		}

		// For sgp4, initialize the integrator
		rec.Xli = rec.Xlamo
		rec.Xni = rec.NoUnkozai
		rec.Atime = 0.0
		rec.Nm = rec.NoUnkozai + rec.Dndt
	}
}

// Dspace provides deep space contributions to mean elements - exact replica from C code
func Dspace(tc float64, rec *ElsetRec) {
	var iretn int
	var delt, ft, theta, x2li, x2omi, xl, xldot, xnddt, xndt, xomi, g22, g32, g44, g52, g54, fasx2, fasx4, fasx6, rptim, step2, stepn, stepp float64

	xndt = 0
	xnddt = 0
	xldot = 0

	fasx2 = 0.13130908
	fasx4 = 2.8843198
	fasx6 = 0.37448087
	g22 = 5.7686396
	g32 = 0.95240898
	g44 = 1.8014998
	g52 = 1.0508330
	g54 = 4.4108898
	rptim = 4.37526908801129966e-3 // this equates to 7.29211514668855e-5 rad/sec
	stepp = 720.0
	stepn = -720.0
	step2 = 259200.0

	// Calculate deep space resonance effects
	rec.Dndt = 0.0
	theta = math.Mod(rec.Gsto+tc*rptim, TWOPI_C)
	rec.Em = rec.Em + rec.Dedt*rec.T

	rec.Inclm = rec.Inclm + rec.Didt*rec.T
	rec.Argpm = rec.Argpm + rec.Domdt*rec.T
	rec.Nodem = rec.Nodem + rec.Dnodt*rec.T
	rec.Mm = rec.Mm + rec.Dmdt*rec.T

	// Update resonances : numerical (euler-maclaurin) integration
	// Epoch restart
	ft = 0.0
	if rec.Irez != 0 {
		// sgp4fix streamline check
		if (rec.Atime == 0.0) || (rec.T*rec.Atime <= 0.0) || (math.Abs(rec.T) < math.Abs(rec.Atime)) {
			rec.Atime = 0.0
			rec.Xni = rec.NoUnkozai
			rec.Xli = rec.Xlamo
		}
		// sgp4fix move check outside loop
		if rec.T > 0.0 {
			delt = stepp
		} else {
			delt = stepn
		}

		iretn = 381 // added for do loop
		for iretn == 381 {
			// Dot terms calculated
			// Near - synchronous resonance terms
			if rec.Irez != 2 {
				xndt = rec.Del1*math.Sin(rec.Xli-fasx2) + rec.Del2*math.Sin(2.0*(rec.Xli-fasx4)) +
					rec.Del3*math.Sin(3.0*(rec.Xli-fasx6))
				xldot = rec.Xni + rec.Xfact
				xnddt = rec.Del1*math.Cos(rec.Xli-fasx2) +
					2.0*rec.Del2*math.Cos(2.0*(rec.Xli-fasx4)) +
					3.0*rec.Del3*math.Cos(3.0*(rec.Xli-fasx6))
				xnddt = xnddt * xldot
			} else {
				// Near - half-day resonance terms
				xomi = rec.Argpo + rec.Argpdot*rec.Atime
				x2omi = xomi + xomi
				x2li = rec.Xli + rec.Xli
				xndt = rec.D2201*math.Sin(x2omi+rec.Xli-g22) + rec.D2211*math.Sin(rec.Xli-g22) +
					rec.D3210*math.Sin(xomi+rec.Xli-g32) + rec.D3222*math.Sin(-xomi+rec.Xli-g32) +
					rec.D4410*math.Sin(x2omi+x2li-g44) + rec.D4422*math.Sin(x2li-g44) +
					rec.D5220*math.Sin(xomi+rec.Xli-g52) + rec.D5232*math.Sin(-xomi+rec.Xli-g52) +
					rec.D5421*math.Sin(xomi+x2li-g54) + rec.D5433*math.Sin(-xomi+x2li-g54)
				xldot = rec.Xni + rec.Xfact
				xnddt = rec.D2201*math.Cos(x2omi+rec.Xli-g22) + rec.D2211*math.Cos(rec.Xli-g22) +
					rec.D3210*math.Cos(xomi+rec.Xli-g32) + rec.D3222*math.Cos(-xomi+rec.Xli-g32) +
					rec.D5220*math.Cos(xomi+rec.Xli-g52) + rec.D5232*math.Cos(-xomi+rec.Xli-g52) +
					2.0*(rec.D4410*math.Cos(x2omi+x2li-g44)+
						rec.D4422*math.Cos(x2li-g44)+rec.D5421*math.Cos(xomi+x2li-g54)+
						rec.D5433*math.Cos(-xomi+x2li-g54))
				xnddt = xnddt * xldot
			}

			// Integrator
			// sgp4fix move end checks to end of routine
			if math.Abs(rec.T-rec.Atime) >= stepp {
				iretn = 381
			} else { // exit here
				ft = rec.T - rec.Atime
				iretn = 0
			}

			if iretn == 381 {
				rec.Xli = rec.Xli + xldot*delt + xndt*step2
				rec.Xni = rec.Xni + xndt*delt + xnddt*step2
				rec.Atime = rec.Atime + delt
			}
		}

		rec.Nm = rec.Xni + xndt*ft + xnddt*ft*ft*0.5
		xl = rec.Xli + xldot*ft + xndt*ft*ft*0.5
		if rec.Irez != 1 {
			rec.Mm = xl - 2.0*rec.Nodem + 2.0*theta
			rec.Dndt = rec.Nm - rec.NoUnkozai
		} else {
			rec.Mm = xl - rec.Nodem - rec.Argpm + theta
			rec.Dndt = rec.Nm - rec.NoUnkozai
		}
		rec.Nm = rec.NoUnkozai + rec.Dndt
	}
}

// SGP4 is the main SGP4 prediction model - exact replica from C code
// This is the combined version that handles both near-earth and deep-space cases
func SGP4(satrec *ElsetRec, tsince float64, r, v []float64) bool {
	var axnl, aynl, betal, cnod, cos2u, coseo1, cosi, cosip, cosisq, cossu, cosu float64
	var delm, delomg, ecose, el2, eo1, esine, argpdf, pl, mrt, mvt, rdotl, rl, rvdot, rvdotl float64
	var sin2u, sineo1, sini, sinip, sinsu, sinu, snod, su, t2, t3, t4, tem5, temp, temp1, temp2, tempa, tempe, templ, u, ux, uy, uz, vx, vy, vz float64
	var xinc, xincp, xl, xlm, xmdf, xmx, xmy, nodedf, xnode, tc, x2o3, vkmpersec, delmtemp float64
	var ktr int

	// Set mathematical constants
	const temp4 = 1.5e-12
	x2o3 = 2.0 / 3.0
	vkmpersec = satrec.Radiusearthkm * satrec.Xke / 60.0

	// Clear sgp4 error flag
	satrec.T = tsince
	satrec.Error = 0

	// Update for secular gravity and atmospheric drag
	xmdf = satrec.Mo + satrec.Mdot*satrec.T
	argpdf = satrec.Argpo + satrec.Argpdot*satrec.T
	nodedf = satrec.Nodeo + satrec.Nodedot*satrec.T
	satrec.Argpm = argpdf
	satrec.Mm = xmdf
	t2 = satrec.T * satrec.T
	satrec.Nodem = nodedf + satrec.Nodecf*t2
	tempa = 1.0 - satrec.Cc1*satrec.T
	tempe = satrec.Bstar * satrec.Cc4 * satrec.T
	templ = satrec.T2cof * t2

	delomg = 0
	delmtemp = 0
	delm = 0
	_ = 0 // temp
	t3 = 0
	t4 = 0
	mrt = 0

	if satrec.Isimp != 1 {
		delomg = satrec.Omgcof * satrec.T
		// sgp4fix use multiply for speed instead of pow
		delmtemp = 1.0 + satrec.Eta*math.Cos(xmdf)
		delm = satrec.Xmcof *
			(delmtemp*delmtemp*delmtemp -
				satrec.Delmo)
		temp = delomg + delm
		satrec.Mm = xmdf + temp
		satrec.Argpm = argpdf - temp
		t3 = t2 * satrec.T
		t4 = t3 * satrec.T
		tempa = tempa - satrec.D2*t2 - satrec.D3*t3 -
			satrec.D4*t4
		tempe = tempe + satrec.Bstar*satrec.Cc5*(math.Sin(satrec.Mm)-satrec.Sinmao)
		templ = templ + satrec.T3cof*t3 + t4*(satrec.T4cof+satrec.T*satrec.T5cof)
	}

	tc = 0
	satrec.Nm = satrec.NoUnkozai
	satrec.Em = satrec.Ecco
	satrec.Inclm = satrec.Inclo
	if satrec.Method == 'd' {
		tc = satrec.T
		Dspace(tc, satrec)
	} // if method = d

	if satrec.Nm <= 0.0 {
		satrec.Error = 2
		return false
	}

	satrec.Am = math.Pow((satrec.Xke/satrec.Nm), x2o3) * tempa * tempa
	satrec.Nm = satrec.Xke / math.Pow(satrec.Am, 1.5)
	satrec.Em = satrec.Em - tempe

	// Fix tolerance for error recognition
	// sgp4fix am is fixed from the previous nm check
	if (satrec.Em >= 1.0) || (satrec.Em < -0.001) {
		satrec.Error = 1
		return false
	}
	// sgp4fix fix tolerance to avoid a divide by zero
	if satrec.Em < 1.0e-6 {
		satrec.Em = 1.0e-6
	}
	satrec.Mm = satrec.Mm + satrec.NoUnkozai*templ
	xlm = satrec.Mm + satrec.Argpm + satrec.Nodem
	satrec.Emsq = satrec.Em * satrec.Em
	_ = 1.0 - satrec.Emsq // temp

	satrec.Nodem = math.Mod(satrec.Nodem, TWOPI_C)
	satrec.Argpm = math.Mod(satrec.Argpm, TWOPI_C)
	xlm = math.Mod(xlm, TWOPI_C)
	satrec.Mm = math.Mod(xlm-satrec.Argpm-satrec.Nodem, TWOPI_C)

	// sgp4fix recover singly averaged mean elements
	satrec.Am = satrec.Am //nolint:ineffassign
	satrec.Em = satrec.Em //nolint:ineffassign
	satrec.Im = satrec.Inclm
	satrec.Om = satrec.Nodem
	satrec.Om2 = satrec.Argpm
	satrec.Mm = satrec.Mm //nolint:ineffassign
	satrec.Nm = satrec.Nm //nolint:ineffassign

	// Compute extra mean quantities
	satrec.Sinim = math.Sin(satrec.Inclm)
	satrec.Cosim = math.Cos(satrec.Inclm)

	// Add lunar-solar periodics
	satrec.Ep = satrec.Em
	xincp = satrec.Inclm
	satrec.Inclp = satrec.Inclm
	satrec.Argpp = satrec.Argpm
	satrec.Nodep = satrec.Nodem
	satrec.Mp = satrec.Mm
	sinip = satrec.Sinim
	cosip = satrec.Cosim
	if satrec.Method == 'd' {
		Dpper(satrec.E3, satrec.Ee2, satrec.Peo, satrec.Pgho,
			satrec.Pho, satrec.Pinco, satrec.Plo, satrec.Se2,
			satrec.Se3, satrec.Sgh2, satrec.Sgh3, satrec.Sgh4,
			satrec.Sh2, satrec.Sh3, satrec.Si2, satrec.Si3,
			satrec.Sl2, satrec.Sl3, satrec.Sl4, satrec.T,
			satrec.Xgh2, satrec.Xgh3, satrec.Xgh4, satrec.Xh2,
			satrec.Xh3, satrec.Xi2, satrec.Xi3, satrec.Xl2,
			satrec.Xl3, satrec.Xl4, satrec.Zmol, satrec.Zmos,
			'n', satrec, satrec.OperationMode)

		xincp = satrec.Inclp
		if xincp < 0.0 {
			xincp = -xincp
			satrec.Nodep = satrec.Nodep + PI_C
			satrec.Argpp = satrec.Argpp - PI_C
		}
		if (satrec.Ep < 0.0) || (satrec.Ep > 1.0) {
			satrec.Error = 3
			return false
		}
	} // if method = d

	// Long period periodics
	if satrec.Method == 'd' {
		sinip = math.Sin(xincp)
		cosip = math.Cos(xincp)
		satrec.Aycof = -0.5 * satrec.J3oj2 * sinip
		// sgp4fix for divide by zero for xincp = 180 deg
		if math.Abs(cosip+1.0) > 1.5e-12 {
			satrec.Xlcof = -0.25 * satrec.J3oj2 * sinip * (3.0 + 5.0*cosip) / (1.0 + cosip)
		} else {
			satrec.Xlcof = -0.25 * satrec.J3oj2 * sinip * (3.0 + 5.0*cosip) / temp4
		}
	}
	axnl = satrec.Ep * math.Cos(satrec.Argpp)
	temp = 1.0 / (satrec.Am * (1.0 - satrec.Ep*satrec.Ep))
	aynl = satrec.Ep*math.Sin(satrec.Argpp) + temp*satrec.Aycof
	xl = satrec.Mp + satrec.Argpp + satrec.Nodep + temp*satrec.Xlcof*axnl

	// Solve kepler's equation
	u = math.Mod(xl-satrec.Nodep, TWOPI_C)
	eo1 = u
	tem5 = 9999.9
	ktr = 1
	sineo1 = 0
	coseo1 = 0
	// sgp4fix for kepler iteration
	for (math.Abs(tem5) >= 1.0e-12) && (ktr <= 10) {
		sineo1 = math.Sin(eo1)
		coseo1 = math.Cos(eo1)
		tem5 = 1.0 - coseo1*axnl - sineo1*aynl
		tem5 = (u - aynl*coseo1 + axnl*sineo1 - eo1) / tem5
		if math.Abs(tem5) >= 0.95 {
			if tem5 > 0.0 {
				tem5 = 0.95
			} else {
				tem5 = -0.95
			}
		}
		eo1 = eo1 + tem5
		ktr = ktr + 1
	}

	// Short period preliminary quantities
	ecose = axnl*coseo1 + aynl*sineo1
	esine = axnl*sineo1 - aynl*coseo1
	el2 = axnl*axnl + aynl*aynl
	pl = satrec.Am * (1.0 - el2)
	if pl < 0.0 {
		satrec.Error = 4
		return false
	} else {
		rl = satrec.Am * (1.0 - ecose)
		rdotl = math.Sqrt(satrec.Am) * esine / rl
		rvdotl = math.Sqrt(pl) / rl
		betal = math.Sqrt(1.0 - el2)
		temp = esine / (1.0 + betal)
		sinu = satrec.Am / rl * (sineo1 - aynl - axnl*temp)
		cosu = satrec.Am / rl * (coseo1 - axnl + aynl*temp)
		su = math.Atan2(sinu, cosu)
		sin2u = (cosu + cosu) * sinu
		cos2u = 1.0 - 2.0*sinu*sinu
		temp = 1.0 / pl
		temp1 = 0.5 * satrec.J2 * temp
		temp2 = temp1 * temp

		// Update for short period periodics
		if satrec.Method == 'd' {
			cosisq = cosip * cosip
			satrec.Con41 = 3.0*cosisq - 1.0
			satrec.X1mth2 = 1.0 - cosisq
			satrec.X7thm1 = 7.0*cosisq - 1.0
		}
		mrt = rl*(1.0-1.5*temp2*betal*satrec.Con41) + 0.5*temp1*satrec.X1mth2*cos2u
		su = su - 0.25*temp2*satrec.X7thm1*sin2u
		xnode = satrec.Nodep + 1.5*temp2*cosip*sin2u
		xinc = xincp + 1.5*temp2*cosip*sinip*cos2u
		mvt = rdotl - satrec.Nm*temp1*satrec.X1mth2*sin2u/satrec.Xke
		rvdot = rvdotl + satrec.Nm*temp1*(satrec.X1mth2*cos2u+1.5*satrec.Con41)/satrec.Xke

		// Orientation vectors
		sinsu = math.Sin(su)
		cossu = math.Cos(su)
		snod = math.Sin(xnode)
		cnod = math.Cos(xnode)
		sini = math.Sin(xinc)
		cosi = math.Cos(xinc)
		xmx = -snod * cosi
		xmy = cnod * cosi
		ux = xmx*sinsu + cnod*cossu
		uy = xmy*sinsu + snod*cossu
		uz = sini * sinsu
		vx = xmx*cossu - cnod*sinsu
		vy = xmy*cossu - snod*sinsu
		vz = sini * cossu

		// Position and velocity (in km and km/sec)
		r[0] = (mrt * ux) * satrec.Radiusearthkm
		r[1] = (mrt * uy) * satrec.Radiusearthkm
		r[2] = (mrt * uz) * satrec.Radiusearthkm
		v[0] = (mvt*ux + rvdot*vx) * vkmpersec
		v[1] = (mvt*uy + rvdot*vy) * vkmpersec
		v[2] = (mvt*uz + rvdot*vz) * vkmpersec
	} // if pl > 0

	// sgp4fix for decaying satellites
	if mrt < 1.0 {
		satrec.Error = 6
		return false
	}

	return true
}

// PropagateSatellite is a wrapper function that uses the C implementation
// This is the combined SGP4/SDP4 model that handles both near-earth and deep-space cases
func PropagateSatellite(tle *TLE, time float64, units ...string) (Vector, Vector, error) {
	// Convert TLE to ElsetRec format - using exact C conversion factors
	xpdotp := 1440.0 / (2.0 * PI_C) // 229.1831180523293 - same as C code

	// Apply BStar exponent like C code does
	bstarWithExp := tle.BStar * math.Pow(10.0, float64(tle.BStarExp))

	satrec := &ElsetRec{
		WhichConst:  WGS72,
		Jdsatepoch:  tle.JulianEpoch,
		JdsatepochF: 0.0,
		Bstar:       bstarWithExp,
		Inclo:       tle.Inclination * DEG2RAD,
		Nodeo:       tle.RightAscension * DEG2RAD,
		Ecco:        tle.Eccentricity,
		Argpo:       tle.ArgumentPerigee * DEG2RAD,
		Mo:          tle.MeanAnomaly * DEG2RAD,
		NoKozai:     tle.MeanMotion / xpdotp,                         // Same as C: tle->n/xpdotp
		Ndot:        tle.MeanMotionDot / (xpdotp * 1440.0),           // Same as C: tle->ndot / (xpdotp*1440.0)
		Nddot:       tle.MeanMotionDDot / (xpdotp * 1440.0 * 1440.0), // Same as C: tle->nddot / (xpdotp*1440.0*1440.0)
		EpochYr:     int(tle.Epoch / 1000),
		EpochDays:   math.Mod(tle.Epoch, 1000),
	}

	// Initialize the satellite record
	if !SGP4Init('a', satrec) {
		return Vector{}, Vector{}, fmt.Errorf("SGP4 initialization failed")
	}

	// Calculate time since epoch in minutes
	timeSinceEpoch := (time - tle.JulianEpoch) * 1440.0

	// Propagate
	var r, v [3]float64
	if !SGP4(satrec, timeSinceEpoch, r[:], v[:]) {
		return Vector{}, Vector{}, fmt.Errorf("SGP4 propagation failed with error code %d", satrec.Error)
	}

	// Convert to Vector format
	pos := Vector{X: r[0], Y: r[1], Z: r[2]}
	vel := Vector{X: v[0], Y: v[1], Z: v[2]}

	// Convert units if requested
	if len(units) > 0 && units[0] == "km" {
		// Already in km from the C implementation
		return pos, vel, nil
	}

	// Default: convert to Earth radii
	pos.X /= 6378.135
	pos.Y /= 6378.135
	pos.Z /= 6378.135
	vel.X /= 6378.135
	vel.Y /= 6378.135
	vel.Z /= 6378.135

	return pos, vel, nil
}
