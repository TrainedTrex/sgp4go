package sgp4go

import (
	"math"
	"time"
)

// JulianDateOfEpoch converts a TLE epoch (year + day of year) to Julian date
func JulianDateOfEpoch(epoch float64) float64 {
	year := int(epoch / 1000.0)
	dayOfYear := epoch - float64(year*1000)

	// Y2K modification from original Pascal code
	// Valid 1957 through 2056
	if year < 57 {
		year += 2000
	} else if year < 100 {
		year += 1900
	}

	// Calculate Julian date for January 1st of the year
	jdJan1 := JulianDate(year, 1, 1)

	// Add the day of year (subtract 1 since day 1 is January 1st)
	return jdJan1 + dayOfYear - 1.0
}

// JulianDate calculates the Julian date for a given year, month, and day
// Based on Astronomical Formulae for Calculators, Jean Meeus, pages 23-25
func JulianDate(year, month, day int) float64 {
	if month <= 2 {
		year -= 1
		month += 12
	}

	a := math.Trunc(float64(year) / 100.0)
	b := 2.0 - a + math.Trunc(a/4.0)

	jd := math.Trunc(365.25*(float64(year)+4716.0)) +
		math.Trunc(30.6001*float64(month+1)) +
		float64(day) + b - 1524.5

	return jd
}

// JulianDateToTime converts a Julian date to a Go time.Time
// Based on Astronomical Formulae for Calculators, Jean Meeus, pages 26-27
func JulianDateToTime(jd float64) time.Time {
	// Convert Julian date to Gregorian date
	jd += 0.5
	z := math.Trunc(jd)
	f := jd - z

	var a float64
	if z < 2299161.0 {
		a = z
	} else {
		alpha := math.Trunc((z - 1867216.25) / 36524.25)
		a = z + 1.0 + alpha - math.Trunc(alpha/4.0)
	}

	b := a + 1524.0
	c := math.Trunc((b - 122.1) / 365.25)
	d := math.Trunc(365.25 * c)
	e := math.Trunc((b - d) / 30.6001)

	day := b - d - math.Trunc(30.6001*e) + f
	month := int(e - 1.0)
	if month > 12 {
		month -= 12
	}
	year := int(c - 4715.0)
	if month > 2 {
		year -= 1
	}

	// Convert to time.Time (assumes UTC)
	return time.Date(year, time.Month(month), int(day), 0, 0, 0, 0, time.UTC)
}

// TimeToJulianDate converts a Go time.Time to Julian date
func TimeToJulianDate(t time.Time) float64 {
	year := t.Year()
	month := int(t.Month())
	day := t.Day()

	return JulianDate(year, month, day) + float64(t.Hour())/24.0 +
		float64(t.Minute())/1440.0 + float64(t.Second())/86400.0
}

// DaysSinceEpoch calculates the number of days since the TLE epoch
func DaysSinceEpoch(julianEpoch, julianDate float64) float64 {
	return julianDate - julianEpoch
}

// MinutesSinceEpoch calculates the number of minutes since the TLE epoch
func MinutesSinceEpoch(julianEpoch, julianDate float64) float64 {
	return (julianDate - julianEpoch) * 1440.0 // 1440 minutes per day
}
