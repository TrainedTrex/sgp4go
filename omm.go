package sgp4go

import (
	"bytes"
	"encoding/json"
	"fmt"
	"strconv"
	"strings"
	"time"
)

// OMM holds the CCSDS Orbit Mean-Elements Message fields needed for SGP4.
// Angles are in degrees, mean motion in rev/day, BStar in 1/earth radii.
type OMM struct {
	ObjectName      string
	ObjectID        string
	NoradCatID      int
	Epoch           time.Time // UTC
	MeanMotion      float64
	Eccentricity    float64
	Inclination     float64
	RaOfAscNode     float64
	ArgOfPericenter float64
	MeanAnomaly     float64
	BStar           float64
	MeanMotionDot   float64
	MeanMotionDDot  float64
	ElementSetNo    int
	RevAtEpoch      int
}

// ParseOMM parses OMM JSON (an array or a single object) served as numeric values or quoted string values
func ParseOMM(data []byte) ([]OMM, error) {
	data = bytes.TrimSpace(data)
	if len(data) > 0 && data[0] == '{' {
		data = append(append([]byte{'['}, data...), ']')
	}
	var raw []map[string]any
	if err := json.Unmarshal(data, &raw); err != nil {
		return nil, fmt.Errorf("error parsing OMM JSON: %v", err)
	}

	out := make([]OMM, len(raw))
	for i, m := range raw {
		var errs []error
		f := func(k string) float64 {
			v, err := ommFloat(m, k)
			if err != nil {
				errs = append(errs, err)
			}
			return v
		}
		n := func(k string) int { return int(f(k)) }
		s := func(k string) string { v, _ := m[k].(string); return v }

		epoch, err := time.Parse("2006-01-02T15:04:05.999999999", strings.TrimSuffix(s("EPOCH"), "Z"))
		if err != nil {
			return nil, fmt.Errorf("OMM record %d: error parsing EPOCH: %v", i, err)
		}
		out[i] = OMM{
			ObjectName:      s("OBJECT_NAME"),
			ObjectID:        s("OBJECT_ID"),
			NoradCatID:      n("NORAD_CAT_ID"),
			Epoch:           epoch,
			MeanMotion:      f("MEAN_MOTION"),
			Eccentricity:    f("ECCENTRICITY"),
			Inclination:     f("INCLINATION"),
			RaOfAscNode:     f("RA_OF_ASC_NODE"),
			ArgOfPericenter: f("ARG_OF_PERICENTER"),
			MeanAnomaly:     f("MEAN_ANOMALY"),
			BStar:           f("BSTAR"),
			MeanMotionDot:   f("MEAN_MOTION_DOT"),
			MeanMotionDDot:  f("MEAN_MOTION_DDOT"),
			ElementSetNo:    n("ELEMENT_SET_NO"),
			RevAtEpoch:      n("REV_AT_EPOCH"),
		}
		if len(errs) > 0 {
			return nil, fmt.Errorf("OMM record %d: %v", i, errs[0])
		}
		if out[i].MeanMotion <= 0 {
			return nil, fmt.Errorf("OMM record %d: MEAN_MOTION must be positive", i)
		}
	}
	return out, nil
}

// ommFloat reads a numeric OMM field that may be a JSON number, a numeric
// string, or missing/null (0).
func ommFloat(m map[string]any, k string) (float64, error) {
	switch v := m[k].(type) {
	case nil:
		return 0, nil
	case float64:
		return v, nil
	case string:
		f, err := ParseFloat(v)
		if err != nil {
			return 0, fmt.Errorf("error parsing %s: %v", k, err)
		}
		return f, nil
	default:
		return 0, fmt.Errorf("unexpected type for %s: %T", k, v)
	}
}

// TLE converts the OMM to a TLE usable by PropagateSatellite.
func (o OMM) TLE() *TLE {
	e := o.Epoch.UTC()
	secs := float64(e.Hour()*3600+e.Minute()*60+e.Second()) + float64(e.Nanosecond())/1e9
	return &TLE{
		SatelliteNumber: strconv.Itoa(o.NoradCatID),
		CatalogNumber:   o.NoradCatID,
		Epoch:           float64(e.Year()%100*1000+e.YearDay()) + secs/86400.0, // display only
		JulianEpoch:     TimeToJulianDate(e),
		GregorianEpoch:  e,
		MeanMotionDot:   o.MeanMotionDot,
		MeanMotionDDot:  o.MeanMotionDDot, // exponent already applied, Exp = 0
		BStar:           o.BStar,          // exponent already applied, Exp = 0
		ElementSet:      fmt.Sprintf("%03d", o.ElementSetNo),
		Inclination:     o.Inclination,
		RightAscension:  o.RaOfAscNode,
		Eccentricity:    o.Eccentricity,
		ArgumentPerigee: o.ArgOfPericenter,
		MeanAnomaly:     o.MeanAnomaly,
		MeanMotion:      o.MeanMotion,
	}
}
