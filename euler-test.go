package ptrack

import (
	"math"
	"math/cmplx"
)

// testTrack track particle to the next space step
func (es *EulerSpace) TestTracktoExit(p *Particle, q *Prism, w VelocityFielder) [][]float64 {
	var aout [][]float64
	nearwell := func() bool {
		wm := w.(*WatMethSoln)
		if !cmplx.IsNaN(wm.zwl[0]) {
			for _, zwell := range wm.zwl {
				zl := (complex(p.X, p.Y) - wm.zc) / wm.r // complex local coordinate
				if cmplx.Abs(zl-zwell) < wellTol/real(wm.r) {
					return true
				}
			}
		}
		return false
	}
	for {
		if !q.Contains(p) || nearwell() {
			break
		}
		vx, vy, vz := w.PointVelocity(p, q, 0.)
		dt := es.Ds / math.Sqrt(math.Pow(vx, 2.)+math.Pow(vy, 2.)+math.Pow(vz, 2.)) // ds/|vn|
		p.X += vx * dt
		p.Y += vy * dt
		p.Z += vz * dt
		p.T += dt
		aout = append(aout, p.State())
	}
	return aout
}

// testTrack track particle to the next time step
func (et *EulerTime) TestTracktoExit(p *Particle, q *Prism, w VelocityFielder) [][]float64 {
	var aout [][]float64
	nearwell := func() bool {
		wm := w.(*WatMethSoln)
		if !cmplx.IsNaN(wm.zwl[0]) {
			for _, zwell := range wm.zwl {
				zl := (complex(p.X, p.Y) - wm.zc) / wm.r // complex local coordinate
				if cmplx.Abs(zl-zwell) < wellTol/real(wm.r) {
					return true
				}
			}
		}
		return false
	}
	for {
		if !q.Contains(p) || nearwell() {
			break
		}
		vx, vy, vz := w.PointVelocity(p, q, 0.)
		p.X += vx * et.Dt
		p.Y += vy * et.Dt
		p.Z += vz * et.Dt
		p.T += et.Dt
		aout = append(aout, p.State())
	}
	return aout
}
