package ptrack

import (
	"math"
	"math/cmplx"
)

// testTracktoExit track particle to the next point
func (rk *RungeKutta) TestTracktoExit(p *Particle, q *Prism, w VelocityFielder) []Particle {
	var aout []Particle
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
		trial(p, q, w, rk.Dt)
		aout = append(aout, *p)
	}
	return aout
}

// testTracktoExit track particle to the next point
func (rk *RungeKuttaAdaptive) TestTracktoExit(p *Particle, q *Prism, w VelocityFielder) []Particle {
	var aout []Particle
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
		p0, p1 := *p, *p

		// 1 full step
		if trial(&p1, q, w, rk.Dt) {
			rk.Dt /= 2.
			continue
		}
		// 2 half steps
		if trial(&p0, q, w, rk.Dt/2.) {
			rk.Dt /= 2.
			continue
		}
		if trial(&p0, q, w, rk.Dt/2.) {
			rk.Dt /= 2.
			continue
		}

		dst := p0.Dist(&p1)
		if dst == 0. {
			rk.Dt *= 2.
		} else {
			rk.Dt *= .9 * math.Pow(rk.Ds/dst, .2) // adaptive timestepping
		}

		if dst > rk.Ds {
			continue // time step too large, repeat calculation
		}

		// update particle state
		p.X = p0.X
		p.Y = p0.Y
		p.Z = p0.Z
		p.T = p0.T
		aout = append(aout, *p)
	}
	return aout
}
