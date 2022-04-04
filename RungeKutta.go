package ptrack

import (
	"math"
)

// RungeKutta particle pathline integration scheme
type RungeKutta struct{ Dt float64 }

// RungeKuttaAdaptive particle pathline integration scheme
type RungeKuttaAdaptive struct{ Ds, Dt float64 }

// Track track particle to the next point
func (rk *RungeKutta) track(done <-chan interface{}, p *Particle, q *Prism, w VelocityFielder) <-chan []float64 {
	chout := make(chan []float64)
	go func() {
		defer close(chout)
		for {
			select {
			case <-done:
				return
			default:
				trial(p, q, w, rk.Dt)
				chout <- p.State()
			}
		}
	}()
	return chout
}

// Track track particle to the next point
func (rk *RungeKuttaAdaptive) track(done <-chan interface{}, p *Particle, q *Prism, w VelocityFielder) <-chan []float64 {
	chout := make(chan []float64)
	go func() {
		defer close(chout)
		for {
			select {
			case <-done:
				return
			default:
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
				chout <- p.State()
			}
		}
	}()
	return chout
}

func trial(p *Particle, q *Prism, w VelocityFielder, dt float64) bool {
	if r, _ := w.Local(p); r > rmax {
		return true
	}
	vx, vy, vz := w.PointVelocity(p, q, 0.)
	k1 := dt * vx
	l1 := dt * vy
	m1 := dt * vz

	p2 := Particle{0, p.X + k1/2., p.Y + l1/2., p.Z + m1/2., p.T + dt/2.}
	if r, _ := w.Local(&p2); r > rmax {
		return true
	}
	vx, vy, vz = w.PointVelocity(&p2, q, 0.)
	k2 := dt * vx
	l2 := dt * vy
	m2 := dt * vz

	p3 := Particle{0, p.X + k2/2., p.Y + l2/2., p.Z + m2/2., p.T + dt/2.}
	if r, _ := w.Local(&p3); r > rmax {
		return true
	}
	vx, vy, vz = w.PointVelocity(&p3, q, 0.)
	k3 := dt * vx
	l3 := dt * vy
	m3 := dt * vz

	p4 := Particle{0, p.X + k3, p.Y + l3, p.Z + m3, p.T + dt}
	if r, _ := w.Local(&p4); r > rmax {
		return true
	}
	vx, vy, vz = w.PointVelocity(&p4, q, 0.)
	k4 := dt * vx
	l4 := dt * vy
	m4 := dt * vz

	p.X += (k1 + 2.*k2 + 2.*k3 + k4) / 6.
	p.Y += (l1 + 2.*l2 + 2.*l3 + l4) / 6.
	p.Z += (m1 + 2.*m2 + 2.*m3 + m4) / 6.
	p.T += dt

	return false
}
