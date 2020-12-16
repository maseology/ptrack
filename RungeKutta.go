package ptrack

import (
	"math"
)

const (
	safetyFactor = .9
	dtx          = 86400.
)

// RungeKutta particle pathline integration scheme
type RungeKutta struct{ Dt float64 }

// RungeKuttaAdaptive particle pathline integration scheme
type RungeKuttaAdaptive struct{ Ds, Dt float64 }

// // TrackToExit tracks particle to the exit point
// func (rk *RungeKutta) TrackToExit(p *Particle, q *Prism, w VelocityFielder) [][]float64 {
// 	output := make([][]float64, 1)
// 	output[0] = p.State()
// nextstep:
// 	rk.Track(p, q, w)
// 	if len(output) == 1 || q.Contains(p) {
// 		output = append(output, p.State())
// 		goto nextstep
// 	}
// 	// somewhere between current point and the last is the prism bound
// 	q.Intersection(p, output[len(output)-1])

// 	return output
// }

// Track track particle to the next point
func (rk *RungeKutta) Track(p *Particle, q *Prism, w VelocityFielder) { trial(p, q, w, rk.Dt) }

// Track track particle to the next point
func (rk *RungeKuttaAdaptive) Track(p *Particle, q *Prism, w VelocityFielder) {
redo:
	p0, p1 := *p, *p

	// 2 half steps
	trial(&p0, q, w, rk.Dt/2.)
	trial(&p0, q, w, rk.Dt/2.)
	// 1 full step
	trial(&p1, q, w, rk.Dt)

	dst := p0.Dist(&p1)
	rk.Dt *= safetyFactor * math.Pow(rk.Ds/dst, .2) // adaptive timestepping
	if math.IsNaN(rk.Dt) {
		print()
	}

	if dst > rk.Ds {
		goto redo // time step too large, repeat calculation
	}

	if rk.Dt > dtx {
		rk.Dt = dtx
	}
	// fmt.Printf(" dt: %15.5f\n", rk.Dt)

	// update particle state
	p.X = p0.X
	p.Y = p0.Y
	p.Z = p0.Z
	p.T = p0.T
}

func trial(p *Particle, q *Prism, w VelocityFielder, dt float64) {

	if p.X > 6000. {
		print()
	}
	vx, vy, vz := w.PointVelocity(p, q, 0.)
	if vx < 0. {
		print()
	}
	k1 := dt * vx
	l1 := dt * vy
	m1 := dt * vz

	p2 := Particle{0, p.X + k1/2., p.Y + l1/2., p.Z + m1/2., p.T + dt/2.}
	if p2.X > 6000. {
		print()
	}
	vx, vy, vz = w.PointVelocity(&p2, q, 0.)
	k2 := dt * vx
	l2 := dt * vy
	m2 := dt * vz

	p3 := Particle{0, p.X + k2/2., p.Y + l2/2., p.Z + m2/2., p.T + dt/2.}
	if p3.X > 6000. {
		print()
	}
	vx, vy, vz = w.PointVelocity(&p3, q, 0.)
	k3 := dt * vx
	l3 := dt * vy
	m3 := dt * vz

	p4 := Particle{0, p.X + k3, p.Y + l3, p.Z + m3, p.T + dt}
	if p4.X > 6000. {
		print()
	}
	vx, vy, vz = w.PointVelocity(&p4, q, 0.)
	k4 := dt * vx
	l4 := dt * vy
	m4 := dt * vz

	p.X += (k1 + 2.*k2 + 2.*k3 + k4) / 6.
	p.Y += (l1 + 2.*l2 + 2.*l3 + l4) / 6.
	p.Z += (m1 + 2.*m2 + 2.*m3 + m4) / 6.
	p.T += dt

}
