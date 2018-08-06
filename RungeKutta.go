package ptrack

import "math"

// RungeKutta particle pathline integration scheme
type RungeKutta struct {
	Dt, Ds   float64
	Adaptive bool
}

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
func (rk *RungeKutta) Track(p *Particle, q *Prism, w VelocityFielder) {
	if rk.Adaptive {
	redo:
		p0, p1 := p.Clone(), p.Clone()

		trial(&p0, q, w, rk.Dt/2.)
		trial(&p0, q, w, rk.Dt/2.)
		trial(&p1, q, w, rk.Dt)

		dst := p0.Dist(&p1)
		rk.Dt *= 0.9 * math.Pow(rk.Ds/dst, 0.2) // adaptive timestepping

		if dst > rk.Ds {
			goto redo // time step too large, repeat calculation
		}

		// update particle state
		p.X = p0.X
		p.Y = p0.Y
		p.Z = p0.Z
		p.T = p0.T
	} else {
		trial(p, q, w, rk.Dt)
	}
}

func trial(p *Particle, q *Prism, w VelocityFielder, dt float64) {

	vx, vy, vz := w.PointVelocity(p, q, 0.)
	k1 := dt * vx
	l1 := dt * vy
	m1 := dt * vz
	p1 := Particle{0, p.X + k1/2., p.Y + l1/2., p.Z + m1/2., p.T + dt/2.}

	vx, vy, vz = w.PointVelocity(&p1, q, 0.)
	k2 := dt * vx
	l2 := dt * vy
	m2 := dt * vz
	p2 := Particle{0, p.X + k2/2., p.Y + l2/2., p.Z + m2/2., p.T + dt/2.}

	vx, vy, vz = w.PointVelocity(&p2, q, 0.)
	k3 := dt * vx
	l3 := dt * vy
	m3 := dt * vz
	p3 := Particle{0, p.X + k3, p.Y + l3, p.Z + m3, p.T + dt}

	vx, vy, vz = w.PointVelocity(&p3, q, 0.)
	k4 := dt * vx
	l4 := dt * vy
	m4 := dt * vz

	p.X += (k1 + 2.*k2 + 2.*k3 + k4) / 6.
	p.Y += (l1 + 2.*l2 + 2.*l3 + l4) / 6.
	p.Z += (m1 + 2.*m2 + 2.*m3 + m4) / 6.
	p.T += dt
}
