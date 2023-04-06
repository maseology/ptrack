package ptrack

import "math"

// EulerSpace constant space step Euler particle pathline integration scheme
type EulerSpace struct{ Ds float64 }

// EulerTime constant time step Euler particle pathline integration scheme
type EulerTime struct{ Dt float64 }

// Track track particle to the next space step
func (es *EulerSpace) track(done <-chan interface{}, p *Particle, q *Prism, w VelocityFielder) <-chan Particle {
	chout := make(chan Particle)
	go func() {
		defer close(chout)
		for {
			select {
			case <-done:
				return
			default:
				vx, vy, vz := w.PointVelocity(p, q, 0.)
				dt := es.Ds / math.Sqrt(math.Pow(vx, 2.)+math.Pow(vy, 2.)+math.Pow(vz, 2.)) // ds/|vn|
				p.X += vx * dt
				p.Y += vy * dt
				p.Z += vz * dt
				p.T += dt
				chout <- *p
			}
		}
	}()
	return chout
}

// Track track particle to the next time step
func (et *EulerTime) track(done <-chan interface{}, p *Particle, q *Prism, w VelocityFielder) <-chan Particle {
	chout := make(chan Particle)
	go func() {
		defer close(chout)
		for {
			select {
			case <-done:
				return
			default:
				vx, vy, vz := w.PointVelocity(p, q, 0.)
				p.X += vx * et.Dt
				p.Y += vy * et.Dt
				p.Z += vz * et.Dt
				p.T += et.Dt
				chout <- *p
			}
		}
	}()
	return chout
}
