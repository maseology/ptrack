package ptrack

import (
	"log"
	"math"
	"math/cmplx"
)

const reallyBig = 9.999e100

// PollockMethod is the solution to the pollock method
// Pollock, D.W., 1989, Documentation of a computer program to compute and display pathlines using results from the U.S. Geological Survey modular three-dimensional finite-difference ground-water flow model: U.S. Geological Survey Open-File Report 89–381.
// Pollock, D.W., 2016, User guide for MODPATH Version 7—A particle-tracking model for MODFLOW: U.S. Geological Survey Open-File Report 2016–1086, 35 p., http://dx.doi.org/10.3133/ofr20161086.
type PollockMethod struct {
	zc, zwl                      complex128
	x0, y0, z0, x1, y1, z1       float64
	vx0, vx1, vy0, vy1, vz0, vz1 float64
	r, ax, ay, az                float64
	dt                           float64
}

// New PollockMethod constructor
func (pm *PollockMethod) New(q *Prism, zwl complex128, Qx0, Qx1, Qy0, Qy1, Qz0, Qz1, dt float64) {
	if len(q.Z) != 4 {
		log.Fatalf("PollockMethod can only be used with rectilinear cells")
	}

	// cell dimensions ' xn:left, xx:right, yn:front, yx:back, zn:bottom, zx:top
	yn, yx, xn, xx := q.getExtentsXY()
	dx, dy, dz := xx-xn, yx-yn, q.Top-q.Bot
	pm.r = math.Max(dy, dx)
	pm.x0 = xn
	pm.y0 = yn
	pm.z0 = q.Bot
	pm.x1 = xx
	pm.y1 = yx
	pm.z1 = q.Top
	pm.dt = dt

	for _, c := range q.Z {
		pm.zc += c
	}
	pm.zc /= complex(float64(4), 0)
	pm.zwl = zwl

	// face velocities (normal, positive right-back-up)
	pm.vx0, pm.vx1 = Qx0/q.Por/dy/dz, Qx1/q.Por/dy/dz
	pm.vy0, pm.vy1 = Qy0/q.Por/dx/dz, Qy1/q.Por/dx/dz
	pm.vz0, pm.vz1 = Qz0/q.Por/dx/dy, Qz1/q.Por/dx/dy

	pm.ax = (pm.vx1 - pm.vx0) / dx
	pm.ay = (pm.vy1 - pm.vy0) / dy
	pm.az = (pm.vz1 - pm.vz0) / dz
}

// PointVelocity returns the velocity vector for a given (x,y,z) coordinate. (must only be used with rectilinear cells)
func (pm *PollockMethod) PointVelocity(p *Particle, d1 *Prism, d2 float64) (float64, float64, float64) {
	return pm.ax*(p.X-pm.x0) + pm.vx0, pm.ay*(p.Y-pm.y0) + pm.vy0, pm.az*(p.Z-pm.z0) + pm.vz0
}

func (pm *PollockMethod) exitTime(p *Particle, vx, vy, vz float64) float64 { // exit time
	texit := func(v0, v1, v, s0, s1, s, a float64) float64 {
		// see figure 3 in Pollock, D.W., 2016, User guide for MODPATH Version 7—A particle-tracking model for MODFLOW: U.S. Geological Survey Open-File Report 2016–1086, 35 p., http://dx.doi.org/10.3133/ofr20161086.
		if v == 0. { // case A will not exit
			return reallyBig
		} else if v0 >= 0. && v1 <= 0. { // case A, will not exit
			return reallyBig
		} else if v0 <= 0. && v1 >= 0 { // case B, flow divide
			if v < 0. {
				return math.Log(v0/v) / a
			} else if v > 0. {
				return math.Log(v1/v) / a
			}
			return reallyBig
		} else if v0 == v1 && v0 != 0. { // case C, constant velocity
			if v < 0. {
				return (s0 - s) / v
			}
			return (s1 - s) / v
		}
		if v < 0. {
			return math.Log(v0/v) / a // equation 9
		} else if v > 0. {
			return math.Log(v1/v) / a
		}
		return reallyBig
	}

	txe := texit(pm.vx0, pm.vx1, vx, pm.x0, pm.x1, p.X, pm.ax)
	tye := texit(pm.vy0, pm.vy1, vy, pm.y0, pm.y1, p.Y, pm.ay)
	tze := texit(pm.vz0, pm.vz1, vz, pm.z0, pm.z1, p.Z, pm.az)

	return math.Min(tze, math.Min(tye, txe))
}

func updatePostition(s0, v0, v, a, t float64) float64 {
	if a == 0. {
		return 0.
	}
	return s0 + (v*math.Exp(a*t)-v0)/a
}

// Local returns whether the point is solvable within the solution space
func (pm *PollockMethod) Local(p *Particle) (float64, bool) {
	azl := cmplx.Abs(complex(p.X, p.Y)-pm.zc) / pm.r // relative coordinate
	return azl, azl <= 1.
}

// track (to exit) particle through prism
func (pm *PollockMethod) track(done <-chan interface{}, p *Particle, q *Prism, vf VelocityFielder) <-chan []float64 {
	chout := make(chan []float64)
	go func() {
		defer close(chout)
		for {
			select {
			case <-done:
				return
			default:
				vx, vy, vz := vf.PointVelocity(p, q, 0.)
				te := pm.exitTime(p, vx, vy, vz) * 1.00001 // adding a little momentum to nudge the particle past the boundary
				xe := updatePostition(pm.x0, pm.vx0, vx, pm.ax, te)
				ye := updatePostition(pm.y0, pm.vy0, vy, pm.ay, te)
				ze := updatePostition(pm.z0, pm.vz0, vz, pm.az, te)

				// update particle locations
				if te > pm.dt { // iterate
					t := 0.
					for {
						t += pm.dt
						if t >= te {
							break
						}
						vx, vy, vz = vf.PointVelocity(p, q, 0.)
						p.X = updatePostition(pm.x0, pm.vx0, vx, pm.ax, pm.dt)
						p.Y = updatePostition(pm.y0, pm.vy0, vy, pm.ay, pm.dt)
						p.Z = updatePostition(pm.z0, pm.vz0, vz, pm.az, pm.dt)
						p.T += pm.dt
						chout <- p.State()
					}
				}

				// exit point
				p.X = xe
				p.Y = ye
				p.Z = ze
				p.T += te
				chout <- p.State()
				return
			}
		}
	}()
	return chout
}
