package ptrack

import (
	"log"
	"math"
	"math/cmplx"
)

// PollockMethod is the solution to the pollock method
// Pollock, D.W., 1989, Documentation of a computer program to compute and display pathlines using results from the U.S. Geological Survey modular three-dimensional finite-difference ground-water flow model: U.S. Geological Survey Open-File Report 89–381.
type PollockMethod struct {
	zc                           complex128
	r, x0, y0, z0                float64
	vx0, vx1, vy0, vy1, vz0, vz1 float64
	ax, ay, az                   float64
	dt                           float64
}

// New PollockMethod constructor
func (pm *PollockMethod) New(q *Prism, Qx0, Qx1, Qy0, Qy1, Qz0, Qz1, Qwell float64) {
	if len(q.Z) != 4 {
		log.Fatalf("PollockMethod can only be used with rectilinear cells")
	}

	// cell dimensions ' xn:left, xx:right, yn:front, yx:back, zn:bottom, zx:top
	yn, yx, xn, xx := q.getExtentsXY()
	dx, dy, dz := xx-xn, yx-yn, q.Top-q.Bot
	pm.r = math.Max(dy, dx)

	for _, c := range q.Z {
		pm.zc += c
	}
	pm.zc /= complex(float64(4), 0)

	// face velocities (normal, positive right-back-up)
	pm.vx0, pm.vx1 = Qx0/q.Por/dy/dz, Qx1/q.Por/dy/dz
	pm.vy0, pm.vy1 = Qy0/q.Por/dx/dz, Qy1/q.Por/dx/dz
	pm.vz0, pm.vz1 = Qz0/q.Por/dx/dy, Qz1/q.Por/dx/dy

	pm.ax = (pm.vx1 - pm.vx0) / dx
	pm.ay = (pm.vy1 - pm.vy0) / dy
	pm.az = (pm.vz1 - pm.vz0) / dz
}

// PointVelocity returns the velocity vector for a given (x,y,z) coordinate. (must only be used with rectilinear cells)
func (pm *PollockMethod) PointVelocity(p *Particle, q *Prism, dummy float64) (float64, float64, float64) {
	return pm.ax*(p.X-pm.x0) + pm.vx0, pm.ay*(p.Y-pm.y0) + pm.vy0, pm.az*(p.Z-pm.z0) + pm.vz0
}

// Contains returns whether the point is solvable within the solution space
func (pm *PollockMethod) Contains(p *Particle) (float64, bool) {
	azl := cmplx.Abs(complex(p.X, p.Y)-pm.zc) / pm.r // relative coordinate
	return azl, azl <= 1.
}

func (pm *PollockMethod) exitTime(vx, vy, vz float64) float64 {
	texit := func(v0, v1, v, a float64) float64 {
		if v != 0. && a > 0. {
			te0, te1 := math.MaxFloat64, math.MaxFloat64
			if v0 != 0. {
				te0 = math.Log(math.Abs(v0/v)) / a
			}
			if v1 != 0. {
				te1 = math.Log(math.Abs(v1/v)) / a
			}
			return math.Min(te0, te1)
		}
		return math.MaxFloat64
	}

	txe := texit(pm.vx0, pm.vx1, vx, pm.ax)
	tye := texit(pm.vy0, pm.vy1, vy, pm.ay)
	tze := texit(pm.vz0, pm.vz1, vz, pm.az)

	return math.Min(tze, math.Min(tye, txe))
}

// Track particle through prism
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
				te := pm.exitTime(vx, vy, vz)
				xe := pm.x0 + (vx*math.Exp(pm.ax*te)-pm.vx0)/pm.ax
				ye := pm.y0 + (vy*math.Exp(pm.ay*te)-pm.vy0)/pm.ay
				ze := pm.z0 + (vz*math.Exp(pm.az*te)-pm.vz0)/pm.az

				// update particle locations
				if te > pm.dt { // iterate
					t := 0.
					for {
						t += pm.dt
						if t >= te {
							break
						}
						vx, vy, vz = vf.PointVelocity(p, q, 0.)
						p.X = pm.x0 + (vx*math.Exp(pm.ax*pm.dt)-pm.vx0)/pm.ax
						p.Y = pm.y0 + (vy*math.Exp(pm.ay*pm.dt)-pm.vy0)/pm.ay
						p.Z = pm.z0 + (vz*math.Exp(pm.az*pm.dt)-pm.vz0)/pm.az
						p.T += pm.dt
						chout <- p.State()
					}
				}
				// exit point (plus a nudge needed to get out of cell)
				p.X = xe //+ vx*pm.dt/1e3
				p.Y = ye //+ vy*pm.dt/1e3
				p.Z = ze //+ vz*pm.dt/1e3
				p.T += te
				chout <- p.State()
				return
			}
		}
	}()
	return chout
}
