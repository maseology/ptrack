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
	xn, xx, yn, yx, r            float64
	x0, x1, y0, y1, z0, z1       float64
	vx0, vx1, vy0, vy1, vz0, vz1 float64
	ax, ay, az                   float64
	dt                           float64
}

// New PollockMethod constructor
func (pm *PollockMethod) New(q *Prism, Qx0, Qx1, Qy0, Qy1, Qz0, Qz1, Qwell float64) {
	if len(q.Z) != 4 {
		log.Fatalf("PollockMethod can only be used with rectilinear cells")
	}

	pm.yn, pm.yx, pm.xn, pm.xx = q.getExtentsXY()
	pm.r = math.Max((pm.yx - pm.yn), (pm.xx - pm.xn))

	for _, c := range q.Z {
		pm.zc += c
	}
	pm.zc /= complex(float64(4), 0)

	// cell dimensions ' xn:left, xx:right, yn:front, yx:back, zn:bottom, zx:top
	yn, yx, xn, xx := q.getExtentsXY()
	dx, dy, dz := xx-xn, yx-yn, q.Top-q.Bot

	// face velocities (normal, positive right-back-up)
	pm.vx0, pm.vx1 = Qx0/q.Por/dy/dz, Qx1/q.Por/dy/dz
	pm.vy0, pm.vy1 = Qy0/q.Por/dx/dz, Qy1/q.Por/dx/dz
	pm.vz0, pm.vz1 = Qz0/q.Por/dx/dy, Qz1/q.Por/dx/dy
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

// Track particle through prism
func (pm *PollockMethod) Track(p *Particle, q *Prism, vf VelocityFielder) {
	vx, vy, vz := pm.PointVelocity(p, q, 0.)

	txe, tye, tze := math.Log(pm.vx1/vx)/pm.ax, math.Log(pm.vy1/vy)/pm.ay, math.Log(pm.vz0/vz)/pm.az
	te := math.Min(tze, math.Min(tye, txe))
	xe := pm.x0 + (vx*math.Exp(pm.ax*te)-pm.vx0)/pm.ax
	ye := pm.y0 + (vy*math.Exp(pm.ay*te)-pm.vy0)/pm.ay
	ze := pm.z0 + (vz*math.Exp(pm.az*te)-pm.vz0)/pm.az

	// update particle locations
	if te > pm.dt {
		t := 0.
		mp := map[float64][]float64{0.: p}
		for {
			t += pm.dt
			if t > te {
				break
			}
			p[0] = x1 + (v[0]*math.Exp(Ax*delt)-vx1)/Ax
			p[1] = y1 + (v[1]*math.Exp(Ay*delt)-vy1)/Ay
			p[2] = z1 + (v[2]*math.Exp(Az*delt)-vz1)/Az

			mp[t] = p

			// update particle velocity
			v = []float64{Ax*(p[0]-x0) + vx0, Ay*(p[1]-y0) + vy0, Az*(p[2]-z0) + vz0}
		}
	}

}

func pollockCube(Qx0, Qx1, Qy0, Qy1, Qz0, Qz1 float64) map[float64][]float64 {

	// velocity gradients
	Ax, Ay, Az := (vx1-vx0)/dx, (vy1-vy0)/dy, (vz1-vz0)/dz

	// particle location
	p := []float64{-0.5, -0.35, 1.0}

	// particle velocity
	v := []float64{Ax*(p[0]-x0) + vx0, Ay*(p[1]-y0) + vy0, Az*(p[2]-z0) + vz0}

	dtExit := math.Log(vz0/v[2]) / Az

	xe := x0 + (v[0]*math.Exp(Ax*dtExit)-vx0)/Ax
	ye := y0 + (v[1]*math.Exp(Ay*dtExit)-vy0)/Ay
	ze := z0 + (v[2]*math.Exp(Az*dtExit)-vz0)/Az

	// update particle locations
	t, delt := .0, 0.0001
	mp := map[float64][]float64{0.: p}
	for {
		t += delt
		if t > dtExit {
			break
		}
		p[0] = x1 + (v[0]*math.Exp(Ax*delt)-vx1)/Ax
		p[1] = y1 + (v[1]*math.Exp(Ay*delt)-vy1)/Ay
		p[2] = z1 + (v[2]*math.Exp(Az*delt)-vz1)/Az

		mp[t] = p

		// update particle velocity
		v = []float64{Ax*(p[0]-x0) + vx0, Ay*(p[1]-y0) + vy0, Az*(p[2]-z0) + vz0}
	}
	mp[dtExit] = []float64{xe, ye, ze}
	return mp
}

// func pollockCubeEndpoint(Qx0, Qx1, Qy0, Qy1, Qz0, Qz1 float64) []float64 {

// 	// '=========================
// 	// ' Pollock method
// 	// '=========================
// 	// ' cell dimensions ' x0:left, x1:right, y0:front, y1:back, z0:bottom, z1:top
// 	// Dim x0 = -0.5, x1 = 0.5, y0 = -0.5, y1 = 0.5, z0 = 0.0, z1 = 1.0
// 	// Dim dx = x1 - x0, dy = y1 - y0, dz = z1 - z0

// 	// ' face velocities (normal, postive right-back-up
// 	// Dim porosity = 0.3
// 	// Dim vx0 = Qx0 / porosity / dy / dz, vx1 = Qx1 / porosity / dy / dz
// 	// Dim vy0 = Qy0 / porosity / dx / dz, vy1 = Qy1 / porosity / dx / dz
// 	// Dim vz0 = Qz0 / porosity / dx / dy, vz1 = Qz1 / porosity / dx / dy

// 	// ' velcity gradients
// 	// Dim Ax = (vx1 - vx0) / dx, Ay = (vy1 - vy0) / dy, Az = (vz1 - vz0) / dz

// 	// ' particle location
// 	// Dim p = {-0.5, -0.35, 1.0}

// 	// ' particle velocity
// 	// Dim v = {Ax * (p(0) - x0) + vx0, Ay * (p(1) - y0) + vy0, Az * (p(2) - z0) + vz0}

// 	// Dim dtx_exit = Math.Log(vx1 / v(0)) / Ax
// 	// Dim dty_exit = Math.Log(vy1 / v(1)) / Ay
// 	// Dim dtz_exit = Math.Log(vz0 / v(2)) / Az
// 	// Dim dt_exit = Math.Min(dtz_exit, Math.Min(dty_exit, dtx_exit))

// 	// Dim xe = x0 + (v(0) * Math.Exp(Ax * dt_exit) - vx0) / Ax
// 	// Dim ye = y0 + (v(1) * Math.Exp(Ay * dt_exit) - vy0) / Ay
// 	// Dim ze = z0 + (v(2) * Math.Exp(Az * dt_exit) - vz0) / Az

// 	// Return {xe, ye, ze}
// 	// End Function
// }
