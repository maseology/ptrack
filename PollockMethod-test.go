package ptrack

import (
	"math"
)

func (pm *PollockMethod) TestTracktoExit(p *Particle, q *Prism, vf VelocityFielder) [][]float64 {
	var aout [][]float64
	vx, vy, vz := vf.PointVelocity(p, q, 0.)
	te := pm.exitTime(p, vx, vy, vz) * 1.00001 // adding a little "momentum" to nudge the particle past the boundary
	xe := updatePostition(pm.x0, pm.vx0, vx, pm.ax, te)
	ye := updatePostition(pm.y0, pm.vy0, vy, pm.ay, te)
	ze := updatePostition(pm.z0, pm.vz0, vz, pm.az, te)
	// xe := pm.x0 + (vx*math.Exp(pm.ax*te)-pm.vx0)/pm.ax
	// ye := pm.y0 + (vy*math.Exp(pm.ay*te)-pm.vy0)/pm.ay
	// ze := pm.z0 + (vz*math.Exp(pm.az*te)-pm.vz0)/pm.az

	// update particle locations
	if te < reallyBig {
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
				// p.X = pm.x0 + (vx*math.Exp(pm.ax*pm.dt)-pm.vx0)/pm.ax
				// p.Y = pm.y0 + (vy*math.Exp(pm.ay*pm.dt)-pm.vy0)/pm.ay
				// p.Z = pm.z0 + (vz*math.Exp(pm.az*pm.dt)-pm.vz0)/pm.az
				p.T += pm.dt
				// fmt.Println(p.State())
				aout = append(aout, p.State())
			}
		}
		p.T += te
	} else {
		p.T = math.Sqrt(math.Pow(p.X-xe, 2)+math.Pow(p.Y-ye, 2)+math.Pow(p.Z-ze, 2)) / math.Sqrt(math.Pow(vx, 2)+math.Pow(vy, 2)+math.Pow(vz, 2))
	}
	// exit point (plus a nudge needed to get out of cell)
	p.X = xe //+ vx*pm.dt/1e3
	p.Y = ye //+ vy*pm.dt/1e3
	p.Z = ze //+ vz*pm.dt/1e3
	// p.Z = (q.Top-q.Bot)/2. + q.Bot
	return aout
}
