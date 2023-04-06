package ptrack

import (
	"fmt"
	"log"
)

// Track a collection of particles through the domain
func (d *Domain) TrackParticles(p Particles) ([][]Particle, int) {
	o, c := make([][]Particle, len(p)), 0
	for k, pp := range p {
		pid := d.findStartingPrism(&pp)
		o[k] = d.trackParticle(&pp, pid)
		c += len(o[k])
	}
	return o, c
}

func (d *Domain) trackParticle(p *Particle, pid int) []Particle {
	var pl pathline

	fmt.Printf("  >>> particle %d start point (x,y,z): %6.3f %6.3f %6.3f released in prism %d\n", p.I, p.X, p.Y, p.Z, pid)
	d.trackRecurse(p, &pl, pid, -1)

	// pl = pl[:len(pl)-1]
	plast := pl[len(pl)-1]
	// fmt.Printf("\tparticle exit point  (x,y,z,t): %6.3f %6.3f %6.3f %6.3es\n", plast[0], plast[1], plast[2], plast[3])
	fmt.Printf("\tparticle exit point  (x,y,z,t): %6.3f %6.3f %6.3f %6.3es\n", plast.X, plast.Y, plast.Z, plast.T)
	return pl
}

func (d *Domain) findStartingPrism(p *Particle) int {
	prismIDs := d.ParticleToPrismIDs(p, -1)
	switch len(prismIDs) {
	case 1:
		return prismIDs[0]
	case 0:
		log.Fatalf(" ** ERROR, particle start point (x,y,z): %6.3f %6.3f %6.3f point not in domain\n", p.X, p.Y, p.Z)
		return -1
	default:
		// var a []int
		// // fmt.Printf(" ** warning, start point (x,y,z): %6.3f %6.3f %6.3f sits on boundary of %d cells, may wish to adjust.\n", p.X, p.Y, p.Z, len(prismIDs))
		// for _, pid := range prismIDs {
		// 	vx, vy, vz := d.VF[pid].PointVelocity(p, d.prsms[pid], 0.)            // velocity
		// 	mag := math.Sqrt(math.Pow(vx, 2) + math.Pow(vy, 2) + math.Pow(vz, 2)) // magnitude
		// 	cdist := d.prsms[pid].Centroid().Dist(p) / 100.                       // for scaling
		// 	newp := Particle{X: p.X + vx/mag*cdist, Y: p.Y + vy/mag*cdist, Z: p.Z + vz/mag*cdist}
		// 	newPrismIDs := d.ParticleToPrismIDs(&newp, -1)
		// 	if len(newPrismIDs) != 1 {
		// 		log.Fatalf(" error when attempting to find appropriate starting cell, may need to adjust scaling above")
		// 	}
		// 	if newPrismIDs[0] == pid {
		// 		a = append(a, pid)
		// 	}
		// }
		// if len(a) != 1 {
		// 	log.Fatalf(" error when attempting to find appropriate starting cell, more than one was found")
		// }
		// return a[0]
		return -1
	}
}
