package ptrack

import (
	"fmt"
	"log"
)

// Track a collection of particles through the domain
func (d *Domain) TrackParticles(p Particles, prnt bool) ([][]Particle, int) {
	o, c := make([][]Particle, len(p)), 0
	for k, pp := range p {
		pid := d.findStartingPrism(&pp)
		o[k] = d.trackParticle(&pp, pid, prnt)
		c += len(o[k])
	}
	return o, c
}

func (d *Domain) trackParticle(p *Particle, pid int, prnt bool) []Particle {
	var pl pathline
	if prnt {
		fmt.Printf("  >>> particle %d start point (x,y,z): %6.3f %6.3f %6.3f released in prism %d\n", p.I, p.X, p.Y, p.Z, pid)
	}
	d.trackRecurse(p, &pl, pid, -1, prnt)

	// pl = pl[:len(pl)-1]
	plast := pl[len(pl)-1]

	if prnt {
		fmt.Printf("\tparticle exit point  (x,y,z,t): %6.3f %6.3f %6.3f %6.3es\n", plast.X, plast.Y, plast.Z, plast.T)
	}

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
		return -1
	}
}
