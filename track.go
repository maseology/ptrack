package ptrack

import (
	"fmt"
)

// Track a collection of particles through the domain
func (d *Domain) Track(p Particles, pt ParticleTracker) [][][]float64 {
	o := make([][][]float64, len(p))
	for k, pp := range p {
		fmt.Printf(" %d: particle start point (x,y,z,t): %6.3f %6.3f %6.3f\n", k, pp.X, pp.Y, pp.Z)
		o[k] = d.trackParticle(&pp, &pt)
	}
	return o
}

func (d *Domain) trackParticle(p *Particle, pt *ParticleTracker) [][]float64 {
	var pathline [][]float64
	for _, pid := range d.ParticleToPrismIDs(p) {
		pathline = make([][]float64, 0)

		d.trackRecurse(p, pt, &pathline, pid, -1)

		if len(pathline) > 10 {
			break // in cases where the startpoint sits on a boundary, return the first that tracked
		}
	}

	plast := pathline[len(pathline)-1]
	fmt.Printf("\tparticle exit point  (x,y,z,t): %6.3f %6.3f %6.3f %6.3es\n", plast[0], plast[1], plast[2], plast[3])

	return pathline
}

func (d *Domain) trackRecurse(p *Particle, pt *ParticleTracker, pl *[][]float64, i, il int) {

	// track within prism
	plt := trackToExit(p, d.prsms[i], d.VF[i], *pt, d.zc[i])
	for _, p := range plt {
		*pl = append(*pl, p) // add tracks
	}

	pids := d.ParticleToPrismIDs(p)
	switch len(pids) {
	case 0:
		fmt.Printf("\tparticle has exited domain %d\n", i)
	case 1:
		if pids[0] == il {
			if i == il {
				fmt.Printf("\ttracking aborted where particle exited well in cell %d\n", i)
			} else {
				fmt.Printf("\ttracking aborted where particle cycle occurred between cells %d-%d\n", i, il)
			}
		} else {
			d.trackRecurse(p, pt, pl, pids[0], i)
		}
	default:
		fmt.Printf(" particle path has diverged %d\n", i)
		d.trackRecurse(p, pt, pl, pids[0], i)
	}
}
