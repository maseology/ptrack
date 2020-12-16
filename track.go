package ptrack

import (
	"fmt"
	"log"
)

// Track a particle through the domain
func (d *Domain) Track(p Particles, pt ParticleTracker) [][][]float64 {
	o := make([][][]float64, len(p))
	for k, pp := range p {
		fmt.Println(pp)
		o[k] = d.trackParticle(&pp, &pt)
	}
	return o
}

func (d *Domain) trackParticle(p *Particle, pt *ParticleTracker) [][]float64 {
	pid := d.ParticleToPrismID(p)
	if pid < 0 {
		log.Fatalln("Track error: particle not found within domain")
	}

	pathline := make([][]float64, 0)
	d.trackRecurse(p, pt, &pathline, pid)

	plast := pathline[len(pathline)-1]
	fmt.Printf(" particle exit point: %6.4f %6.4f %6.4f %6.4f\n", plast[0], plast[1], plast[2], plast[3])

	return pathline
}

func (d *Domain) trackRecurse(p *Particle, pt *ParticleTracker, pl *[][]float64, i int) {
	fmt.Printf(" > prism %d\n", i)
	if _, ok := d.VF[i]; !ok {
		var wm WatMethSoln
		ql, qb, qt := d.flx[i].LatBotTop()
		wm.New(d.prsms[i], ql, d.zc[i], -qt, qb, d.flx[i].qw, d.m, d.n)
		d.VF[i] = &wm
	}

	ec, plt := TrackToExit(p, d.prsms[i], d.VF[i], *pt, d.zc[i], i)
	for _, p := range plt {
		*pl = append(*pl, p)
	}

	switch {
	case ec == -9999:
		if pid := d.ParticleToPrismID(p); pid < 0 {
			fmt.Println(" particle has exited domain at start point")
		} else {
			yn, yx, xn, xx := d.prsms[pid].ExtentsXY()
			fmt.Printf(" error: particle has not appeared to exit prism %d [%.1f:%.1f, %.1f:%.1f]\n", pid, xn, xx, yn, yx)
		}
	case ec < 0:
		fmt.Printf(" particle has exited at well %d\n", -ec)
	case d.conn[i][ec] < 0:
		fmt.Printf(" particle has exited prism %d\n", i)
	default:
		// fmt.Printf(" particle has exited at face %d\n", ec)
		d.trackRecurse(p, pt, pl, d.conn[i][ec])
	}
}
