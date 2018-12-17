package ptrack

import (
	"fmt"
	"log"
	"math/cmplx"
)

// Domain is a set of cells that constitute a model
type Domain struct {
	prsms map[int]*Prism
	flx   map[int]*PrismFlux
	VF    map[int]VelocityFielder
	conn  map[int][]int
	zc    map[int]complex128
	m, n  int
}

// New Domain constructor
func (d *Domain) New(pset PrismSet, pflxs map[int]*PrismFlux) {
	d.prsms = pset.P
	d.conn = pset.Conn
	d.flx = pflxs
	d.VF = make(map[int]VelocityFielder)
	d.zc = make(map[int]complex128)
	d.m = 80
	d.n = 30
	for i, q := range d.prsms {
		zw := q.CentroidXY()
		if d.flx[i].qw == 0. {
			zw = cmplx.NaN()
		}
		d.zc[i] = zw
	}
}

// ParticleToPrismID determines the prism for which a particle is located
func (d *Domain) ParticleToPrismID(p *Particle) int {
	// brute force solution
	for i, r := range d.prsms {
		if r.Contains(p) {
			return i
		}
	}
	return -1
}

// Track a particle through the domain
func (d *Domain) Track(p *Particle) [][]float64 {
	pid := d.ParticleToPrismID(p)
	if pid < 0 {
		log.Fatalln("Track error: particle not found within domain")
	}

	pathline := make([][]float64, 0)
	d.trackRecurse(p, &pathline, pid)

	plast := pathline[len(pathline)-1]
	fmt.Printf(" particle exit point: %6.4f %6.4f %6.4f %6.4f\n", plast[0], plast[1], plast[2], plast[3])

	return pathline
}

func (d *Domain) trackRecurse(p *Particle, pl *[][]float64, i int) {
	fmt.Println(i)
	if _, ok := d.VF[i]; !ok {
		var wm WatMethSoln
		ql, qb, qt := d.flx[i].LatBotTop()
		wm.New(d.prsms[i], ql, d.zc[i], -qt, qb, d.flx[i].qw, d.m, d.n)
		d.VF[i] = &wm
	}

	// pt := RungeKutta{Dt: 0.0001, Ds: 0.0001, Adaptive: false}
	pt := EulerSpace{Ds: 0.001}
	// pt := EulerTime{Dt: 0.0001}

	ec, plt := TrackToExit(p, d.prsms[i], d.VF[i], &pt, d.zc[i])
	for _, p := range plt {
		*pl = append(*pl, p)
	}

	switch {
	case ec == -9999:
		if pid := d.ParticleToPrismID(p); pid < 0 {
			fmt.Println(" particle has exited domain at start point")
		} else {
			fmt.Printf(" error: particle has not appeared to exit prism %d\n", pid)
		}
	case ec < 0:
		fmt.Printf(" particle has exited at well %d\n", -(ec + 1))
	default:
		fmt.Printf(" particle has exited at face %d\n", ec)
		d.trackRecurse(p, pl, d.conn[i][ec])
	}
}
