package ptrack

import (
	"fmt"
	"math"
	"math/cmplx"
)

const (
	ds = 1.
	dt = 1.
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

// Nprism returns the prisms (cells) in the domain
func (d *Domain) Nprism() int { return len(d.prsms) }

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

// Print properties of the domain
func (d *Domain) Print() {
	zn, zx, yn, yx, xn, xx := d.Extent()
	fmt.Printf("  Nprism: %d\n  Extent X = [%.1f, %.1f]; Y = [%.1f, %.1f]; Z = [%.1f, %.1f]\n", d.Nprism(), xn, xx, yn, yx, zn, zx)
}

// Extent returns the global extent of all coordinates [yn, yx, xn, xx]
func (d *Domain) Extent() (zn, zx, yn, yx, xn, xx float64) {
	zn, zx, yn, yx, xn, xx = math.MaxFloat64, -math.MaxFloat64, math.MaxFloat64, -math.MaxFloat64, math.MaxFloat64, -math.MaxFloat64
	for _, v := range d.prsms {
		yn1, yx1, xn1, xx1 := v.ExtentsXY()
		yn = math.Min(yn, yn1)
		yx = math.Max(yx, yx1)
		xn = math.Min(xn, xn1)
		xx = math.Max(xx, xx1)
		zn = math.Min(zn, v.Bot)
		zx = math.Max(zx, v.Top)
	}
	return
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
