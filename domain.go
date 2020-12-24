package ptrack

import (
	"fmt"
	"math"
	"math/cmplx"
)

// Domain is a set of cells that constitute a model
type Domain struct {
	prsms  map[int]*Prism
	flx    map[int]*PrismFlux
	VF     map[int]VelocityFielder
	conn   map[int][]int
	zc     map[int]complex128
	extent []float64
}

// Nprism returns the prisms (cells) in the domain
func (d *Domain) Nprism() int { return len(d.prsms) }

// New Domain constructor
func (d *Domain) New(pset PrismSet, pflxs map[int]*PrismFlux) {
	d.prsms = pset.P
	d.conn = pset.Conn
	d.flx = pflxs
	d.zc = make(map[int]complex128, len(d.prsms))
	for i, q := range d.prsms {
		zw := complex(q.CentroidXY())
		if d.flx[i].qw == 0. {
			zw = cmplx.NaN()
		}
		d.zc[i] = zw
	}
	zn, zx, yn, yx, xn, xx := d.getExtent()
	d.extent = []float64{zn, zx, yn, yx, xn, xx}
}

// MakeWatMeth crates velocity field using the Waterloo Method
func (d *Domain) MakeWatMeth() {
	d.VF = make(map[int]VelocityFielder, len(d.prsms))
	for i, q := range d.prsms {
		var wm WatMethSoln
		ql, qb, qt := d.flx[i].LatBotTop()
		wm.New(q, ql, d.zc[i], -qt, qb, d.flx[i].qw)
		d.VF[i] = &wm
	}
}

// // MakePollock crates velocity field using the Waterloo Method
// func (d *Domain) MakePollock() {
// 	d.VF = make(map[int]VelocityFielder, len(d.prsms))
// 	for i, q := range d.prsms {
// 		var pm PollockMethod
// 		ql, qb, qt := d.flx[i].LatBotTop()
// 		pm.New(q, ql[0], ql[0], ql[0], ql[0], qb, -qt, d.flx[i].qw)
// 		d.VF[i] = &pm
// 	}
// }

// Print properties of the domain
func (d *Domain) Print() {
	zn, zx, yn, yx, xn, xx := d.getExtent()
	fmt.Printf("  Nprism: %d\n  Extent X = [%.1f, %.1f]; Y = [%.1f, %.1f]; Z = [%.1f, %.1f]\n", d.Nprism(), xn, xx, yn, yx, zn, zx)
}

func (d *Domain) getExtent() (zn, zx, yn, yx, xn, xx float64) {
	zn, zx, yn, yx, xn, xx = math.MaxFloat64, -math.MaxFloat64, math.MaxFloat64, -math.MaxFloat64, math.MaxFloat64, -math.MaxFloat64
	for _, q := range d.prsms {
		yn1, yx1, xn1, xx1 := q.getExtentsXY()
		yn = math.Min(yn, yn1)
		yx = math.Max(yx, yx1)
		xn = math.Min(xn, xn1)
		xx = math.Max(xx, xx1)
		zn = math.Min(zn, q.Bot)
		zx = math.Max(zx, q.Top)
	}
	return
}

// ParticleToPrismIDs returns a set of prisms for which a particle is located
func (d *Domain) ParticleToPrismIDs(p *Particle) []int {
	// brute force solution
	var ii []int
	for i, r := range d.prsms {
		if r.Contains(p) {
			ii = append(ii, i)
		}
	}
	return ii
}

// ParticleToPrismID returns the first prism for which a particle is located
func (d *Domain) ParticleToPrismID(p *Particle) int {
	// brute force solution
	for i, r := range d.prsms {
		if r.Contains(p) {
			return i
		}
	}
	return -1
}
