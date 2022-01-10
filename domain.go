package ptrack

import (
	"fmt"
	"math"
	"math/cmplx"
	"sort"

	"github.com/maseology/mmio"
)

// Domain is a set of cells that constitute a model
type Domain struct {
	prsms  map[int]*Prism
	flx    map[int]*PrismFlux
	VF     map[int]VelocityFielder
	conn   map[int][]int
	zw     map[int]complex128
	extent []float64
}

// Nprism returns the prisms (cells) in the domain
func (d *Domain) Nprism() int { return len(d.prsms) }

// New Domain constructor
func (d *Domain) New(pset PrismSet, pflxs map[int]*PrismFlux) {
	d.prsms = pset.P
	d.conn = pset.Conn
	d.flx = pflxs
	d.zw = make(map[int]complex128, len(d.prsms))
	for i, q := range d.prsms {
		zwt := complex(q.CentroidXY())
		if d.flx[i].qw == 0. {
			zwt = cmplx.NaN()
		}
		d.zw[i] = zwt
	}
	zn, zx, yn, yx, xn, xx := d.getExtent()
	d.extent = []float64{zn, zx, yn, yx, xn, xx}
}

// MakeWaterloo creates velocity field using the Waterloo Method
func (d *Domain) MakeWaterloo() {
	fmt.Println(" building Waterloo method flow field..")
	d.VF = make(map[int]VelocityFielder, len(d.prsms))
	for i, q := range d.prsms {
		var wm WatMethSoln
		ql, qb, qt := d.flx[i].LatBotTop() // left-up-right-down-bottom-top
		wm.New(i, q, ql, d.zw[i], -qt, qb, d.flx[i].qw, 80, 30, false)
		d.VF[i] = &wm
	}
}

// MakePollock creates velocity field using the Pollock (MODPATH) Method
func (d *Domain) MakePollock(dt float64) *PollockMethod {
	fmt.Println(" building Pollock method flow field..")
	d.VF = make(map[int]VelocityFielder, len(d.prsms))
	for i, q := range d.prsms {
		var pm PollockMethod
		ql, qb, qt := d.flx[i].LatBotTop()                            // left-up-right-down-bottom-top
		pm.New(q, d.zw[i], ql[0], -ql[2], ql[3], -ql[1], qb, -qt, dt) // q (prism), well (assumed centroid), Qx0, Qx1, Qy0, Qy1, Qz0, Qz1,  dt
		d.VF[i] = &pm
	}
	return &PollockMethod{}
}

// MakeVector creates velocity field based on a uniform prism velocity vector
func (d *Domain) MakeVector() *VectorMethSoln {
	fmt.Println(" building vector-based flow field..")
	d.VF = make(map[int]VelocityFielder, len(d.prsms))
	for i, q := range d.prsms {
		var vm VectorMethSoln
		r, zc := 0., 0+0.i
		for _, c := range q.Z {
			zc += c
		}
		zc /= complex(float64(len(q.Z)), 0)
		for j := 0; j < len(q.Z); j++ {
			r = math.Max(r, cmplx.Abs(q.Z[j]-zc))
		}
		vm.New(zc, d.flx[i].q, q.Por, r)
		d.VF[i] = &vm
	}
	return &VectorMethSoln{}
}

// Print properties of the domain
func (d *Domain) Print() {
	zn, zx, yn, yx, xn, xx := d.getExtent()
	fmt.Printf("  Nprism: %d\n  Extent X = [%.1f, %.1f]; Y = [%.1f, %.1f]; Z = [%.1f, %.1f]\n", d.Nprism(), xn, xx, yn, yx, zn, zx)
}

func (d *Domain) ResetTops(v float64) {
	for _, p := range d.prsms {
		p.Top = v
	}
}

func (d *Domain) PrismTopsToSaturated() {
	for _, p := range d.prsms {
		p.Top = p.Bn
	}
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
func (d *Domain) ParticleToPrismIDs(p *Particle, prsmLast int) []int {
	var pids []int
	if prsmLast < 0 { // brute force solution
		for i, r := range d.prsms {
			if r.Contains(p) {
				pids = append(pids, i)
			}
		}
	} else {
		for _, pid := range d.conn[prsmLast] {
			if pid < 0 { //left-up-right-down-bottom-top
				continue
			}
			if d.prsms[pid].Contains(p) {
				pids = append(pids, pid)
			}
		}
	}
	return pids
}

// // ParticleToPrismID returns the first prism for which a particle is located
// func (d *Domain) ParticleToPrismID(p *Particle) int {
// 	// brute force solution
// 	for i, r := range d.prsms {
// 		if r.Contains(p) {
// 			return i
// 		}
// 	}
// 	return -1
// }

func (d *Domain) PrintToCSV(fp string) {
	csvw := mmio.NewCSVwriter(fp)
	csvw.WriteHead("pid,conn[left-up-right-down-bottom-top],dim[range(x y z)],maxQ")
	pi := make([]int, 0, len(d.prsms))
	for i := range d.prsms {
		pi = append(pi, i)
	}
	sort.Ints(pi)
	for _, i := range pi {
		yn, yx, xn, xx := d.prsms[i].getExtentsXY()
		csvw.WriteLine(i, fmt.Sprintf("%v", d.conn[i]), fmt.Sprintf("[(%.1f %.1f) (%.1f %.1f) (%.1f %.1f)]", xn, xx, yn, yx, d.prsms[i].Bot, d.prsms[i].Top), d.flx[i].MaxAbsFlux())
	}

	csvw.Close()
}
