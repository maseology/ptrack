package ptrack

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/maseology/mmaths/vector"
)

const ncheck, xuniq, prcsn = 100, 10, .01

type pathline []Particle

var cycl map[int]int

func (d *Domain) trackRecurse(p *Particle, pl *pathline, i, il int, prnt bool) {
	if il < 0 {
		cycl = map[int]int{}
	}
	cycl[i]++
	if cycl[i] > 1 {
		if prnt {
			fmt.Printf("\ttracking aborted where particle cycle occurred at cell %d (%6.3f,%6.3f,%6.3f,%6.3e)\n", i, p.X, p.Y, p.Z, p.T)
		}
		return
	}

	p.C = i
	*pl = append(*pl, *p)

	// check for cycles
	unique := func(s []complex128) []complex128 {
		keys := make(map[complex128]bool)
		list := []complex128{}
		for _, ss := range s {
			ss = complex(math.Round(real(ss)/prcsn)*prcsn, math.Round(imag(ss)/prcsn)*prcsn)
			if _, ok := keys[ss]; !ok {
				keys[ss] = true
				list = append(list, ss)
			}
		}
		return list
	}
	if len(*pl) > ncheck {
		lpl := len(*pl)
		a := (*pl)[lpl-ncheck : lpl]
		xy := make([]complex128, len(a))
		for i, aa := range a {
			// xy[i] = complex(aa[0], aa[1])
			xy[i] = complex(aa.X, aa.Y)
		}
		uxy := unique(xy)
		// fmt.Println(lpl, len(a), len(uxy))

		if len(uxy) <= xuniq {
			if prnt {
				fmt.Printf("\tparticle has exited domain at prism %d (%6.3f,%6.3f,%6.3f,%6.3e)\n\tcycle found involving %d prisms\n", i, p.X, p.Y, p.Z, p.T, len(uxy))
			}
			return
		}
	}

	// check for well
	switch d.VF[i].(type) {
	case *PollockMethod:
		d.pt = d.VF[i].(*PollockMethod)              // type assertion, analytical solution (no tracking needed)
		if v, ok := d.zw[i]; ok && !cmplx.IsNaN(v) { //!cmplx.IsNaN(d.zw[i]) {
			if prnt {
				fmt.Printf("\tparticle has exited domain at BC prism %d (%6.3f,%6.3f,%6.3f,%6.3e)\n", i, p.X, p.Y, p.Z, p.T)
			}
			return
		}
	case *VectorMethSoln:
		d.pt = d.VF[i].(*VectorMethSoln)             // type assertion, geometrical solution (no tracking needed)
		if v, ok := d.zw[i]; ok && !cmplx.IsNaN(v) { //!cmplx.IsNaN(d.zw[i]) {
			if prnt {
				fmt.Printf("\tparticle has exited domain at BC prism %d (%6.3f,%6.3f,%6.3f,%6.3e)\n", i, p.X, p.Y, p.Z, p.T)
			}
			return
		}
	default:
		wm := d.VF[i].(*WatMethSoln)
		if !cmplx.IsNaN(wm.zwl[0]) { // wells are solved internal to the prism
			for _, zwell := range wm.zwl {
				zl := (complex(p.X, p.Y) - wm.zc) / wm.r // complex local coordinate
				if cmplx.Abs(zl-zwell) < wellTol/real(wm.r) {
					if prnt {
						fmt.Printf("\tparticle has exited by well at prism %d (%6.3f,%6.3f,%6.3f,%6.3e)\n", i, p.X, p.Y, p.Z, p.T)
					}
					return
				}
			}
		}
	}

	// track within prism
	plt := trackToPrismExit(p, d.prsms[i], d.VF[i], d.pt)
	// plt := pt.(*PollockMethod).TestTracktoExit(p, d.prsms[i], d.VF[i]) //  for testing (not concurrent)

	if len(plt) > 2 {
		*pl = append(*pl, plt...) // add tracks
		panic("if never encountered, can remove this")
	}

	// func() {
	// 	n1 := len(*pl) - 1
	// 	for i := 0; i < 4; i++ {
	// 		if plt[0][i] != (*pl)[n1][i] {
	// 			return
	// 		}
	// 	}
	// 	plt = plt[1:]
	// }()

	pids := d.ParticleToPrismIDs(p, i)
	// fmt.Println(i, pids, p.X, p.Y, p.Z, p.T)
	switch len(pids) {
	case 0:
		if prnt {
			fmt.Printf("\tparticle has exited domain at prism %d\n", i)
		}
	case 1:
		if pids[0] == il {
			if i == il {
				panic("todo")
			}
			if prnt {
				fmt.Printf("\ttracking aborted where particle cycle occurred between cells %d-%d (%6.3f,%6.3f,%6.3f,%6.3e)\n", i, il, p.X, p.Y, p.Z, p.T)
			}
		} else if pids[0] == i {
			panic("todo")
		} else {
			d.trackRecurse(p, pl, pids[0], i, prnt)
		}
	default:
		if prnt {
			fmt.Println(" particle likely at edge/vertex")
		}
		// selecting based on closest centroid-to-centroid trajectory
		dsv, isv := math.MaxFloat64, -1
		x0, y0, z0 := d.prsms[i].Centroid()
		pxyz := [3]float64{p.X, p.Y, p.Z}
		for i, pid := range pids {
			x1, y1, z1 := d.prsms[pid].Centroid()
			d, _, _ := vector.PointToLine(pxyz, [3]float64{x0, y0, z0}, [3]float64{x1, y1, z1})
			if d < dsv {
				dsv = d
				isv = i
			}
		}
		d.trackRecurse(p, pl, pids[isv], i, prnt)
	}
}
