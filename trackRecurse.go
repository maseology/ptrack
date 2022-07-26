package ptrack

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/maseology/mmaths/vector"
)

func (d *Domain) trackRecurse(p *Particle, pl *[][]float64, i, il int) {
	pt := d.pt // var pt ParticleTracker
	*pl = append(*pl, p.State())

	// check for loops
	const ncheck, xuniq, prcsn = 100, 10, .01
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
			xy[i] = complex(aa[0], aa[1])
		}
		uxy := unique(xy)
		// fmt.Println(lpl, len(a), len(uxy))

		if len(uxy) <= xuniq {
			fmt.Printf("\tparticle has exited domain at prism %d (%6.3f,%6.3f,%6.3f,%6.3e)\n\tcycle found involving %d prisms\n", i, p.X, p.Y, p.Z, p.T, len(uxy))
			return
		}
	}

	// check for well
	switch d.VF[i].(type) {
	case *PollockMethod:
		pt = d.VF[i].(*PollockMethod)                // type assertion, analytical solution (no tracking needed)
		if v, ok := d.zw[i]; ok && !cmplx.IsNaN(v) { //!cmplx.IsNaN(d.zw[i]) {
			fmt.Printf("\tparticle has exited domain at BC prism %d (%6.3f,%6.3f,%6.3f,%6.3e)\n", i, p.X, p.Y, p.Z, p.T)
			return
		}
	case *VectorMethSoln:
		pt = d.VF[i].(*VectorMethSoln)               // type assertion, analytical solution (no tracking needed)
		if v, ok := d.zw[i]; ok && !cmplx.IsNaN(v) { //!cmplx.IsNaN(d.zw[i]) {
			fmt.Printf("\tparticle has exited domain at BC prism %d (%6.3f,%6.3f,%6.3f,%6.3e)\n", i, p.X, p.Y, p.Z, p.T)
			return
		}
	default:
		wm := d.VF[i].(*WatMethSoln)
		if !cmplx.IsNaN(wm.zwl[0]) { // wells are solved internal to the prism
			for _, zwell := range wm.zwl {
				zl := (complex(p.X, p.Y) - wm.zc) / wm.r // complex local coordinate
				if cmplx.Abs(zl-zwell) < wellTol/real(wm.r) {
					fmt.Printf("\tparticle has exited by well at prism %d (%6.3f,%6.3f,%6.3f,%6.3e)\n", i, p.X, p.Y, p.Z, p.T)
					return
				}
			}
		}
	}

	// track within prism
	// plt := trackToPrismExit(p, d.prsms[i], d.VF[i], pt)
	plt := pt.(*PollockMethod).TestTracktoExit(p, d.prsms[i], d.VF[i]) //  for testing (not concurrent)

	*pl = append(*pl, plt...) // add tracks

	pids := d.ParticleToPrismIDs(p, i)
	// fmt.Println(i, pids, p.X, p.Y, p.Z, p.T)
	switch len(pids) {
	case 0:
		fmt.Printf("\tparticle has exited domain at prism %d\n", i)
	case 1:
		if pids[0] == il {
			if i == il {
				panic("todo")
			}
			fmt.Printf("\ttracking aborted where particle cycle occurred between cells %d-%d (%6.3f,%6.3f,%6.3f,%6.3e)\n", i, il, p.X, p.Y, p.Z, p.T)
		} else if pids[0] == i {
			panic("todo")
		} else {
			d.trackRecurse(p, pl, pids[0], i)
		}
	default: // particle likely at edge/vertex
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
		d.trackRecurse(p, pl, pids[isv], i)
	}

	// switch len(pids) {
	// case 0:
	// 	fmt.Printf("\tparticle has exited domain at prism %d\n", i)
	// case 1:
	// 	if pids[0] == i {
	// 		fmt.Printf("\ttracking aborted where particle exited cell %d\n", i)
	// 		*pl = append(*pl, p.State()) // adding final point
	// 	} else if pids[0] == il {
	// 		if i == il {
	// 			fmt.Printf("\ttracking aborted where particle exited cell %d\n", i)
	// 			*pl = append(*pl, p.State()) // adding final point
	// 		} else {
	// 			// q0, q1 := d.flx[il], d.flx[i]
	// 			// fmt.Println(q0, q1)
	// 			// fmt.Println(pos) //left-up-right-down-bottom-top

	// 			fmt.Printf("\ttracking aborted where particle cycle occurred between cells %d-%d (%6.3f,%6.3f,%6.3f,%6.3e)\n", i, il, p.X, p.Y, p.Z, p.T)
	// 		}
	// 	} else {
	// 		fmt.Printf("\tparticle path from prism %d entering %d (%6.3f,%6.3f,%6.3f,%6.3e)\n", i, pids[0], p.X, p.Y, p.Z, p.T)
	// 		d.trackRecurse(p, pl, pids[0], i)
	// 	}
	// case 2:
	// 	if (pids[0] == i && pids[1] == il) || (pids[0] == il && pids[1] == i) {
	// 		// q0, q1 := d.flx[il], d.flx[i]
	// 		// fmt.Println(q0, q1)
	// 		// fmt.Println(pos) //left-up-right-down-bottom-top

	// 		fmt.Printf("\ttracking aborted where particle cycle occurred between cells %d-%d.\n", i, il)

	// 	} else if pids[0] == i {
	// 		fmt.Printf("\tparticle path from prism %d entering %d (%6.3f,%6.3f,%6.3f,%6.3e)\n", i, pids[1], p.X, p.Y, p.Z, p.T)
	// 		d.trackRecurse(p, pl, pids[1], i)
	// 	} else if pids[1] == i {
	// 		fmt.Printf("\tparticle path from prism %d entering %d (%6.3f,%6.3f,%6.3f,%6.3e)\n", i, pids[0], p.X, p.Y, p.Z, p.T)
	// 		d.trackRecurse(p, pl, pids[0], i)
	// 	} else {
	// 		fmt.Printf("\tparticle path has diverged at %d, going with pids=%v (%6.3f,%6.3f,%6.3f,%6.3e)\n", i, pids, p.X, p.Y, p.Z, p.T)
	// 		d.trackRecurse(p, pl, pids[0], i)
	// 	}
	// default:
	// 	for _, pid := range pids {
	// 		if pid != i {
	// 			fmt.Printf("\tparticle path has diverged at %d, going with pids=%v; (randomly) attempting prism %d\n", i, pids, pid)
	// 			d.trackRecurse(p, pl, pid, i)
	// 			break
	// 		}
	// 	}
	// }
}
