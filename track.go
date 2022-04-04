package ptrack

import (
	"fmt"
	"log"
	"math/cmplx"
)

// Track a collection of particles through the domain
func (d *Domain) TrackParticles(p Particles) ([][][]float64, int) {
	o, c := make([][][]float64, len(p)), 0
	for k, pp := range p {
		o[k] = d.trackParticle(&pp)
		c += len(o[k])
	}
	return o, c
}

func (d *Domain) trackParticle(p *Particle) [][]float64 {
	var pathline [][]float64
	pid := d.findStartingPrism(p)

	fmt.Printf("  >>> particle %d start point (x,y,z): %6.3f %6.3f %6.3f released in prism %d\n", p.I, p.X, p.Y, p.Z, pid)
	d.trackRecurse(p, &pathline, pid, -1)

	// pathline = pathline[:len(pathline)-1]
	plast := pathline[len(pathline)-1]
	fmt.Printf("\tparticle exit point  (x,y,z,t): %6.3f %6.3f %6.3f %6.3es\n", plast[0], plast[1], plast[2], plast[3])
	return pathline
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

func (d *Domain) trackRecurse(p *Particle, pl *[][]float64, i, il int) {
	pt := d.pt // var pt ParticleTracker
	*pl = append(*pl, p.State())

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
	plt := trackToPrismExit(p, d.prsms[i], d.VF[i], pt)
	*pl = append(*pl, plt...) // add tracks

	pids := d.ParticleToPrismIDs(p, i)
	fmt.Println(pids)
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
	default:
		panic("todo")
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
