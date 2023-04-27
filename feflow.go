package ptrack

import (
	"fmt"
	"log"

	"github.com/maseology/goHydro/mesh"
)

// currently reads only 1 state.
func ReadHSTRAT(hstratFP string) (Domain, *mesh.Slice) {

	h, err := mesh.ReadHSTRAT(hstratFP, true)
	if err != nil {
		panic(err)
	}

	// get prisms
	fmt.Printf(" Model structure read: %s nodes, %s elements, %d layers\n", big(h.Nn), big(h.Ne), h.Nly)
	pset := func() map[int]*Prism {
		prsms := make(map[int]*Prism, h.Ne)
		//  p1--p2
		//   | /          clockwise, left-top-right-bottom
		//   p0
		for eid, nds := range h.Exr {
			nh2 := len(nds) / 2
			fnh := float64(len(nds))
			z := make([]complex128, nh2)
			top, bot := 0., 0.
			for i := 0; i < nh2; i++ {
				xyz := h.Nxyz[nds[i]]
				z[nh2-i-1] = complex(xyz[0], xyz[1]) // reverse order
				// z[i] = complex(xyz[0], xyz[1])
				top += xyz[2] // averaging top elevation
			}
			// checkConsistentOrder(z)
			for i := nh2; i < 2*nh2; i++ {
				xyz := h.Nxyz[nds[i]]
				bot += xyz[2] // averaging bottom elevation
			}
			prsms[eid] = &Prism{
				Z:   z,
				Top: top / fnh,
				Bot: bot / fnh,
				Por: float64(h.Hgeo[eid].N),
				Bn:  top / fnh,
				Tn:  0.,
			}
			prsms[eid].computeArea()
		}
		return prsms
	}()

	// get saturation
	if len(h.Nh) == h.Nn {
		for eid, nids := range h.Exr {
			sh := 0.
			for _, nid := range nids {
				sh += h.Nh[nid]
			}
			sh /= float64(len(nids))
			if sh < pset[eid].Top { // Bn defaulted to top
				if sh < pset[eid].Bot {
					pset[eid].Bn = pset[eid].Bot // dry cell
				} else {
					pset[eid].Bn = sh
				}
			}
		}
	} else {
		log.Fatalln("ReadHSTRAT todo")
	}

	// get flux
	pflx := make(map[int][]float64, h.Ne)
	if len(h.Vxyz) == h.Nn {
		for eid, nids := range h.Exr {
			q := make([]float64, 3)
			for _, nid := range nids {
				for j := 0; j < 3; j++ {
					q[j] += h.Vxyz[nid][j]
				}
			}
			f := float64(h.Hgeo[eid].N) / float64(len(nids))
			for j := 0; j < 3; j++ {
				q[j] *= f
			}
			pflx[eid] = q
		}
	} else {
		log.Fatalln("ReadHSTRAT todo")
	}

	fmt.Printf("  results collected at end of simulation: %s vectors and %s scalars collected\n", big(len(pflx)), big(len(h.Nh)))

	var d Domain
	d.New(pset, h.BuildElementalConnectivity(false), pflx, nil)
	d.Nly = h.Nly
	d.Minthick = h.MinThick
	return d, h.TopSlice()
}

// func checkConsistentOrder(zs []complex128) {
// 	s, nfaces := 0., len(zs)
// 	for j := range zs {
// 		jj := (j + 1) % nfaces
// 		s += (real(zs[jj]) - real(zs[j])) * (imag(zs[jj]) + imag(zs[j]))
// 	}
// 	if s >= 0. {
// 		fmt.Println("cw")
// 	} else {
// 		fmt.Println("ccw")
// 	}
// 	// Dim db1 As Double = 0.0
// 	// For i = 0 To _v.Count - 2
// 	// 	db1 += (_v(i + 1).X - _v(i).X) * (_v(i + 1).Y + _v(i).Y)
// 	// Next
// 	// db1 += (_v(0).X - _v.Last.X) * (_v(0).Y + _v.Last.Y)
// 	// Return db1 >= 0.0
// }
