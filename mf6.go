package ptrack

import (
	gomf6 "github.com/maseology/goMF6"
)

// only works for an "unstructured" quadralinear grid
func ReadMF6(mf6Prfx string) Domain {
	// read flux
	mf6 := gomf6.ReadMF6(mf6Prfx)

	// geometry
	nc := len(mf6.Prsms)
	pset := make(map[int]*Prism, nc)
	for nid, mfprsm := range mf6.Prsms {
		getbn := func(bn float64) float64 {
			if bn < mfprsm.Bot {
				bn = mfprsm.Bot
			} else if bn > mfprsm.Top {
				bn = mfprsm.Top
			}
			return bn
		}
		pset[nid] = &Prism{
			Z:    mfprsm.Z,
			Top:  mfprsm.Top,
			Bot:  mfprsm.Bot,
			Area: mfprsm.A,
			Por:  mfprsm.Por,
			Bn:   getbn(mfprsm.H0), // saturated thickness at time step tn
			Tn:   0.,               // initial time step (both bn and tn will adjust in transient cases)
		}
	}

	// flux/transfers
	conn := make(map[int][]int, nc)
	pflx := make(map[int][]float64, nc)
	extentsToDirection := func(pFrom, pTo *Prism) int {
		cx, cy, cz := pTo.Centroid()
		if pFrom.ContainsXY(cx, cy) {
			if cz < pFrom.Bot {
				return 4 // bottom
			} else if cz > pFrom.Top {
				return 5 // top
			} else {
				panic("ReadMF6cbc.extentsToDirection: point in prism 001")
			}
		}

		if cz < pFrom.Bot || cz > pFrom.Top {
			panic("ReadMF6cbc.extentsToDirection: non-cardinality")
		}

		yn, yx, xn, xx := pFrom.getExtentsXY()
		if cx >= xn && cx <= xx {
			if cy > yx {
				return 1 // up
			} else if cy < yn {
				return 3 // down
			} else {
				panic("ReadMF6cbc.extentsToDirection: point in prism 002")
			}
		} else if cy >= yn && cy <= yx {
			if cx > xx {
				return 2 // right
			} else if cx < xn {
				return 0 // left
			} else {
				panic("ReadMF6cbc.extentsToDirection: point in prism 003")
			}
		} else {
			panic("ReadMF6cbc.extentsToDirection: unknown")
		}
	}

	for nid, mfprsm := range mf6.Prsms {
		conn[nid] = []int{-1, -1, -1, -1, -1, -1} // convert unstructured flux to left-up-right-down-bottom-top
		pflx[nid] = make([]float64, 6)
		for i, cnid := range mfprsm.Conn {
			pos := extentsToDirection(pset[nid], pset[cnid])
			conn[nid][pos] = cnid
			pflx[nid][pos] += mfprsm.Q[i+1]
		}

		// // add recharge as top flux (it appears this is not the case with MODPATH)
		// nq := len(mfprsm.Q) - len(mfprsm.Conn) - 1
		// if nq == 1 {
		// 	pflx[nid][5] += mfprsm.Q[len(mfprsm.Q)-1] // top
		// } else if nq != 0 {
		// 	panic("ReadMF6cbc mf6.Prsms flux read")
		// }
	}

	var d Domain
	d.New(pset, conn, pflx, mf6.Qw)
	// d.Minthick = hstrat.MinThick
	return d
}
