package ptrack

import (
	"encoding/binary"
	"fmt"
	"log"
	"strings"

	"github.com/maseology/goHydro/grid"
	"github.com/maseology/mmio"
)

// currently needs an accompanying .gdef; also not reading wells
func ReadFLX(hstratFP, flxFP string) (Domain, *grid.Definition) {

	// get geometry
	gd, err := grid.ReadGDEF(mmio.RemoveExtension(hstratFP)+".gdef", true)
	if err != nil {
		gd, err = grid.ReadGDEF(mmio.RemoveExtension(flxFP)+".gdef", true)
		if err != nil {
			log.Fatalf(" ReadFLX gdef read error: %v", err)
		}
	}
	fmt.Printf(" GDEF read: %s cells (%d rows, %d columns), %s actives\n", big(gd.Ncells()), gd.Nrow, gd.Ncol, big(gd.Nact))
	hstrat, err := grid.ReadHSTRAT(hstratFP, true)
	if err != nil {
		log.Fatalf(" ReadFLX hstrat read error: %v", err)
	}
	if len(hstrat.Cells) != gd.Nact*hstrat.Nlay {
		if len(hstrat.Cells) != gd.Nrow*gd.Ncol*hstrat.Nlay {
			log.Fatalf(" ReadFLX hstrat count error, nact = %d, nlay = %d, nhs = %d", gd.Nact, hstrat.Nlay, len(hstrat.Cells))
		}
		print("")
	}
	fmt.Printf(" HSTRAT read: %s cells (%d layers), %s elements per layer\n\n", big(len(hstrat.Cells)), hstrat.Nlay, big(len(hstrat.Cells)/hstrat.Nlay))

	// read flux
	qR, qF, qL, nlay := readMF2005(flxFP, true)
	// if nlay != hstrat.Nlay {
	// 	log.Fatalf(" ReadFLX hstrat layer error: %v", err)
	// }
	fmt.Printf("\n MODFLOW CBC read: flux across %s faces,  %d layers\n", big(len(qR)+len(qF)+len(qL)), nlay)

	// build connectivity
	fmt.Print(" building connectivity")
	conn := func() map[int][]int {
		c, nc := make(map[int][]int, gd.Nact*nlay), gd.Ncells()
		ret1 := func(x, b int) int {
			if x < 0 {
				return x
			}
			return x + b
		}
		for ly := 0; ly < nlay; ly++ {
			print(".")
			lync := ly * nc
			for _, cid := range gd.Sactives {
				lcid := cid + lync
				b := gd.Buffer(cid, true, true)                                                                 // up-left-right-down
				c[lcid] = []int{ret1(b[1], lync), ret1(b[0], lync), ret1(b[2], lync), ret1(b[3], lync), -1, -1} // left-up-right-down-bottom-top
				if ly > 0 {
					c[lcid][5] = cid + lync - nc // top
				}
				if ly < nlay-1 {
					c[lcid][4] = cid + lync + nc // bottom
				}
			}
		}
		return c
	}()

	fmt.Println("\n building prisms..")
	pset := func() map[int]*Prism {
		prsms := make(map[int]*Prism, gd.Nact*nlay)
		// p1---p2   y       0---nc
		//  | c |    |       |       clockwise, left-top-right-bottom
		// p0---p3   0---x   nr
		nc := gd.Ncells()
		for _, cid := range gd.Sactives {
			// cw2 := gd.Cwidth / 2 // uniform cells
			// for _, cid := range gd.Sactives {
			// 	c := gd.Coord[cid]
			// 	p0 := complex(c.X-cw2, c.Y-cw2)
			// 	p1 := complex(c.X-cw2, c.Y+cw2)
			// 	p2 := complex(c.X+cw2, c.Y+cw2)
			// 	p3 := complex(c.X+cw2, c.Y-cw2)
			// 	zs := []complex128{p0, p1, p2, p3}

			zs := func() []complex128 {
				p := gd.CellPerimeter(cid)
				return []complex128{complex(p[0][0], p[0][1]), complex(p[1][0], p[1][1]), complex(p[2][0], p[2][1]), complex(p[3][0], p[3][1])}
			}()
			for k := 0; k < nlay; k++ {
				lcid := cid + k*nc
				t := float64(hstrat.Cells[lcid].Top)
				b := float64(hstrat.Cells[lcid].Bottom)
				// bn := t // fully saturated
				bn := float64(hstrat.Cells[lcid].H0)
				if bn < b {
					bn = b
				} else if bn > t {
					bn = t
				}

				prsms[lcid] = &Prism{
					Z:   zs,
					Top: t,
					Bot: b,
					Por: defaultPorosity,
					Bn:  bn,
					Tn:  0.,
				}
				prsms[lcid].computeArea()
			}
		}
		return prsms
	}()

	// get flux
	fmt.Println(" collecting flux..")
	pflx := func() map[int][]float64 {
		q := make(map[int][]float64, len(conn))
		for cid, cc := range conn {
			a := make([]float64, 6) // left-up-right-down-bottom-top
			// a[0] = -qR[cc[0]]
			// a[1] = -qF[cc[1]]
			// a[2] = qR[cid]
			// a[3] = qF[cid]
			// a[4] = qL[cid]
			// a[5] = -qL[cc[5]]
			a[0] = qR[cc[0]]
			a[1] = qF[cc[1]]
			a[2] = -qR[cid]
			a[3] = -qF[cid]
			a[4] = qL[cid]
			a[5] = -qL[cc[5]]
			q[cid] = a
		}
		return q
	}()

	var d Domain
	d.New(pset, conn, pflx, nil)
	d.Minthick = hstrat.MinThick
	return d, gd
}

func readMF2005(fp string, prnt bool) (qRight, qFront, qLower map[int]float64, nlay int) {
	bflx := mmio.OpenBinary(fp)
	dat1D := make(map[string]map[int]float64)
	dat2D := make(map[string]map[int]map[int]float64)

	for {
		KSTP, ok := mmio.ReadInt32check(bflx)
		if !ok {
			break
		}
		KPER := mmio.ReadInt32(bflx)
		PNAME := mmio.ReadBytes(bflx, 16)
		NC := mmio.ReadInt32(bflx)
		NR := mmio.ReadInt32(bflx)
		NL := mmio.ReadInt32(bflx)
		nlay = func(x32 int32) int {
			x := int(x32)
			if x < 0 {
				return -x
			}
			return x
		}(NL)

		ICODE := int32(0)
		if NL < 0 {
			ICODE = mmio.ReadInt32(bflx)
			DELTS := mmio.ReadFloat32(bflx)
			PERTIMS := mmio.ReadFloat32(bflx)
			TOTIMS := mmio.ReadFloat32(bflx)
			_ = DELTS
			_ = PERTIMS
			_ = TOTIMS
			NL = -NL
		}
		nc2 := int(NC) * int(NR)
		nc3 := nc2 * nlay

		txt := strings.TrimSpace(string(PNAME[:]))
		if prnt {
			fmt.Printf("%15s: ICODE %d; KSTP %d; KPER %d; NC %d; NR %d; NL %d\n", txt, ICODE, KPER, KSTP, NC, NR, NL)
		}

		switch ICODE {
		case 0, 1: // Read 1D array of size NDIM1*NDIM2*NDIM3
			m1 := make(map[int]float64, nc3)
			for i := 0; i < nc3; i++ {
				m1[i] = float64(mmio.ReadFloat32(bflx))
			}
			dat1D[txt] = m1
		case 2: // by cell index
			NLST := mmio.ReadInt32(bflx)
			m1 := make(map[int]float64, NLST)
			for i := 0; i < int(NLST); i++ {
				ICELL := mmio.ReadInt32(bflx)
				VAL1 := mmio.ReadFloat32(bflx)
				m1[int(ICELL)] = float64(VAL1)
			}
			dat1D[txt] = m1
		case 3: // 2D scaler to layer
			// get layer receiving scaler
			i1 := make(map[int]float64, nc2)
			for i := 0; i < nc2; i++ {
				i1[i] = float64(mmio.ReadInt32(bflx))
			}
			dat1D[txt+"-layer"] = i1

			// get values
			m1 := make(map[int]float64, nc2)
			for i := 0; i < nc2; i++ {
				m1[i] = float64(mmio.ReadFloat32(bflx))
			}
			dat1D[txt] = m1

		case 5: // by boundary index
			NAUX := mmio.ReadInt32(bflx) - 1
			AUXTEXT := make([]string, NAUX)
			for i := 0; i < int(NAUX); i++ {
				ba := mmio.ReadBytes(bflx, 16)
				AUXTEXT[i] = strings.TrimSpace(string(ba[:]))
				fmt.Println("               - AUXTXT: " + AUXTEXT[i])
			}

			NLST := mmio.ReadInt32(bflx)
			if NAUX == 0 {
				m1 := make(map[int]float64, NLST)
				for i := 0; i < int(NLST); i++ {
					ICELL := mmio.ReadInt32(bflx)
					VAL1 := mmio.ReadFloat32(bflx)
					m1[int(ICELL)] = float64(VAL1)
				}
				dat1D[txt] = m1
			} else {
				for n := 0; n < int(NAUX); n++ {
					xface := make([]float32, NLST) // ISTRM in SFR (see line 4220 in gwf2sfr7.f)
					m1 := make(map[int]float64, NLST)
					for i := 0; i < int(NLST); i++ {
						ICELL := mmio.ReadInt32(bflx)
						VAL1 := mmio.ReadFloat32(bflx)
						m1[int(ICELL)] = float64(VAL1)
						xface[i] = mmio.ReadFloat32(bflx)
					}
					dat1D[AUXTEXT[n]] = m1
				}
			}

		case 6: // Read text identifiers, auxiliary text labels, and list of information.
			a := cbcAuxReader{}
			a.cbcAuxRead(bflx)
			auxtext := make([]string, int(a.NDAT))
			for i := 0; i < int(a.NDAT)-1; i++ {
				var b1 [16]byte
				if err := binary.Read(bflx, binary.LittleEndian, &b1); err != nil {
					log.Fatalln("Fatal error: AUXTEXT read failed: ", err)
				}
				auxtext[i] = string(b1[:])
			}
			var nlist int32
			if err := binary.Read(bflx, binary.LittleEndian, &nlist); err != nil {
				log.Fatalln("Fatal error: NLIST read failed: ", err)
			}
			d2D := make(map[int]map[int]float64)
			for i := 0; i < int(nlist); i++ {
				var id1, id2 int32
				if err := binary.Read(bflx, binary.LittleEndian, &id1); err != nil {
					log.Fatalln("Fatal error: ID1 read failed: ", err)
				}
				if err := binary.Read(bflx, binary.LittleEndian, &id2); err != nil {
					log.Fatalln("Fatal error: ID2 read failed: ", err)
				}
				m1 := make(map[int]float64)
				for j := 0; j < int(a.NDAT); j++ {
					m1[j] = mmio.ReadFloat64(bflx)
				}
				d2D[int(id1)-1] = m1
			}
			dat2D[txt] = d2D
		default:
			log.Fatalf("MODFLOW CBC read error: IMETH=%d not supported", ICODE)
		}
	}

	// print available outputs
	if prnt {
		fmt.Println("  CBC: 2D")
		for i := range dat2D {
			fmt.Printf("      %s\n", i)
		}
		fmt.Println("  CBC: 1D")
		for i := range dat1D {
			fmt.Printf("      %s\n", i)
		}
	}

	return dat1D["FLOW RIGHT FACE"], dat1D["FLOW FRONT FACE"], dat1D["FLOW LOWER FACE"], nlay
}
