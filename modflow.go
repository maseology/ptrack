package ptrack

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"log"
	"strconv"
	"strings"

	"github.com/maseology/mmio"
)

// ReadMODFLOW reads a MODFLOW6 output file
func ReadMODFLOW(fprfx string) Domain {
	grbfp := fmt.Sprintf("%s.disu.grb", fprfx)
	if _, ok := mmio.FileExists(grbfp); !ok {
		grbfp = fmt.Sprintf("%s.dis.grb", fprfx)
		if _, ok := mmio.FileExists(grbfp); !ok {
			panic("ReadMODFLOW: no grb found")
		}
	}

	pset, conn, jaxr := readGRB(grbfp)
	fpcbc := fmt.Sprintf("%s.cbc", fprfx)
	if _, ok := mmio.FileExists(fpcbc); !ok {
		fpcbc = fmt.Sprintf("%s.flx", fprfx)
	}
	pflx, pqw := readCBC(fpcbc, jaxr)
	func() {
		for m, v := range readDependentVariable(fmt.Sprintf("%s.hds", fprfx)) {
			fmt.Printf("  DV: %s\n", m)
			if m == "HEAD" {
				for i, vv := range v {
					if vv < pset[i].Top {
						pset[i].Bn = vv
					}
				}
			}
		}
	}()

	var d Domain
	d.New(pset, conn, pflx, pqw)
	return d
}

func readGRB(fp string) (map[int]*Prism, map[int][]int, map[int]jaxr) {
	buf := mmio.OpenBinary(fp)
	var btyp, bver [50]byte
	if err := binary.Read(buf, binary.LittleEndian, &btyp); err != nil {
		log.Fatalln("Fatal error: readGRB read 001 failed: ", err)
	}
	if err := binary.Read(buf, binary.LittleEndian, &bver); err != nil {
		log.Fatalln("Fatal error: readGRB read 002 failed: ", err)
	}
	ttyp, tver := strings.TrimSpace(string(btyp[:])), strings.TrimSpace(string(bver[:]))
	if tver != "VERSION 1" {
		log.Fatalf("Error:\n GRB %s version not supported: '%s'", fp, tver)
	}

	switch ttyp {
	case "GRID DIS":
		// fmt.Println(ttyp, tver)
		readGRBheader(buf)
		return readGRBgrid(buf)
	case "GRID DISU":
		// fmt.Println(ttyp, tver)
		readGRBheader(buf)
		return readGRBU(buf)
	default:
		log.Fatalf("GRB type '%s' currently not supported", ttyp)
		var p map[int]*Prism
		var c map[int][]int
		var jaxr map[int]jaxr
		return p, c, jaxr
	}
}

func readGRBheader(b *bytes.Reader) {
	// read past *.grb header
	var bntxt, blentxt [50]byte
	if err := binary.Read(b, binary.LittleEndian, &bntxt); err != nil {
		log.Fatalln("Fatal error: readGRBheader read 001 failed: ", err)
	}
	if err := binary.Read(b, binary.LittleEndian, &blentxt); err != nil {
		log.Fatalln("Fatal error: readGRBheader read 002 failed: ", err)
	}
	ntxt, err := strconv.Atoi(strings.TrimSpace(string(bntxt[:])[5:]))
	if err != nil {
		log.Fatalln("Fatal error: readGRBheader read 003 failed: ", err)
	}
	lentxt, err := strconv.Atoi(strings.TrimSpace(string(blentxt[:])[7:]))
	if err != nil {
		log.Fatalln("Fatal error: readGRBheader read 004 failed: ", err)
	}
	for {
		ln := make([]byte, lentxt, lentxt)
		if err := binary.Read(b, binary.LittleEndian, ln); err != nil {
			log.Fatalln("Fatal error: readGRBheader read 005 failed: ", err)
		}
		// fmt.Println(strings.TrimSpace(string(ln[:])))
		ntxt--
		if ntxt <= 0 {
			break
		}
	}
}

func readGRBgrid(buf *bytes.Reader) (map[int]*Prism, map[int][]int, map[int]jaxr) {
	g := grbGridHreader{}
	g.read(buf)

	delr, delc := make(map[int]float64), make(map[int]float64)
	for j := 0; j < int(g.NCOL); j++ {
		delr[j] = mmio.ReadFloat64(buf) // cell width
	}
	for i := 0; i < int(g.NROW); i++ {
		delc[i] = mmio.ReadFloat64(buf) // cell height
		g.YORIGIN += delc[i]            // adjusting origin from lower-left to upper-left
	}

	top, botm := make([]float64, int(g.NROW*g.NCOL)), make([]float64, int(g.NCELLS))
	for i := 0; i < int(g.NROW*g.NCOL); i++ {
		top[i] = mmio.ReadFloat64(buf)
	}
	for i := 0; i < int(g.NCELLS); i++ {
		botm[i] = mmio.ReadFloat64(buf)
	}

	ia, ja := make([]int, int(g.NCELLS)+1), make([]int, int(g.NJA))
	for i := 0; i <= int(g.NCELLS); i++ {
		ia[i] = int(mmio.ReadInt32(buf)) - 1
	}
	for j := 0; j < int(g.NJA); j++ {
		ja[j] = int(mmio.ReadInt32(buf)) - 1
	}

	idomain, icelltype := make([]int, int(g.NCELLS)), make([]int, int(g.NCELLS))
	for i := 0; i < int(g.NCELLS); i++ {
		idomain[i] = int(mmio.ReadInt32(buf))
	}
	for i := 0; i < int(g.NCELLS); i++ {
		icelltype[i] = int(mmio.ReadInt32(buf)) //  specifies how saturated thickness is treated
	}

	if !mmio.ReachedEOF(buf) {
		log.Fatalln("Fatal error: readGRB read 003 failed: have not reached EOF")
	}

	// fmt.Printf("  nl,nr,nc: %v,%v,%v; UL-origin: (%v, %v)\n", g.NLAY, g.NROW, g.NCOL, g.XORIGIN, g.YORIGIN)
	c, cpl, prsms := 0, int(g.NROW*g.NCOL), make(map[int]*Prism)
	for k := 0; k < int(g.NLAY); k++ {
		cl, o := 0, complex(g.XORIGIN, g.YORIGIN) // converted to upper-left (above)
		for i := 0; i < int(g.NROW); i++ {
			dy := -delc[i]
			for j := 0; j < int(g.NCOL); j++ {
				dx := delr[j]
				// p1---p2   y       0---nc
				//  | c |    |       |       clockwise, left-top-right-bottom
				// p0---p3   0---x   nr
				z := []complex128{o + complex(0., dy), o, o + complex(dx, 0.), o + complex(dx, dy)}
				if idomain[c] >= 0 {
					var p Prism
					if k == 0 {
						p.New(z, top[c], botm[c], top[cl], 0., defaultPorosity)
					} else {
						for kk := k - 1; kk >= 0; kk-- {
							c0 := kk*cpl + cl
							if idomain[c0] > 0 {
								p.New(z, botm[c0], botm[c], top[cl], 0., defaultPorosity)
								break
							} else if idomain[c0] == 0 {
								p.New(z, botm[c0], botm[c], botm[c0], 0., defaultPorosity)
								break
							} else if kk == 0 {
								p.New(z, top[cl], botm[c], top[cl], 0., defaultPorosity)
							}
						}
					}
					prsms[c] = &p
					// fmt.Println(c, k, i, j, p.Z, p.Top, p.Bot)
				}
				o += complex(dx, 0.)
				c++
				cl++
			}
			o = complex(g.XORIGIN, imag(o)+dy)
		}
	}

	conn := g.buildTopology() // make(map[int][]int)
	// for i := 0; i < int(g.NCELLS); i++ {
	// 	c1 := make([]int, ia[i+1]-ia[i])
	// 	for j := ia[i]; j < ia[i+1]; j++ {
	// 		c1[j-ia[i]] = ja[j]
	// 	}
	// 	conn[i] = c1
	// }

	// check connections
	jaxrOut, jaxrcnt := make(map[int]jaxr), 0
	for i := 0; i < int(g.NCELLS); i++ {
		i1, c1 := make([]int, ia[i+1]-ia[i]), make([]int, ia[i+1]-ia[i]) // MF6 order (looks to be) above-up-left-right-down-below
		for j := ia[i]; j < ia[i+1]; j++ {
			c1[j-ia[i]] = ja[j]
			i1[j-ia[i]] = j
		}

		connkey := make(map[int]bool) // temporary map for list checking
		for _, v := range conn[i] {
			if v >= 0 {
				connkey[v] = true
			}
		}
		if c1[0] != i {
			log.Fatalf("Fatal error: readGRB cell id check 004 failed:\nCreated: %v\nFound: %v\n", i, c1[0])
		}
		if len(c1)-1 != len(connkey) {
			log.Fatalf("Fatal error: readGRB connectivity check 005 failed, cell %d\n:\nCreated: %v\nFound: %v\n", i, conn[i], c1[1:])
		}
		for _, c := range c1[1:] {
			if !connkey[c] {
				log.Fatalf("Fatal error: readGRB connectivity check 006 failed, cell %d\n:\nCreated: %v\nFound: %v\n", i, conn[i], c1[1:])
			}
		}

		// fmt.Println(conn[i], c1[1:])
		for j, v := range conn[i] {
			for j2, v2 := range c1[1:] {
				if v == v2 {
					jaxrOut[jaxrcnt] = jaxr{
						f: c1[0],
						t: v,
						p: j,
						i: i1[j2+1],
					} // JA to prsm.conn cross-reference
					jaxrcnt++
				}
			}
		}
	}

	if len(jaxrOut) != int(g.NJA-g.NCELLS) {
		log.Fatalf("Fatal error: readGRB connectivity check 007 failed, number of connections created (%d) not equal to NJA (less number of cells)\n", len(jaxrOut))
	}

	// fmt.Println("left-up-right-down-bottom-top")
	// fmt.Println("\nJAXR [from to pos ia]:")
	// for _, v := range jaxrOut {
	// 	fmt.Println(v)
	// }

	return prsms, conn, jaxrOut
}

func readGRBU(buf *bytes.Reader) (map[int]*Prism, map[int][]int, map[int]jaxr) {
	g := grbuHreader{}
	g.read(buf)

	nc := int(g.NODES)
	top, botm := make([]float64, nc), make([]float64, nc)
	if err := binary.Read(buf, binary.LittleEndian, top); err != nil {
		panic(err)
	}
	if err := binary.Read(buf, binary.LittleEndian, botm); err != nil {
		panic(err)
	}

	ia, ja, icelltype := make([]int32, nc+1), make([]int32, int(g.NJA)), make([]int32, nc)
	if err := binary.Read(buf, binary.LittleEndian, ia); err != nil {
		panic(err)
	}
	if err := binary.Read(buf, binary.LittleEndian, ja); err != nil {
		panic(err)
	}
	if err := binary.Read(buf, binary.LittleEndian, icelltype); err != nil {
		panic(err)
	}

	type ddd struct {
		NVERT, NJAVERT int32
	}
	var dddd ddd
	if err := binary.Read(buf, binary.LittleEndian, &dddd); err != nil {
		panic(err)
	}
	fmt.Println(dddd)

	// if !mmio.ReachedEOF(buf) {
	// 	log.Fatalln("Fatal error: readGRB read 003 failed: have not reached EOF")
	// }

	// // fmt.Printf("  nl,nr,nc: %v,%v,%v; UL-origin: (%v, %v)\n", g.NLAY, g.NROW, g.NCOL, g.XORIGIN, g.YORIGIN)
	// c, cpl, prsms := 0, nc, make(map[int]*Prism)
	// for k := 0; k < int(g.NLAY); k++ {
	// 	cl, o := 0, complex(g.XORIGIN, g.YORIGIN) // converted to upper-left (above)
	// 	for i := 0; i < int(g.NROW); i++ {
	// 		dy := -delc[i]
	// 		for j := 0; j < int(g.NCOL); j++ {
	// 			dx := delr[j]
	// 			// p1---p2   y       0---nc
	// 			//  | c |    |       |       clockwise, left-top-right-bottom
	// 			// p0---p3   0---x   nr
	// 			z := []complex128{o + complex(0., dy), o, o + complex(dx, 0.), o + complex(dx, dy)}
	// 			if idomain[c] >= 0 {
	// 				var p Prism
	// 				if k == 0 {
	// 					p.New(z, top[c], botm[c], top[cl], 0., defaultPorosity)
	// 				} else {
	// 					for kk := k - 1; kk >= 0; kk-- {
	// 						c0 := kk*cpl + cl
	// 						if idomain[c0] > 0 {
	// 							p.New(z, botm[c0], botm[c], top[cl], 0., defaultPorosity)
	// 							break
	// 						} else if idomain[c0] == 0 {
	// 							p.New(z, botm[c0], botm[c], botm[c0], 0., defaultPorosity)
	// 							break
	// 						} else if kk == 0 {
	// 							p.New(z, top[cl], botm[c], top[cl], 0., defaultPorosity)
	// 						}
	// 					}
	// 				}
	// 				prsms[c] = &p
	// 				// fmt.Println(c, k, i, j, p.Z, p.Top, p.Bot)
	// 			}
	// 			o += complex(dx, 0.)
	// 			c++
	// 			cl++
	// 		}
	// 		o = complex(g.XORIGIN, imag(o)+dy)
	// 	}
	// }

	// conn := g.buildTopology() // make(map[int][]int)
	// // for i := 0; i < int(g.NCELLS); i++ {
	// // 	c1 := make([]int, ia[i+1]-ia[i])
	// // 	for j := ia[i]; j < ia[i+1]; j++ {
	// // 		c1[j-ia[i]] = ja[j]
	// // 	}
	// // 	conn[i] = c1
	// // }

	// // check connections
	// jaxrOut, jaxrcnt := make(map[int]jaxr), 0
	// for i := 0; i < int(g.NCELLS); i++ {
	// 	i1, c1 := make([]int, ia[i+1]-ia[i]), make([]int, ia[i+1]-ia[i]) // MF6 order (looks to be) above-up-left-right-down-below
	// 	for j := ia[i]; j < ia[i+1]; j++ {
	// 		c1[j-ia[i]] = ja[j]
	// 		i1[j-ia[i]] = j
	// 	}

	// 	connkey := make(map[int]bool) // temporary map for list checking
	// 	for _, v := range conn[i] {
	// 		if v >= 0 {
	// 			connkey[v] = true
	// 		}
	// 	}
	// 	if c1[0] != i {
	// 		log.Fatalf("Fatal error: readGRB cell id check 004 failed:\nCreated: %v\nFound: %v\n", i, c1[0])
	// 	}
	// 	if len(c1)-1 != len(connkey) {
	// 		log.Fatalf("Fatal error: readGRB connectivity check 005 failed, cell %d\n:\nCreated: %v\nFound: %v\n", i, conn[i], c1[1:])
	// 	}
	// 	for _, c := range c1[1:] {
	// 		if !connkey[c] {
	// 			log.Fatalf("Fatal error: readGRB connectivity check 006 failed, cell %d\n:\nCreated: %v\nFound: %v\n", i, conn[i], c1[1:])
	// 		}
	// 	}

	// 	// fmt.Println(conn[i], c1[1:])
	// 	for j, v := range conn[i] {
	// 		for j2, v2 := range c1[1:] {
	// 			if v == v2 {
	// 				jaxrOut[jaxrcnt] = jaxr{
	// 					f: c1[0],
	// 					t: v,
	// 					p: j,
	// 					i: i1[j2+1],
	// 				} // JA to prsm.conn cross-reference
	// 				jaxrcnt++
	// 			}
	// 		}
	// 	}
	// }

	// if len(jaxrOut) != int(g.NJA-g.NCELLS) {
	// 	log.Fatalf("Fatal error: readGRB connectivity check 007 failed, number of connections created (%d) not equal to NJA (less number of cells)\n", len(jaxrOut))
	// }

	// // fmt.Println("left-up-right-down-bottom-top")
	// // fmt.Println("\nJAXR [from to pos ia]:")
	// // for _, v := range jaxrOut {
	// // 	fmt.Println(v)
	// // }

	// return prsms, conn, jaxrOut
	return nil, nil, nil
}

func (g *grbGridHreader) buildTopology() map[int][]int {
	cid, nl, nr, nc := 0, int(g.NLAY), int(g.NROW), int(g.NCOL)
	// fmt.Println(nl, nr, nc)
	tp := make(map[int][]int)
	for k := 0; k < nl; k++ {
		for i := 0; i < nr; i++ {
			for j := 0; j < nc; j++ {
				c1 := []int{-1, -1, -1, -1, -1, -1} // initialize, left-up-right-down-bottom-top

				// left
				if j > 0 {
					c1[0] = cid - 1
				}

				// up
				if i > 0 {
					c1[1] = cid - nc
				}

				// right
				if j < nc-1 {
					c1[2] = cid + 1
				}

				// down
				if i < nr-1 {
					c1[3] = cid + nc
				}

				// bottom/below
				if k < nl-1 {
					c1[4] = cid + nc*nr
				}

				// top/above
				if k > 0 {
					c1[5] = cid - nc*nr
				}

				tp[cid] = c1
				cid++
			}
		}
	}
	return tp
}

type jaxr struct {
	f, t, p, i int
}

type grbGridHreader struct {
	NCELLS, NLAY, NROW, NCOL, NJA int32
	XORIGIN, YORIGIN, ANGROT      float64
}

func (g *grbGridHreader) read(b *bytes.Reader) bool {
	err := binary.Read(b, binary.LittleEndian, g)
	// fmt.Println(*g)
	if err != nil {
		if err == io.EOF {
			return true
		}
		log.Fatalln("Fatal error: grbGridHread failed: ", err)
	}
	return false
}

type grbuHreader struct {
	NODES, NJA               int32
	XORIGIN, YORIGIN, ANGROT float64
}

func (g *grbuHreader) read(b *bytes.Reader) bool {
	err := binary.Read(b, binary.LittleEndian, g)
	// fmt.Println(*g)
	if err != nil {
		if err == io.EOF {
			return true
		}
		log.Fatalln("Fatal error: grbGridHread failed: ", err)
	}
	return false
}

func readCBC(fp string, jaxr map[int]jaxr) (pflx map[int][]float64, pqw map[int]float64) {
	bflx := mmio.OpenBinary(fp)
	dat1D := make(map[string]map[int]float64)
	dat2D := make(map[string]map[int]map[int]float64)

	for {
		h := cbcHreader{}
		if h.cbcHread(bflx) {
			break // EOF
		}

		txt := strings.TrimSpace(string(h.TEXT[:]))
		// fmt.Printf("KSTP %d; KPER %d: %s\n", h.KPER, h.KSTP, txt)
		switch h.IMETH {
		case 1: // Read 1D array of size NDIM1*NDIM2*NDIM3
			m1, n := make(map[int]float64), int(-h.NDIM1*h.NDIM2*h.NDIM3)
			for i := 0; i < n; i++ {
				m1[i] = mmio.ReadFloat64(bflx)
			}
			dat1D[txt] = m1
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
			log.Fatalf("MODFLOW CBC read error: IMETH=%d not supported", h.IMETH)
		}
	}

	// print available outputs
	fmt.Println("  CBC: 2D")
	for i := range dat2D {
		fmt.Printf("      %s\n", i)
	}
	fmt.Println("  CBC: 1D")
	for i := range dat1D {
		fmt.Printf("      %s\n", i)
	}

	pflx = make(map[int][]float64)
	if val, ok := dat1D["FLOW-JA-FACE"]; ok {
		// fmt.Printf("\nFLOW-JA-FACE data (%d):\n", len(val))
		for _, ja := range jaxr { // initialize
			pflx[ja.f] = []float64{0., 0., 0., 0., 0., 0.} // left-up-right-down-bottom-top
		}
		for _, ja := range jaxr {
			// fmt.Printf("from %d to %d flux %v\n", ja.f, ja.t, val[ja.i])
			pflx[ja.f][ja.p] = val[ja.i]
		}
	}

	pqw = make(map[int]float64)
	if val, ok := dat2D["WEL"]; ok {
		// fmt.Println("\nWEL data:")
		for i, v := range val {
			if len(v) > 1 {
				log.Fatalln("MODFLOW CBC read error: WEL given with greater than 1 NDAT")
			}
			// fmt.Println(i, v[0])
			pqw[i] += v[0]
		}
	}
	if val, ok := dat2D["RCH"]; ok {
		// fmt.Println("\nWEL data:")
		for i, v := range val {
			if len(v) > 1 {
				log.Fatalln("MODFLOW CBC read error: RCH given with greater than 1 NDAT")
			}
			// fmt.Println(i, v[0])
			pflx[i][5] = v[0]
		}
	}
	if val, ok := dat2D["CHD"]; ok {
		// fmt.Println("\nCHD data:")
		for i, v := range val {
			if len(v) > 1 {
				log.Fatalln("MODFLOW CBC read error: CHD given with greater than 1 NDAT")
			}
			// fmt.Println(i, v[0])
			pqw[i] += v[0]
		}
	}

	// fmt.Println("\nflux summary [left-up-right-down-bottom-top] well")
	// for i := 0; i < len(pflx); i++ {
	// 	fmt.Println(i, *pflx[i])
	// }

	return
}

type cbcHreader struct {
	KSTP, KPER          int32
	TEXT                [16]byte
	NDIM1               int32 // NCOL (DIS); NCPL (DISV); NODES (DISU); NJA (FLOW-JA-FACE, IMETH=1)
	NDIM2               int32 // NROW (DIS);    1 (DISV);     1 (DISU);   1 (FLOW-JA-FACE, IMETH=1)
	NDIM3               int32 // NLAY (DIS); NLAY (DISV);     1 (DISU);   1 (FLOW-JA-FACE, IMETH=1)
	IMETH               int32
	DELT, PERTIM, TOTIM float64
}

func (h *cbcHreader) cbcHread(b *bytes.Reader) bool {
	err := binary.Read(b, binary.LittleEndian, h)
	if err != nil {
		if err == io.EOF {
			return true
		}
		log.Fatalln("Fatal error: cbcHread failed: ", err)
	}
	return false
}

type cbcAuxReader struct {
	TXT1ID1, TXT2ID1, TXT1ID2, TXT2ID2 [16]byte
	NDAT                               int32
}

func (a *cbcAuxReader) cbcAuxRead(b *bytes.Reader) {
	err := binary.Read(b, binary.LittleEndian, a)
	if err != nil {
		log.Fatalln("Fatal error: cbcAuxRead failed: ", err)
	}
}

func readDependentVariable(fp string) map[string]map[int]float64 {
	bflx, m1 := mmio.OpenBinary(fp), make(map[string]map[int]float64)
	for {
		h := dvarHreader{}
		if h.dvarHread(bflx) {
			break // EOF
		}

		txt := strings.TrimSpace(string(h.TEXT[:]))
		// fmt.Printf("Layer %d; KSTP %d; KPER %d: %s\n", h.ILAY, h.KPER, h.KSTP, txt)
		m2, c := make(map[int]float64), int(h.ILAY-1)*int(h.NROW*h.NCOL)
		for i := 0; i < int(h.NROW); i++ {
			for j := 0; j < int(h.NCOL); j++ {
				m2[c] = mmio.ReadFloat64(bflx)
				c++
			}
		}
		if m1[txt] == nil {
			m1[txt] = make(map[int]float64)
		}
		for i, v := range m2 {
			m1[txt][i] = v
		}
	}
	return m1

	// 	Using br As New BinaryReader(New FileStream(_filepath, FileMode.Open), System.Text.Encoding.Default)
	// 	Dim cnt As Integer = 1
	// 100:            Dim KSTP = br.ReadInt32
	// 	Dim KPER = br.ReadInt32
	// 	Dim PERTIM = br.ReadDouble
	// 	Dim TOTIM = br.ReadDouble
	// 	Dim TEXT = mmIO.HexToString(BitConverter.ToString(br.ReadBytes(16))).Trim
	// 	Dim NCOL = br.ReadInt32 ' NCPL (DISV); NODES (DISU)
	// 	Dim NROW = br.ReadInt32 ' if DISV or DISU, NROW=1
	// 	Dim ILAY = br.ReadInt32 ' if DISU, ILAY=1
	// 	Dim dic1 As New Dictionary(Of Integer, Double), cnt2 = 0
	// 	For i = 1 To NROW
	// 		For j = 1 To NCOL
	// 			For k = 1 To ILAY
	// 				dic1.Add(cnt2, br.ReadDouble)
	// 				cnt2 += 1
	// 			Next
	// 		Next
	// 	Next
	// 	_t.Add(cnt, dic1)
	// 	_tnam.Add(cnt, String.Format("{0}_SP{1:00000}_TS{2:00000}", TEXT, KPER, KSTP))
	// 	cnt += 1
	// 	If br.PeekChar <> -1 Then GoTo 100
	// End Using
}

type dvarHreader struct {
	KSTP, KPER    int32
	PERTIM, TOTIM float64
	TEXT          [16]byte
	NCOL          int32 // NCOL (DIS); NCPL (DISV); NODES (DISU)
	NROW          int32 // NROW (DIS);    1 (DISV);     1 (DISU)
	ILAY          int32 // NLAY (DIS); NLAY (DISV);     1 (DISU)
}

func (h *dvarHreader) dvarHread(b *bytes.Reader) bool {
	err := binary.Read(b, binary.LittleEndian, h)
	if err != nil {
		if err == io.EOF {
			return true
		}
		log.Fatalln("Fatal error: dvarHread failed: ", err)
	}
	return false
}
