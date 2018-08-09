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

// ReadMODFLOW reads a MODFLOW output file and builds its velocity field
func ReadMODFLOW(fprfx string) {
	readGRB(fmt.Sprintf("%s.dis.grb", fprfx))
	// readCBC(fmt.Sprintf("%s.cbc", fprfx))
}

func readGRB(fp string) PrismSet {
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
		fmt.Println(ttyp, tver)
		readGRBheader(buf)
		return readGridGRB(buf)
	default:
		log.Fatalf("GRB type '%s' currently not supported", ttyp)
		var p PrismSet
		return p
	}
}

func readGRBheader(b *bytes.Reader) {
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

func readGridGRB(buf *bytes.Reader) PrismSet {
	g := grbGridHreader{}
	g.grbGridHread(buf)

	delr, delc := make(map[int]float64), make(map[int]float64)
	for j := 0; j < int(g.NCOL); j++ {
		delr[j] = mmio.ReadFloat64(buf) // cell width
	}
	for i := 0; i < int(g.NROW); i++ {
		delc[i] = mmio.ReadFloat64(buf) // cell height
		g.YOFFSET += delc[i]            // adjusting origin from lower-left to upper-left
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

	c, cpl, prsms := 0, int(g.NROW*g.NCOL), make(map[int]Prism)
	for k := 0; k < int(g.NLAY); k++ {
		cl, o := 0, complex(g.XOFFSET, g.YOFFSET) // origin upper-left (converted above)
		for i := 0; i < int(g.NROW); i++ {
			for j := 0; j < int(g.NCOL); j++ {
				dx, dy := delr[j], -delc[i]
				z := []complex128{o, o + complex(dx, 0.), o + complex(dx, dy), o + complex(0., dy)}
				if idomain[c] >= 0 {
					var p Prism
					if k == 0 {
						p.New(z, top[c], botm[c], top[cl], 0., 0.3)
					} else {
						for kk := k - 1; kk >= 0; kk-- {
							c0 := kk*cpl + cl
							if idomain[c0] > 0 {
								p.New(z, botm[c0], botm[c], top[cl], 0., 0.3)
								break
							} else if idomain[c0] == 0 {
								p.New(z, botm[c0], botm[c], botm[c0], 0., 0.3)
								break
							} else if kk == 0 {
								p.New(z, top[cl], botm[c], top[cl], 0., 0.3)
							}
						}
					}
					prsms[c] = p
				}
				o += complex(delr[j], 0.)
				c++
				cl++
			}
			o = complex(g.XOFFSET, g.YOFFSET-delc[i])
		}
	}

	conn := make(map[int][]int)
	for i := 0; i < int(g.NCELLS); i++ {
		c1 := make([]int, ia[i+1]-ia[i])
		for j := ia[i]; j < ia[i+1]; j++ {
			c1[j-ia[i]] = ja[j]
		}
		conn[i] = c1
	}

	return PrismSet{
		P:    prsms,
		Conn: conn,
	}
}

type grbGridHreader struct {
	NCELLS, NLAY, NROW, NCOL, NJA int32
	XOFFSET, YOFFSET, ANGROT      float64
}

func (g *grbGridHreader) grbGridHread(b *bytes.Reader) bool {
	err := binary.Read(b, binary.LittleEndian, g)
	if err != nil {
		if err == io.EOF {
			return true
		}
		log.Fatalln("Fatal error: grbGridHread failed: ", err)
	}
	return false
}

func readCBC(fp string) {
	bflx := mmio.OpenBinary(fp)
	dat1D := make(map[string]map[int]float64)
	dat2D := make(map[string]map[int]map[int]float64)

	for {
		h := cbcHreader{}
		if h.cbcHread(bflx) {
			break // EOF
		}

		txt := string(h.TEXT[:])
		fmt.Printf("KSTP %d; KPER %d: %s\n", h.KPER, h.KSTP, txt)
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
			id1, id2, d2D := make(map[int]int32), make(map[int]int32), make(map[int]map[int]float64)
			for i := 0; i < int(nlist); i++ {
				if err := binary.Read(bflx, binary.LittleEndian, id1[i]); err != nil {
					log.Fatalln("Fatal error: ID1 read failed: ", err)
				}
				if err := binary.Read(bflx, binary.LittleEndian, id2[i]); err != nil {
					log.Fatalln("Fatal error: ID2 read failed: ", err)
				}
				m1 := make(map[int]float64)
				for i := 0; i < int(a.NDAT); i++ {
					m1[i] = mmio.ReadFloat64(bflx)
				}
				d2D[i] = m1
			}
			dat2D[txt] = d2D
		default:
			log.Fatalf("MODFLOW CBC read error: IMETH=%d not supported", h.IMETH)
		}
	}
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
