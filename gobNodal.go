package ptrack

import (
	"encoding/gob"
	"fmt"
	"os"

	"github.com/maseology/goHydro/grid"
	"github.com/maseology/mmio"
)

type Vert struct {
	X, Y, Z, T                                               float64
	VertID, PathID, PrsmID, CellID, Layer, Order, USid, DSid int
}

func (d *Domain) ExportNetworkGob(fp string, gd *grid.Definition, pl [][][]float64, xr []int) error {

	// build topology
	var nds []Vert
	nc, k := gd.Ncells(), 0
	for i, pln := range pl {
		pnds := make([]Vert, len(pln))
		// pidlast := -1
		for j, v := range pln {
			pid, cid, ly := func() (int, int, int) {
				cid := gd.PointToCellID(v[0], v[1])
				ly := 0
				for {
					if p, ok := d.prsms[cid+ly*nc]; ok {
						if v[2] <= p.Top && v[2] >= p.Bot {
							return cid + ly*nc, cid, ly + 1
						}
					} else {
						break
					}
					ly++
				}
				return -1, cid, ly
			}()
			// pid, ly := func() (int, int) {
			// 	p := &Particle{X: v[0], Y: v[1], Z: v[2]}
			// 	cids := d.ParticleToPrismIDs(p, pidlast)
			// 	if len(cids) != 1 {
			// 		print("")
			// 	}
			// 	ly := -1
			// 	if _, ok := d.prsms[cids[0]]; !ok {
			// 		print("")
			// 	}
			// 	if !d.prsms[cids[0]].ContainsXYZ(v[0], v[1], v[2]) {
			// 		print("")
			// 	}
			// 	pidlast = cids[0]
			// 	return cids[0], ly
			// }()
			pnds[j] = Vert{
				X:      v[0],
				Y:      v[1],
				Z:      v[2],
				T:      v[3],
				VertID: k,
				PathID: xr[i],
				PrsmID: pid, // the prism/cell it's in
				CellID: cid, // the 2d grid cell it's in
				Layer:  ly,
				Order:  j,
				USid:   -1,
				DSid:   -1,
			}
			k++
		}

		// build topology
		switch len(pnds) {
		case 0:
			// do nothing
		case 1:
			pnds[0].USid = -1
			pnds[0].DSid = -1
		default:
			for j, p := range pnds {
				if j == 0 {
					p.USid = -1
					p.DSid = pnds[j+1].VertID
				} else if j == len(pnds)-1 {
					p.USid = pnds[j-1].VertID
					p.DSid = -1
				} else {
					p.USid = pnds[j-1].VertID
					p.DSid = pnds[j+1].VertID
				}
				pnds[j] = p
			}
			nds = append(nds, pnds...) // ignoring 1-point vertices
		}
		// nds = append(nds, pnds...)
	}
	fmt.Printf(" Saving %d nodes..\n", len(nds))
	f, err := os.Create(fp)
	if err != nil {
		return err
	}
	enc := gob.NewEncoder(f)
	err = enc.Encode(nds)
	if err != nil {
		return err
	}
	f.Close()

	gd.SaveAs(mmio.RemoveExtension(fp) + ".gdef")
	return nil
}

func LoadNetworkGob(fp string) ([]*Vert, error) {
	var nds []*Vert
	f, err := os.Open(fp)
	if err != nil {
		return nil, err
	}
	enc := gob.NewDecoder(f)
	err = enc.Decode(&nds)
	if err != nil {
		return nil, err
	}
	f.Close()
	return nds, nil
}
