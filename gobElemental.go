package ptrack

import (
	"encoding/gob"
	"fmt"
	"os"
)

func (d *Domain) ExportMeshNetworkGob(fp string, pl [][]Particle, xr []int) error {

	// build topology
	var nds []Vert
	nc, k := len(d.Prisms())/d.Nly, 0
	for i, pln := range pl {
		pnds := make([]Vert, len(pln))
		// pidlast := -1
		for j, v := range pln {
			// pid, cid, ly := func() (int, int, int) {
			// 	eid := mh.PointToElementID(v.X, v.Y, searchradius) // SLOW
			// 	ly := 0
			// 	for {
			// 		if p, ok := prsms[eid+ly*nc]; ok {
			// 			if v.Z <= p.Top && v.Z >= p.Bot {
			// 				return eid + ly*nc, eid, ly + 1
			// 			}
			// 		} else {
			// 			break
			// 		}
			// 		ly++
			// 	}
			// 	return -1, eid, ly
			// }()
			// pnds[j] = Vert{
			// 	X:      v.X,
			// 	Y:      v.Y,
			// 	Z:      v.Z,
			// 	T:      v.T,
			// 	VertID: k,
			// 	PathID: xr[i],
			// 	PrsmID: pid, // the prism/cell it's in
			// 	CellID: cid, // the 2d grid cell/element it's in
			// 	Layer:  ly,
			// 	Order:  j,
			// 	USid:   -1,
			// 	DSid:   -1,
			// }

			pnds[j] = Vert{
				X:      v.X,
				Y:      v.Y,
				Z:      v.Z,
				T:      v.T,
				VertID: k,
				PathID: xr[i],
				PrsmID: v.C,      // the prism/cell it's in
				CellID: v.C % nc, // the 2d grid cell/element it's in
				Layer:  v.C/nc + 1,
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

	fmt.Printf(" Saving %s nodes..\n", big(len(nds)))
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

	// mh.SaveAs(strings.TrimSuffix(fp, filepath.Ext(fp)) + ".ah2") // RemoveExtension
	return nil
}
