package ptrack

import (
	"fmt"
	"log"
	"strings"

	"github.com/maseology/goHGS/structure"
)

func DirPrfx(fp string) string {
	return fp[:strings.Index(fp, "o.")]
}

// currently reads only 1 state.
func ReadHGS(q_pmFP string) Domain {
	prfx := q_pmFP
	if strings.Contains(q_pmFP, "o.") {
		prfx = q_pmFP[:strings.Index(q_pmFP, "o.")]
	}

	// get prisms
	h := structure.Read(prfx)
	fmt.Printf(" Model structure read: %d nodes, %d elements, %d layers\n", h.Nn, h.Ne, h.Nly)
	pset := func() PrismSet {
		prsms := make(map[int]*Prism, h.Ne)
		// p1---p2   y       0---nc
		//  | c |    |       |       clockwise, left-top-right-bottom
		// p0---p3   0---x   nr
		for eid, nds := range h.Exr {
			nh := len(nds) / 2
			z := make([]complex128, nh)
			top, bot := 0., 0.
			for i := 0; i < nh; i++ {
				xyz := h.Nxyz[nds[i]]
				z[nh-i-1] = complex(xyz[0], xyz[1]) // reverse order
				top += xyz[2]                       // averaging top elevation
			}
			for i := nh; i < 2*nh; i++ {
				xyz := h.Nxyz[nds[i]]
				bot += xyz[2] // averaging bottom elevation
			}
			prsms[eid] = &Prism{
				Z:   z,
				Top: top / float64(nh),
				Bot: bot / float64(nh),
				Por: defaultPorosity,
				Bn:  top / float64(nh),
				Tn:  0.,
			}
			prsms[eid].computeArea()
		}

		return PrismSet{
			P:    prsms,
			Conn: h.BuildTopology(),
		}
	}()

	// get flux
	t, v := h.ReadElementalVectors(q_pmFP)
	pflx := make(map[int]*PrismFlux)
	if len(v) == h.Ne {
		for i, vv := range v {
			pflx[i] = &PrismFlux{q: []float64{float64(vv[0]), float64(vv[1]), float64(vv[2])}} // centroid
			// fmt.Println(i, vv)
		}
	} else {
		log.Fatalln("ReadHGS todo")
	}

	// get saturation
	_, s := h.ReadNodelScalars(strings.Replace(q_pmFP, "q_pm", "head_pm", -1))
	if len(s) == h.Nn {
		for eid, nids := range h.Exr {
			ss := 0.
			for _, nid := range nids {
				ss += s[nid]
			}
			ss /= float64(len(nids))
			if ss < pset.P[eid].Top {
				if ss < pset.P[eid].Bot {
					pset.P[eid].Bn = pset.P[eid].Bot // dry cell
				} else {
					pset.P[eid].Bn = ss
				}
			}
		}
	} else {
		log.Fatalln("ReadHGS todo")
	}

	fmt.Printf(" results collected at time %f: %d vectors and %d scalars collected\n", t, len(v), len(s))

	var d Domain
	d.New(pset, pflx)
	return d
}
