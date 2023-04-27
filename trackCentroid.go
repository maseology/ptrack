package ptrack

import (
	"math"
	"sort"
)

// Track a collection of particles through the centroid of at least 1 model cell
func (d *Domain) TrackCentroidalParticles(excl map[int]bool) ([][]Particle, int, []int) {

	chknan := func(a []Particle) ([]Particle, bool) {
		rm, fxd := []int{}, false
		for i, aa := range a {
			// for _, aaa := range aa {
			// 	if math.IsNaN(aaa) {
			// 		fxd = true
			// 		rm = append(rm, i)
			// 		break
			// 	}
			// }
			if math.IsNaN(aa.X) || math.IsNaN(aa.Y) || math.IsNaN(aa.Z) || math.IsNaN(aa.T) {
				fxd = true
				rm = append(rm, i)
				break
			}
		}
		if fxd {
			sort.Ints(rm)
			for i, j := 0, len(rm)-1; i < j; i, j = i+1, j-1 {
				rm[i], rm[j] = rm[j], rm[i]
			}
			for _, i := range rm {
				if i == len(a)-1 {
					a = a[:i]
				} else {
					a = append(a[:i], a[i+1:]...)
				}
			}
		}
		return a, fxd
	}

	np := d.Nprism() - len(excl)
	o, pxr, c, k := make([][]Particle, np), make([]int, np), 0, 0
	for pid, p := range d.prsms {
		if excl[pid] {
			continue
		}
		a := d.trackParticle(p.CentroidParticle(pid), pid)
		if x, ok := chknan(a); ok {
			a = x
		}
		o[k] = a
		pxr[k] = pid
		c += len(o[k])
		k++
		// if k == 1000 {
		// 	break
		// }
	}

	// reverse tracks
	println(" reversing flux filed..")
	d.ReverseVectorField()
	for k, pid := range pxr {
		ar := d.trackParticle(d.prsms[pid].CentroidParticle(pid), pid)
		if x, ok := chknan(ar); ok {
			ar = x
		}
		c += len(ar)
		for i, j := 0, len(ar)-1; i < j; i, j = i+1, j-1 {
			ar[i], ar[j] = ar[j], ar[i] // reverse array
		}
		// for i := 0; i < len(o); i++ {
		// 	for j := 0; j < len(o[i]); j++ {
		// 		o[i][j][3] = -o[i][j][3] // reverse tracking time
		// 	}
		// }
		o[k] = append(ar[:len(ar)-1], o[k]...)
		// if k == 1000 {
		// 	break
		// }
	}

	// func() {
	// 	f, err := os.Create("o.gob")
	// 	if err != nil {
	// 		panic(err)
	// 	}
	// 	enc := gob.NewEncoder(f)
	// 	err = enc.Encode(o)
	// 	if err != nil {
	// 		panic(err)
	// 	}
	// 	f.Close()
	// }()
	// func() {
	// 	f, err := os.Create("pxr.gob")
	// 	if err != nil {
	// 		panic(err)
	// 	}
	// 	enc := gob.NewEncoder(f)
	// 	err = enc.Encode(pxr)
	// 	if err != nil {
	// 		panic(err)
	// 	}
	// 	f.Close()
	// }()
	// func() {
	// 	f, err := os.Create("c.gob")
	// 	if err != nil {
	// 		panic(err)
	// 	}
	// 	enc := gob.NewEncoder(f)
	// 	err = enc.Encode(c)
	// 	if err != nil {
	// 		panic(err)
	// 	}
	// 	f.Close()
	// }()

	return o, c, pxr
}
