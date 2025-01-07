package ptrack

// Track a collection of particles through the centroid of at least 1 model cell
func (d *Domain) TrackCentroidalParticlesReverse(prnt bool) ([][]Particle, int, []int) {
	o, pxr, c, k := make([][]Particle, d.Nprism()), make([]int, d.Nprism()), 0, 0
	d.ReverseVectorField()
	for pid, p := range d.prsms {
		o[k] = d.trackParticle(p.CentroidParticle(pid), pid, prnt)
		pxr[k] = pid
		c += len(o[k])
		k++
	}

	// // reverse tracks
	// d.ReverseVectorField()
	// for k, pid := range pxr {
	// 	ar := d.trackParticle(d.prsms[pid].CentroidParticle(pid), pid)
	// 	c += len(ar)
	// 	for i, j := 0, len(ar)-1; i < j; i, j = i+1, j-1 {
	// 		ar[i], ar[j] = ar[j], ar[i] // reverse array
	// 	}
	// 	o[k] = append(ar, o[k]...)
	// }

	return o, c, pxr
}
