package ptrack

// Track a collection of particles through the centroid of at least 1 model cell
func (d *Domain) TrackCentroidalParticles() ([][][]float64, int, []int) {
	o, pxr, c, k := make([][][]float64, d.Nprism()), make([]int, d.Nprism()), 0, 0

	for pid, p := range d.prsms {
		o[k] = d.trackParticle(p.CentroidParticle(pid), pid)
		pxr[k] = pid
		c += len(o[k])
		k++
	}

	// pid := 447380 //1506355
	// o[k] = d.trackParticle(d.prsms[pid].CentroidParticle(pid), pid)
	// c += len(o[k])

	return o, c, pxr
}
