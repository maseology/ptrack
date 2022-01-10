package ptrack

// ParticleTracker interface to particle tracking methods
type ParticleTracker interface {
	track(done <-chan interface{}, p *Particle, q *Prism, vf VelocityFielder) <-chan []float64
	// TracktoExit(p *Particle, q *Prism, vf VelocityFielder) [][]float64
}

// TrackToExit tracks particle to the exit point
func trackToExit(p *Particle, q *Prism, w VelocityFielder, pt ParticleTracker) [][]float64 {
	// return pt.testTracktoExit(p, q, w) //  for testing (not concurrent)
	pcoll := [][]float64{p.State()}
	done := make(chan interface{})
	for pstate := range pt.track(done, p, q, w) {
		pcoll = append(pcoll, pstate)
		if !q.Contains(p) {
			break
		}
	}
	close(done)
	return pcoll
}
