package ptrack

// ParticleTracker interface to particle tracking methods
type ParticleTracker interface {
	track(done <-chan interface{}, p *Particle, q *Prism, vf VelocityFielder) <-chan Particle
}

// trackToPrismExit tracks particle to the exit point of a prism
func trackToPrismExit(p *Particle, q *Prism, w VelocityFielder, pt ParticleTracker) []Particle {
	//panic("to check")
	// pcoll := [][]float64{p.State()}
	pcoll := []Particle{*p}
	done := make(chan interface{})
	for pstate := range pt.track(done, p, q, w) {
		pcoll = append(pcoll, pstate)
		if !q.Contains(p) {
			break // left prism
		}
	}
	close(done)
	return pcoll
}
