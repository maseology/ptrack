package ptrack

import (
	"fmt"
	"math"
	"math/cmplx"
)

// ParticleTracker interface to particle tracking methods
type ParticleTracker interface {
	track(done <-chan interface{}, p *Particle, q *Prism, vf VelocityFielder) <-chan []float64
}

// TrackToExit tracks particle to the exit point
func trackToExit(p *Particle, q *Prism, w VelocityFielder, pt ParticleTracker, zwell complex128) [][]float64 {
	pcoll, iswell := [][]float64{p.State()}, !cmplx.IsNaN(zwell)
	wx, wy := real(zwell), imag(zwell)
	done := make(chan interface{})
	switch pt.(type) {
	case *PollockMethod:
		pt = w.(*PollockMethod) // type assertion
		if iswell {
			fmt.Println("caution particle has been place in a cell with well, may not be compatible with Pollock Method")
		}
	}
	for pstate := range pt.track(done, p, q, w) {
		pcoll = append(pcoll, pstate)
		if !q.Contains(p) {
			break
		}
		if iswell {
			if wDist(wx, wy, pstate[0], pstate[1]) < wellTol {
				// fmt.Println("particle exited at well")
				break
			}
		}
	}
	close(done)
	return pcoll
}

func wDist(zwx, zwy, px, py float64) float64 {
	return math.Sqrt(math.Pow(zwx-px, 2.) + math.Pow(zwy-py, 2.))
}
