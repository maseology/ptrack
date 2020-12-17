package ptrack

const rmax = 1.5 // threshold of Contains()

// VelocityFielder interface for flux field scheme
type VelocityFielder interface {
	PointVelocity(p *Particle, q *Prism, v float64) (float64, float64, float64)
	Contains(p *Particle) (float64, bool)
}
