package ptrack

// VelocityFielder interface for flux field scheme
type VelocityFielder interface {
	PointVelocity(p *Particle, q *Prism, v float64) (float64, float64, float64)
	Local(p *Particle) (float64, bool)
	ReverseVectorField()
}
