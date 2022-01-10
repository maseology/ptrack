package ptrack

import (
	"fmt"
	"math"
)

// Particles is a collection of Particle
type Particles []Particle

// Particle struct
type Particle struct {
	I          int
	X, Y, Z, T float64
}

// State returns the particles current state in CSV format
func (p *Particle) State() []float64 {
	return []float64{p.X, p.Y, p.Z, p.T}
}

// ResetState returns the particles current state in CSV format
func (p *Particle) ResetState(s []float64) {
	p.X = s[0]
	p.Y = s[1]
	p.Z = s[2]
	p.T = s[3]
}

// PrintState returns the particles current state in CSV format
func (p *Particle) PrintState() string {
	return fmt.Sprintf("%d,%v,%v,%v,%v", p.I, p.X, p.Y, p.Z, p.T)
}

// Dist returns the Euclidian distance between to points
func (p *Particle) Dist(p1 *Particle) float64 {
	return math.Sqrt(math.Pow(p.X-p1.X, 2.) + math.Pow(p.Y-p1.Y, 2.) + math.Pow(p.Z-p1.Z, 2.))
}
