package ptrack

import (
	"fmt"
	"math"
	"math/cmplx"
)

var mingtzero = 1.e-10

// Particle struct
type Particle struct {
	I          int
	X, Y, Z, T float64
}

// State returns the particles current state in CSV format
func (p *Particle) State() []float64 {
	return []float64{p.X, p.Y, p.Z, p.T}
}

// PrintState returns the particles current state in CSV format
func (p *Particle) PrintState() string {
	return fmt.Sprintf("%d,%v,%v,%v,%v", p.I, p.X, p.Y, p.Z, p.T)
}

// Dist returns the Euclidian distance between to points
func (p *Particle) Dist(p1 *Particle) float64 {
	return math.Sqrt(math.Pow(p.X-p1.X, 2.) + math.Pow(p.Y-p1.Y, 2.) + math.Pow(p.Z-p1.Z, 2.))
}

// Clone creates an exact copy of Particle
func (p *Particle) Clone() Particle {
	return Particle{I: p.I, X: p.X, Y: p.Y, Z: p.Z, T: p.T}
}

// ParticleTracker interface to particle tracking methods
type ParticleTracker interface {
	Track(p *Particle, q *Prism, vf VelocityFielder)
	// TrackToExit(p *Particle, q *Prism, vf VelocityFielder) []Particle
}

// VelocityFielder interface for flux field scheme
type VelocityFielder interface {
	PointVelocity(p *Particle, q *Prism, v float64) (float64, float64, float64)
}

// TrackToExit tracks particle to the exit point
func TrackToExit(p *Particle, q *Prism, w VelocityFielder, pt ParticleTracker, zw complex128) (int, [][]float64) {
	output, wx, wy, frst, ecode := make([][]float64, 1), real(zw), imag(zw), true, 0
	output[0] = p.State()
	for {
		pt.Track(p, q, w)
		if frst {
			frst = false
			if !q.ContainsVert(p) {
				// fmt.Println("early vertical exit")
				break
			}
			output = append(output, p.State())
		} else if q.Contains(p) {
			output = append(output, p.State())
		} else {
			// fmt.Println("particle exited cell")
			break
		}
		if !cmplx.IsNaN(zw) {
			if wDist(wx, wy, p.X, p.Y) < 0.1 {
				// fmt.Println("particle exited at well")
				ecode = -1
				break
			}
		}
	}
	if ecode == 0 {
		ecode = q.Intersection(p, output[len(output)-1])
	}
	output = append(output, p.State())

	return ecode, output
	// exit code (ecode) - nf: number of lateral faces:
	//  0..nf-1 - exited laterally, ecode = lateral face id
	//  nf      - exited bottom face
	//  nf+1    - exited top face
	//  <0      - exited from well (-well id)
	//  -9999   - error
}

func wDist(zwx, zwy, px, py float64) float64 {
	return math.Sqrt(math.Pow(zwx-px, 2.) + math.Pow(zwy-py, 2.))
}
