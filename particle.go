package ptrack

import (
	"fmt"
	"log"
	"math"

	"github.com/maseology/mmio"
	geojson "github.com/paulmach/go.geojson"
)

// Particles is a collection of Particle
type Particles []Particle

// Particle struct
type Particle struct {
	I, C       int
	X, Y, Z, T float64
}

// func (p *Particle) State() []float64 {
// 	return []float64{p.X, p.Y, p.Z, p.T}
// }

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

func SaveGeojson(fp string, apl [][]Particle, epl int) {
	fc := geojson.NewFeatureCollection()

	// // as polylines
	// for i, pln := range apl {
	// 	f := geojson.NewLineStringFeature(pln)
	// 	f.SetProperty("pid", pxr[i])
	// 	fc.AddFeature(f)
	// }

	// as points
	for _, pln := range apl {
		for j, p := range pln {
			f := geojson.NewPointFeature([]float64{p.X, p.Y, p.Z})
			f.SetProperty("pid", p.I)
			f.SetProperty("time", p.T/86400/365.24)
			f.SetProperty("cid", p.C)
			f.SetProperty("eid", p.C%epl)
			f.SetProperty("vid", j)
			fc.AddFeature(f)
		}
	}

	rawJSON, err := fc.MarshalJSON()
	if err != nil {
		log.Fatalf("MarshalJSON error: %v", err)
	}
	mmio.WriteString(fp, string(rawJSON)+"\n")
}
