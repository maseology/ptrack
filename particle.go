package ptrack

import (
	"encoding/json"
	"fmt"
	"log"
	"math"
	"os"

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
	mmio.DeleteFile("pl.geojson")
	fc := geojson.NewFeatureCollection()

	// as polylines
	for _, pln := range apl {
		pl := make([][]float64, len(pln))
		for j, p := range pln {
			pl[j] = []float64{p.X, p.Y, p.Z}
		}
		f := geojson.NewLineStringFeature(pl)
		f.SetProperty("pid", pln[0].I)
		fc.AddFeature(f)
	}

	// // as points
	// for _, pln := range apl {
	// 	for j, p := range pln {
	// 		f := geojson.NewPointFeature([]float64{p.X, p.Y, p.Z})
	// 		f.SetProperty("pid", p.I)
	// 		f.SetProperty("time", p.T/86400/365.24)
	// 		f.SetProperty("cid", p.C)
	// 		f.SetProperty("eid", p.C%epl)
	// 		f.SetProperty("vid", j)
	// 		fc.AddFeature(f)
	// 	}
	// }

	rawJSON, err := fc.MarshalJSON()
	if err != nil {
		log.Fatalf("MarshalJSON error: %v", err)
	}
	mmio.WriteString(fp, string(rawJSON)+"\n")
}

type pjson struct {
	X float64 `json:"x"`
	Y float64 `json:"y"`
	Z float64 `json:"z"`
	T float64 `json:"time"`
	K float64 `json:"k"`
	I int     `json:"particleid"`
}

// SaveJson saves a format that can be integrated with flopy's cross section plotter:
// xsect = flopy.plot.PlotCrossSection(model=gwf, line={"line": [(0, 5031.25), (10000, 5031.25)]})
// ptlns = pd.read_json('filename.json')
// ptlns = ptlns.sort_values(["particleid", "time"])
// xsect.plot_pathline(ptlns, ax=ax, colors='blue', linewidths=1, facecolors='none')
func SaveJson(fp string, apl [][]Particle) {
	var ptlns []pjson
	for j, pln := range apl {
		for _, p := range pln {
			pt := pjson{
				X: p.X,
				Y: p.Y,
				Z: p.Z,
				T: p.T,
				K: 0,
				I: j, //p.I,
			}
			ptlns = append(ptlns, pt)
		}
		// jsonData, err := json.MarshalIndent(ptl, "", "  ")
		// if err != nil {
		// 	panic(err)
		// }
		// fmt.Println(string(jsonData))

		// f, err := os.Create(fp)
		// if err != nil {
		// 	panic(err)
		// }
		// defer f.Close()

		// encoder := json.NewEncoder(f)
		// encoder.SetIndent("", "")

		// err = encoder.Encode(ptl)
		// if err != nil {
		// 	panic(err)
		// }
	}
	f, err := os.Create(fp)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	encoder := json.NewEncoder(f)
	encoder.SetIndent("", "")

	err = encoder.Encode(ptlns)
	if err != nil {
		panic(err)
	}
}
