package ptrack

import (
	"math"
	"math/cmplx"
)

// ParticleTracker interface to particle tracking methods
type ParticleTracker interface {
	Track(p *Particle, q *Prism, vf VelocityFielder)
	// TrackToExit(p *Particle, q *Prism, vf VelocityFielder) []Particle
}

// TrackToExit tracks particle to the exit point
func TrackToExit(p *Particle, q *Prism, w VelocityFielder, pt ParticleTracker, zwell complex128, wellID int) (int, [][]float64) {
	pcoll := [][]float64{p.State()}
	wx, wy, exitcode := real(zwell), imag(zwell), 0
	for {
		pt.Track(p, q, w)
		pcoll = append(pcoll, p.State())
		if !q.Contains(p) {
			exitcode = q.Intersection(p, pcoll[len(pcoll)-2])
			break
		}
		if !cmplx.IsNaN(zwell) {
			if wDist(wx, wy, p.X, p.Y) < wellTol {
				// fmt.Println("particle exited at well")
				exitcode = -wellID // currently defaulted to centroid wells id'd by their cellID
				break
			}
		}
	}
	return exitcode, pcoll
	// exit code (exitcode) - nf: number of lateral faces:
	//  0..nf-1 - exited laterally, ecode = lateral face id
	//  nf      - exited bottom face
	//  nf+1    - exited top face
	//  <0      - exited from well (-well id)
	//  -9999   - error

	// output, wx, wy, frst, ecode := make([][]float64, 1), real(zwell), imag(zwell), true, 0
	// output[0] = p.State()
	// for {
	// 	pt.Track(p, q, w)
	// 	if frst {
	// 		frst = false
	// 		if !(p.Z <= q.Top && p.Z >= q.Bot) {
	// 			// fmt.Println("early vertical exit")
	// 			break
	// 		}
	// 		output = append(output, p.State())
	// 	} else if q.Contains(p) {
	// 		output = append(output, p.State())
	// 	} else {
	// 		// fmt.Println("particle exited cell")
	// 		output = append(output, p.State())
	// 		break
	// 	}
	// 	if !cmplx.IsNaN(zwell) {
	// 		if wDist(wx, wy, p.X, p.Y) < wellTol {
	// 			// fmt.Println("particle exited at well")
	// 			ecode = -wellID // currently defaulted to centroid wells id'd by their cellID
	// 			break
	// 		}
	// 	}
	// }
	// if ecode == 0 {
	// 	ecode = q.Intersection(p, output[len(output)-2])
	// }
	// // output = append(output, p.State())

	// return ecode, output
	// // exit code (ecode) - nf: number of lateral faces:
	// //  0..nf-1 - exited laterally, ecode = lateral face id
	// //  nf      - exited bottom face
	// //  nf+1    - exited top face
	// //  <0      - exited from well (-well id)
	// //  -9999   - error
}

func wDist(zwx, zwy, px, py float64) float64 {
	return math.Sqrt(math.Pow(zwx-px, 2.) + math.Pow(zwy-py, 2.))
}
