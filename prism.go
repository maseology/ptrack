package ptrack

import (
	"math"

	"github.com/maseology/mmaths"
)

// const prsmTol1 = .01 // distance from cell top or bottom to remain "contained"

// PrismSet struct represents the model domain
type PrismSet struct {
	P    map[int]*Prism
	Conn map[int][]int
}

// Prism struct represents a singular model prism
type Prism struct {
	Z                        []complex128
	Top, Bot, A, Bn, Por, Tn float64
}

// New prism constructor
func (q *Prism) New(z []complex128, top, bot, bn, tn, porosity float64) {
	q.Z = z // complex coordinates
	q.Top = top
	q.Bot = bot
	q.Por = porosity
	q.Bn = bn // saturated thickness at time step tn
	q.Tn = tn // initial time step (both bn and tn will adjust in transient cases)
	q.computeArea()
}

func (q *Prism) computeArea() {
	q.A = 0.
	nfaces := len(q.Z)
	for j := range q.Z {
		jj := (j + 1) % nfaces
		q.A += real(q.Z[j])*imag(q.Z[jj]) - real(q.Z[jj])*imag(q.Z[j])
	}
	q.A /= -2. // negative used here because vertices are entered in clockwise order
	if q.A <= 0. {
		panic("vertex error, may be given in counter-clockwise order")
	}
}

// CentroidXY returns the complex-coordinates (x,y)=(real,imag) of the prism centroid
func (q *Prism) CentroidXY() complex128 {
	sc, c := 0.+0.i, 0
	for _, v := range q.Z {
		sc += v
		c++
	}
	return sc / complex(float64(c), 0.)
}

// getExtentsXY returns the XY-extents of the prism
func (q *Prism) getExtentsXY() (yn, yx, xn, xx float64) {
	yn, yx, xn, xx = math.MaxFloat64, -math.MaxFloat64, math.MaxFloat64, -math.MaxFloat64
	for _, v := range q.Z {
		yn = math.Min(yn, imag(v))
		yx = math.Max(yx, imag(v))
		xn = math.Min(xn, real(v))
		xx = math.Max(xx, real(v))
	}
	return
}

const tol = 1e-10

// ContainsXY returns true if the given (x,y) coordinates are contained by the prism planform bounds
func (q *Prism) ContainsXY(x, y float64) bool {
	return mmaths.PnPolyC(q.Z, complex(x, y), tol)
}

// Contains returns true if the given particle is contained by the prism bounds
func (q *Prism) Contains(p *Particle) bool {
	if !mmaths.PnPolyC(q.Z, complex(p.X, p.Y), tol) {
		return false
	}
	return p.Z <= q.Top && p.Z >= q.Bot
}

// // ContainsVert returns true if the given particle is contained by the prism's vertical bounds
// func (q *Prism) ContainsVert(p *Particle) bool {
// 	return p.Z <= q.Top && p.Z >= q.Bot
// }

// // Intersection determines the particles exit point relative to its last position
// func (q *Prism) Intersection(p *Particle, lastpos []float64) int {
// 	ecode, nfaces, f := 0, len(q.Z), -1.
// 	x1, y1, x2, y2 := lastpos[0], lastpos[1], p.X, p.Y
// 	for j := range q.Z {
// 		jj := (j + 1) % nfaces
// 		x3, y3, x4, y4 := real(q.Z[j]), imag(q.Z[j]), real(q.Z[jj]), imag(q.Z[jj])
// 		f = intersection2D(x1, y1, x2, y2, x3, y3, x4, y4)
// 		if !math.IsNaN(f) {
// 			ecode = j
// 			break
// 		}
// 	}
// 	if p.Z <= q.Top && p.Z >= q.Bot { // still contained vertically
// 		if math.IsNaN(f) {
// 			if math.Abs(p.Z-q.Top) < prsmTol1 {
// 				return nfaces + 1
// 			} else if math.Abs(p.Z-q.Bot) < prsmTol1 {
// 				return nfaces
// 			} else {
// 				fmt.Println(*p)
// 				return -9999
// 				// panic("Prism.Intersection error: particle has not appeared to exit prism.")
// 			}
// 		}
// 	} else { // possibly exited from a vertical face
// 		var f2 float64
// 		ec2 := nfaces
// 		if p.Z > q.Top {
// 			f2 = (q.Top - lastpos[2]) / (p.Z - lastpos[2])
// 			ec2++
// 		} else {
// 			f2 = (q.Bot - lastpos[2]) / (p.Z - lastpos[2])
// 		}
// 		if f2 < 0. || f2 > 1. {
// 			if math.Abs(p.Z-q.Top) < prsmTol1 {
// 				ec2++
// 			} else if math.Abs(p.Z-q.Bot) < prsmTol1 {
// 				// do nothing, exiting bottom
// 			} else {
// 				fmt.Println(p.Z, q.Top, q.Bot, lastpos[2], f)
// 				return -9999
// 				// panic("Prism.Intersection error 001")
// 			}
// 		}
// 		if math.IsNaN(f) || f2 < f { // particle exited top or bottom face
// 			f = f2
// 			ecode = ec2
// 		}
// 	}
// 	// f *= 1.0001 // slightly over project particle to ensure it's found in the following cell
// 	// p.X = lastpos[0] + f*(p.X-lastpos[0])
// 	// p.Y = lastpos[1] + f*(p.Y-lastpos[1])
// 	// p.Z = lastpos[2] + f*(p.Z-lastpos[2])
// 	// p.T = lastpos[3] + f*(p.T-lastpos[3])

// 	return ecode
// 	// exit code (ecode) - nf: number of lateral faces:
// 	//  0..nf-1 - exited laterally, ecode = lateral face id
// 	//  nf      - exited bottom face
// 	//  nf+1    - exited top face
// 	//  <0      - exited from well (-well id)
// 	//  -9999   - error
// }

// func intersection2D(x1, y1, x2, y2, x3, y3, x4, y4 float64) float64 {

// 	sx1 := x2 - x1
// 	sy1 := y2 - y1
// 	sx2 := x4 - x3
// 	sy2 := y4 - y3

// 	s := (-sy1*(x1-x3) + sx1*(y1-y3)) / (-sx2*sy1 + sx1*sy2)
// 	t := (sx2*(y1-y3) - sy2*(x1-x3)) / (-sx2*sy1 + sx1*sy2)

// 	if s >= 0 && s <= 1 && t >= 0 && t <= 1 {
// 		// Collision detected
// 		// return []float64{x1 + (t * sx1), y1 + (t * sy1)}
// 		return t
// 	}

// 	// // taken from LineSegment.Intersection2D
// 	// // first degree Bézier parameter
// 	// d := ((x1-x2)*(y3-y4) - (y1-y2)*(x3-x4))
// 	// t := ((x1-x3)*(y3-y4) - (y1-y3)*(x3-x4))
// 	// u := ((x1-x2)*(y1-y3) - (y1-y2)*(x1-x3))
// 	// t /= d
// 	// u /= -d
// 	// if t >= 0. && t <= 1. && u >= 0. && u <= 1. {
// 	// 	return t
// 	// }
// 	return math.NaN() // No collision/line segments do not intersect
// }
