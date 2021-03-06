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

// CentroidXY returns the coordinates of the prism centroid
func (q *Prism) CentroidXY() (x, y float64) {
	sc, c := 0.+0.i, 0
	for _, v := range q.Z {
		sc += v
		c++
	}
	ccxy := sc / complex(float64(c), 0.)
	return real(ccxy), imag(ccxy)
}

// // CentroidXY returns the complex-coordinates (x,y)=(real,imag) of the prism centroid
// func (q *Prism) CentroidXY() complex128 {
// 	sc, c := 0.+0.i, 0
// 	for _, v := range q.Z {
// 		sc += v
// 		c++
// 	}
// 	return sc / complex(float64(c), 0.)
// }

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
