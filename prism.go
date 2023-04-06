package ptrack

import (
	"fmt"
	"math"

	"github.com/maseology/mmaths"
)

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

// Function IsClockwise() As Boolean
// Dim db1 As Double = 0.0
// For i = 0 To _v.Count - 2
// 	db1 += (_v(i + 1).X - _v(i).X) * (_v(i + 1).Y + _v(i).Y)
// Next
// db1 += (_v(0).X - _v.Last.X) * (_v(0).Y + _v.Last.Y)
// Return db1 >= 0.0
// End Function

func (q *Prism) computeArea() {
	nfaces := len(q.Z)
	q.A = 0.
	xo, yo := real(q.Z[0]), imag(q.Z[0])
	for j, z := range q.Z {
		jj := (j + 1) % nfaces
		q.A += (real(z)-xo)*(imag(q.Z[jj])-yo) - (real(q.Z[jj])-xo)*(imag(z)-yo)
	}
	q.A /= -2. // negative used here because vertices are entered in clockwise order
	if q.A <= 0. {
		for i, z := range q.Z {
			fmt.Printf("%d,%f,%f\n", i+1, real(z), imag(z))
		}
		panic("vertex error, may be given in counter-clockwise order")
	}
}

func (q *Prism) Saturation() float64 {
	return (q.Bn - q.Bot) / (q.Top - q.Bot)
}

// Centroid returns the coordinates of the prism centroid
func (q *Prism) Centroid() (x, y, z float64) {
	x, y = q.CentroidXY()
	z = (q.Top-q.Bot)/2. + q.Bot
	return
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

// CentroidParticle returns the coordinates of the prism centroid
func (q *Prism) CentroidParticle(i int) *Particle {
	x, y, z := q.Centroid()
	return &Particle{I: i, C: i, X: x, Y: y, Z: z}
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

// ContainsXY returns true if the given (x,y) coordinates are contained by the prism planform bounds
func (q *Prism) ContainsXY(x, y float64) bool {
	return mmaths.PnPolyC(q.Z, complex(x, y), tol)
}

// ContainsXYZ returns true if the given (x,y) coordinates are contained by the prism planform bounds
func (q *Prism) ContainsXYZ(x, y, z float64) bool {
	if !mmaths.PnPolyC(q.Z, complex(x, y), tol) {
		return false
	}
	return z <= q.Top && z >= q.Bot
}

// Contains returns true if the given particle is contained by the prism bounds
func (q *Prism) Contains(p *Particle) bool {
	if !mmaths.PnPolyC(q.Z, complex(p.X, p.Y), tol) {
		return false
	}
	return p.Z <= q.Top && p.Z >= q.Bot
}
