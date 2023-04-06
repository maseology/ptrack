package ptrack

import (
	"math"
	"math/cmplx"

	"github.com/maseology/mmaths"
)

// VectorMethSoln uses a uniform vector within a prism
type VectorMethSoln struct {
	zc                  complex128
	r, vx, vy, vzt, vzb float64
}

func (vm *VectorMethSoln) New(zc complex128, q []float64, por, r float64) { // left-up-right-down-bottom-top
	vm.zc = zc
	vm.r = r
	switch len(q) {
	case 3:
		vm.vx = q[0] / por
		vm.vy = q[1] / por
		vm.vzb = q[2] / por
		vm.vzt = q[2] / por
	case 6:
		vm.vx = (q[0] - q[2]) / 2. / por
		vm.vy = (q[3] - q[1]) / 2. / por
		vm.vzb = q[4] / por
		vm.vzt = -q[5] / por
	default:
		panic("New VectorMethSoln error")
	}
}

// PointVelocity returns the velocity vector for a given (x,y,z) coordinate
func (vm *VectorMethSoln) PointVelocity(d1 *Particle, d2 *Prism, d3 float64) (float64, float64, float64) {
	x := (d1.Z - d2.Bot) / (d2.Top - d2.Bot)
	return vm.vx, vm.vy, x*vm.vzt + (x-1.)*vm.vzb
}

// Local returns whether the point is solvable within the solution space
func (vm *VectorMethSoln) Local(p *Particle) (float64, bool) {
	azl := cmplx.Abs(complex(p.X, p.Y)-vm.zc) / vm.r // relative coordinate
	return azl, azl <= 1.
}

func (vm *VectorMethSoln) track(done <-chan interface{}, p *Particle, q *Prism, vf VelocityFielder) <-chan Particle {
	chout := make(chan Particle)
	go func() {
		defer close(chout)
		for {
			select {
			case <-done:
				return
			default:

				// get relative angles [-pi,pi]
				cp := complex(p.X, p.Y)

				vx, vy, vz := vf.PointVelocity(p, q, 0.)
				cv := complex(vx, vy)
				av := cmplx.Abs(cv)
				// vtheta := cmplx.Phase(cv)
				cv *= complex(3*vm.r/av, 0.)
				cv += cp

				// check lateral exit
				xi, yi := func() (float64, float64) {
					for i := range q.Z {
						ii := (i + 1) % len(q.Z)
						p0 := mmaths.Point{X: real(q.Z[ii]), Y: imag(q.Z[ii])}
						p1 := mmaths.Point{X: real(q.Z[i]), Y: imag(q.Z[i])}
						ls := mmaths.LineSegment{P0: p0, P1: p1}
						pA := mmaths.Point{X: p.X, Y: p.Y}
						pB := mmaths.Point{X: real(cv), Y: imag(cv)}
						pr := mmaths.LineSegment{P0: pA, P1: pB}
						xy, f := ls.Intersection2D(&pr)
						if math.IsNaN(f) {
							continue
						}
						return xy.X, xy.Y
					}
					return math.NaN(), math.NaN()
				}()
				d := math.Sqrt(math.Pow(xi-p.X, 2) + math.Pow(yi-p.Y, 2))
				te := d / av * 1.00001 // adding a little momentum to nudge the particle past the boundary
				xe := p.X + vx*te
				ye := p.Y + vy*te
				ze := p.Z + vz*te

				// check vertical exit
				if ze > q.Top {
					tez := (q.Top - p.Z) / vz
					// if tez < te {
					xe = tez/te*(xe-p.X) + p.X
					ye = tez/te*(ye-p.Y) + p.Y
					ze = q.Top
					te = tez
					// }
				} else if ze < q.Bot {
					tez := (q.Bot - p.Z) / vz
					// if tez < te {
					xe = tez/te*(xe-p.X) + p.X
					ye = tez/te*(ye-p.Y) + p.Y
					ze = q.Bot
					te = tez
					// }
				}

				// exit point
				p.X = xe
				p.Y = ye
				p.Z = ze
				p.T += te
				chout <- *p //.State()
				return
			}
		}
	}()
	return chout
}

func (vm *VectorMethSoln) ReverseVectorField() {
	vm.vx *= -1.
	vm.vy *= -1.
	vm.vzb *= -1.
	vm.vzt *= -1.
}
