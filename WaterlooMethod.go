package ptrack

import (
	"fmt"
	"log"
	"math"
	"math/cmplx"

	"github.com/maseology/mmio"
	"gonum.org/v1/gonum/mat"
)

// WatMethSoln (The Waterloo Method Solution): is a struct that contains a grid cell's internal flow field
// Ramadhan, M., 2015, A Semi-Analytical Particle Tracking Algorithm for Arbitrary Unstructured Grids. Unpublished MASc. Thesis. University of Waterloo, Waterloo Ontario.
type WatMethSoln struct {
	aT, zwl        []complex128
	zc, r          complex128
	qv, qw, ql, qb float64
	m, n, nf       int // n control points; order of approximiation
}

// New WatMethSoln constructor
func (w *WatMethSoln) New(prismID int, p *Prism, Qj []float64, zw complex128, Qtop, Qbot, Qwell float64, m, n int, prnt bool) {
	w.m = m      // Total control points
	w.n = n      // Order of Approximation
	if m < 2*n { // constraint
		log.Fatal("control point specification error: m < 2*n")
	}
	w.ql = 0.
	w.nf = len(p.Z) // n faces
	w.zwl = make([]complex128, 1)
	wbal, Qvert := 0., Qtop-Qbot
	for j := 0; j < w.nf; j++ {
		wbal += Qj[j] // flux along face j; positive in, negative out (left-up-right-down)
		w.zc += p.Z[j]
	}
	w.ql = wbal // cumulative lateral inflows
	wbal += Qwell - Qvert
	if math.Abs(wbal)/w.ql > mingtzero { // step 1: check mass balance (eq 3.7)
		fmt.Printf("Qwell: %6.3f  Qbot: %6.3f  Qtop: %6.3f  Qlat %6.3f\n", Qwell, Qbot, Qtop, Qj)
		log.Fatalf("cell mass-balance error in prism %d: %v\n", prismID, wbal)
	}
	w.zc /= complex(float64(w.nf), 0.) // cell centroid
	w.qv = Qvert / p.Area              // Qvert = Qtop-Qbot, positive up, negative down
	w.qw = Qwell / 2. / math.Pi

	w.buildCoefTaylor(p.Z, Qj, prnt) // (p.Z, Qj, Qvert == 0. && Qwell == 0.)

	w.qb = Qbot / p.Area
	w.zwl[0] = (zw - w.zc) / w.r // local well coordinate
	if cmplx.IsNaN(zw) && Qwell != 0. {
		log.Panicln("need to specify well coordinate")
	}
}

func (w *WatMethSoln) buildCoefTaylor(zj []complex128, qj []float64, prnt bool) {
	// The Waterloo method
	// Muhammad Ramadhan, 2015. a Semi-Analytic Particle Tracking Algorithm for Arbitrary Unstructured Grids. M.A.Sc thesis. University of Waterloo.

	// step 2: compute cell face angles, lengths, (planform) area, perimeter, maximum radius
	alpha, lj, p, r := make([]float64, w.nf), make([]float64, w.nf), 0., 0.
	for j := 0; j < w.nf; j++ {
		jj := (j + 1) % w.nf
		d := zj[jj] - zj[j]
		lj[j] = cmplx.Abs(d) // length of each side of the cell [L]
		p += lj[j]           // total length = perimeter
		r = math.Max(r, cmplx.Abs(w.zc-zj[j]))
		alpha[j] = cmplx.Phase(complex(imag(d), -real(d))) // angle for each side of the polygon calculated in radian (table 3.6)
	}
	w.r = complex(r, 0.)

	// step 3: build control points
	sCtrl, zCtrl, ijx := make([]float64, w.m), make([]complex128, w.m), make([]int, w.m)
	pPrev, pNext, j, sc := 0., lj[0], 0, p/float64(w.m)
	for i := 0; i < w.m; i++ {
		sCtrl[i] = (float64(i) + 0.5) * sc // eq. 3.6a
		if sCtrl[i] > pNext {
			j++
			pPrev = pNext
			pNext += lj[j]
		}
		ijx[i] = j // control point to face cross-reference
		jnext := (j + 1) % w.nf
		zCtrlGlobal := complex((sCtrl[i]-pPrev)/lj[j], 0.)*(zj[jnext]-zj[j]) + zj[j] // eq. 3.6b
		zCtrl[i] = (zCtrlGlobal - w.zc) / w.r                                        // local control point (Z_ctrl) eq. 3.6c
	}
	if prnt {
		w.saveControlPoints(zCtrl)
	}

	// step 4: determine normalized cell flows: Qcell, qVert and qWell, and build Qtaylor
	qTaylorni := make([]float64, w.m)
	for i := 0; i < w.m; i++ {
		j := ijx[i]                                                        // side ith control point on the jth face
		ca := complex(math.Cos(alpha[j]), math.Sin(alpha[j]))              // complex angle
		qCellni := qj[j] / lj[j]                                           // eq. 3.8 normalized cell flows
		qVertni := real(complex(-w.qv*r*real(zCtrl[i]), 0.) * ca)          // eq. 3.11a
		qWellni := real(-ca * complex(w.qw/r, 0.) / (zCtrl[i] - w.zwl[0])) // eq. 3.11b
		qTaylorni[i] = qCellni - qVertni - qWellni                         // eq. 3.10
	}

	// step 5: build Phi_Taylor
	psiTaylor := make([]float64, w.m)
	psiTaylor[0] = qTaylorni[0] * sCtrl[0]
	for i := 1; i < w.m; i++ {
		psiTaylor[i] = 0.5*(qTaylorni[i]+qTaylorni[i-1])*(sCtrl[i]-sCtrl[i-1]) + psiTaylor[i-1] // eq.3.12b
	}

	// step 6: build matrices to solve [phiU]*[a]=[phyTaylor] (eq.3.15)
	psiTaylorU := make([]float64, 2*w.m*w.n)
	for i := 0; i < w.m; i++ {
		for j := 0; j < w.n; j++ {
			jf := complex(float64(j), 0.)
			c := cmplx.Pow(zCtrl[i], jf)
			psiTaylorU[i*2*w.n+j] = imag(c)
			psiTaylorU[i*2*w.n+j+w.n] = real(c)
		}
	}
	psiU := mat.NewDense(w.m, 2*w.n, psiTaylorU)
	psiT := mat.NewVecDense(w.m, psiTaylor)
	aTaylorV := svdSolve(psiU, psiT)

	aTaylor := make([]complex128, w.n)
	for j := 0; j < w.n; j++ {
		aTaylor[j] = complex(aTaylorV.At(j, 0), aTaylorV.At(j+w.n, 0))
	}
	w.aT = aTaylor

	if prnt {
		w.plotPerimeterFlux(zj, qj, lj, p, 500) //, w.m) //, 500) //
	}
}

// PointVelocity returns the velocity vector for a given (x,y,z) coordinate. ** Set dbdt = 0. for steady-state cases
func (w *WatMethSoln) PointVelocity(p *Particle, q *Prism, dbdt float64) (float64, float64, float64) {
	// q.Bn: saturated thickness at beginning of time step at time q.Tn; dbdt: rate of change in saturated thickness
	bl, bz := math.Min(q.Bn, q.Top-q.Bot), math.Min(q.Bn+(p.T-q.Tn)*dbdt, q.Top-q.Bot) // corrected saturated thickness for lateral and vertical computations
	zl := (complex(p.X, p.Y) - w.zc) / w.r                                             // complex local coordinate
	o := w.cmplxVelFlow(zl)
	if w.qv != 0. {
		o += w.cmplxVelVert(zl)
	}
	if w.qw != 0. && zl != 0 {
		o += w.cmplxVelWell(zl, 0)
	}
	// vz := (w.qb + (p.Z-q.Bot)*w.ql/q.Area/bl) / q.Por     // eq. 3.19 (steady-state case)
	vz := (w.qb + (p.Z-q.Bot)*w.ql/q.Area/bz) / q.Por // eq. 3.18 (transient case)
	vx := real(-o) / bl / q.Por
	vy := imag(o) / bl / q.Por
	// fmt.Println("  vel", vx, vy)
	return vx, vy, vz // eq. 3.17
}

// Local returns whether the point is solvable within the solution space
func (w *WatMethSoln) Local(p *Particle) (float64, bool) {
	zl := (complex(p.X, p.Y) - w.zc) / w.r // complex local coordinate
	azl := cmplx.Abs(zl)                   // relative coordinate
	return azl, azl <= 1.
}

func (w *WatMethSoln) ReverseVectorField() {
	println("TODO WatMethSoln.ReverseVectorField")
}

func (w *WatMethSoln) saveControlPoints(zCtrl []complex128) {
	txtw, _ := mmio.NewTXTwriter("surfer-controlpoints.bln")
	for i := 0; i < w.m; i++ {
		zCtrlGlobal := zCtrl[i]*w.r + w.zc
		txtw.WriteLine("1")
		txtw.WriteLine(fmt.Sprintf("%v %v", real(zCtrlGlobal), imag(zCtrlGlobal)))
	}
	txtw.Close()
}

// ExportComplexPotentialField creates a *.csv file containing the distribution of the resulting complex potential field for viewing
func (w *WatMethSoln) ExportComplexPotentialField(q *Prism, pointDensity int) {
	yn, yx, xn, xx := q.getExtentsXY()
	txtw, _ := mmio.NewTXTwriter("surfer-cell.bln")
	txtw.WriteLine(fmt.Sprint(len(q.Z) + 1))
	for _, v := range q.Z {
		txtw.WriteLine(fmt.Sprintf("%v %v", real(v), imag(v)))
	}
	txtw.WriteLine(fmt.Sprintf("%v %v", real(q.Z[0]), imag(q.Z[0])))
	txtw.Close()

	csvw := mmio.NewCSVwriter("surfer-hs.csv")
	csvw.WriteHead("x,y,h,s")
	for i := 0; i < pointDensity; i++ {
		fy := float64(i) / float64(pointDensity-1)
		fy *= yx - yn
		fy += yn
		for j := 0; j < pointDensity; j++ {
			fx := float64(j) / float64(pointDensity-1)
			fx *= xx - xn
			fx += xn
			if q.ContainsXY(fx, fy) {
				zl := (complex(fx, fy) - w.zc) / w.r // local coordinate
				o := w.cmplxPotFlow(zl)
				if w.qv != 0. {
					o += w.cmplxPotVert(zl)
				}
				if w.qw != 0. {
					o += w.cmplxPotWell(zl, 0)
				}
				h, s := real(o), imag(o)
				csvw.WriteLine(fx, fy, h, s)
			}
		}
	}
	csvw.Close()
}
