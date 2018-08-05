package ptrack

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/maseology/mmaths"
	"github.com/maseology/mmio"
	plt "github.com/maseology/mmplt"
	"gonum.org/v1/gonum/mat"
)

// CoefTaylor : returns the Taylor coefficients solved for a given cell geometry
func CoefTaylor(zj []complex128, qj []float64, zw complex128, qVert, qWell float64, m, n int) ([]complex128, float64) {
	// The Waterloo method
	// Muhammad Ramadhan, 2015. a Semi-Analytic Particle Tracking Algorithm for Arbitrary Unstructured Grids. M.A.Sc thesis. University of Waterloo.

	if m < 2*n {
		panic("total control points is less than twice the order of approximation, no solution can be found")
	}

	nfaces := len(zj)
	var zc complex128 // cell centroid
	for j := 0; j < nfaces; j++ {
		zc += zj[j]
	}
	zc /= complex(float64(nfaces), 0.)

	// step 1: check mass balance (eq3.7)
	// step 2: compute cell face angles, lengths, (planform) area, perimeter, maximum radius
	alpha, lj, a, p, r, wbal := make([]float64, nfaces), make([]float64, nfaces), 0., 0., 0., 0.
	for j := 0; j < nfaces; j++ {
		wbal += qj[j]
		jj := (j + 1) % nfaces
		d := zj[jj] - zj[j]
		lj[j] = cmplx.Abs(d)
		a += real(zj[j])*imag(zj[jj]) - real(zj[jj])*imag(zj[j])
		p += lj[j]
		r = math.Max(r, cmplx.Abs(zc-zj[j]))
		alpha[j] = cmplx.Phase(complex(real(d), -imag(d))) // table 3.6 (note error in original thesis)
	}
	wbal += qVert + qWell
	if wbal > mingtzero {
		panic("cell mass-balance error")
	}
	a /= -2. // negative used here because vertices are entered in clockwise order
	if a <= 0. {
		panic("vertex error, may be given in counter-clockwise order")
	}

	// step 3: build control points
	sCtrl, zCtrl, ijx := make([]float64, m), make([]complex128, m), make([]int, m)
	pPrev, pNext, j, sc := 0., lj[0], 0, p/float64(m)
	for i := 0; i < m; i++ {
		sCtrl[i] = (float64(i) + 0.5) * sc // eq. 3.6a
		if sCtrl[i] > pNext {
			pPrev = pNext
			j++
			pNext += lj[j]
		}
		ijx[i] = j // control point to face cross-reference
		jj := (j + 1) % nfaces
		zCtrlGlobal := complex((sCtrl[i]-pPrev)/lj[j], 0.)*(zj[jj]-zj[j]) + zj[j] // eq. 3.6b
		zCtrl[i] = (zCtrlGlobal - zc) / complex(r, 0.)                            // local control point (Z_ctrl) eq. 3.6c
	}
	saveControlPoints(zCtrl, zc, r, m)

	// step 4: determine normalized cell flows: Qcell, qVert and qWell, and build Qtaylor
	zw = (zw - zc) / complex(r, 0.) // local well coordinate
	qTaylorni := make([]float64, m)
	for i := 0; i < m; i++ {
		j := ijx[i]
		ca := complex(math.Cos(alpha[j]), math.Sin(alpha[j]))
		qCellni := qj[j] / lj[j]                                                 // eq. 3.8
		qVertni := real(complex(-qVert*r*real(zCtrl[i])/a, 0.) * ca)             // eq. 3.11a
		qWellni := real(-ca * complex(qWell/2./math.Pi/r, 0.) / (zCtrl[i] - zw)) // eq. 3.11b
		qTaylorni[i] = qCellni - qVertni - qWellni                               // eq. 3.10
	}

	// step 5: build Phi_Taylor
	psiTaylor := make([]float64, m)
	psiTaylor[0] = qTaylorni[0] * sCtrl[0]
	for i := 1; i < m; i++ {
		psiTaylor[i] = 0.5*(qTaylorni[i]+qTaylorni[i-1])*(sCtrl[i]-sCtrl[i-1]) + psiTaylor[i-1] // eq.3.12b
	}

	// step 6: build matrices to solve [phiU]*[a]=[phyTaylor] (eq.3.15)
	psiTaylorU := make([]float64, 2*m*n)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			jf := complex(float64(j), 0.)
			c := cmplx.Pow(zCtrl[i], jf)
			psiTaylorU[i*2*n+j] = imag(c)
			psiTaylorU[i*2*n+j+n] = real(c)
		}
	}
	psiU := mat.NewDense(m, 2*n, psiTaylorU)
	psiT := mat.NewVecDense(m, psiTaylor)
	aTaylorV := svdSolve(psiU, psiT)

	aTaylor := make([]complex128, n)
	for j := 0; j < n; j++ {
		aTaylor[j] = complex(aTaylorV.At(j, 0), aTaylorV.At(j+n, 0))
	}

	plotPerimeterFlux(zj, aTaylor, qj, lj, zc, p, r, nfaces, n, 500)

	return aTaylor, r
}

func svdSolve(a *mat.Dense, b *mat.VecDense) *mat.VecDense {
	// following https://www.youtube.com/watch?v=oTCLm-WnX9Y
	// svdSolve(mat.NewDense(3, 2, []float64{1., 0., 0., 2., 0., 1.}), mat.NewVecDense(3, []float64{0., 1., 0.}))
	// Solve x in Ax=b
	ar, ac := a.Dims()

	var svd mat.SVD
	if !svd.Factorize(a, mat.SVDFull) {
		panic("SVD solver error")
	}
	u := svd.UTo(nil)
	v := svd.VTo(nil)
	sv := svd.Values(nil) // sigma vectors
	for i := 0; i < len(sv); i++ {
		if sv[i] != 0. {
			sv[i] = 1. / sv[i]
		}
	}
	s := mat.NewDiagonalRect(ar, ac, sv)
	si := mat.DenseCopyOf(s.T()) // pseudo-inverse

	z := mat.NewVecDense(ar, nil)
	z.MulVec(u.T(), b)

	y := mat.NewVecDense(ac, nil)
	y.MulVec(si, z)

	x := mat.NewVecDense(ac, nil)
	x.MulVec(v, y)

	// fmt.Printf("\nA matrix: \n%2.1f\n\n",
	// 	mat.Formatted(a, mat.Squeeze()))
	// fmt.Printf("\nb vector: \n%2.1f\n\n",
	// 	mat.Formatted(b, mat.Squeeze()))

	// fmt.Printf("\nu matrix: \n%6.3f\n\n",
	// 	mat.Formatted(u, mat.Squeeze()))
	// fmt.Printf("\nv matrix: \n%6.3f\n\n",
	// 	mat.Formatted(v, mat.Squeeze()))
	// fmt.Printf("\ns matrix: \n%6.3f\n\n",
	// 	mat.Formatted(s, mat.Squeeze()))
	// fmt.Printf("\ns' matrix: \n%6.3f\n\n",
	// 	mat.Formatted(si, mat.Squeeze()))
	// fmt.Printf("\nz matrix: \n%6.3f\n\n",
	// 	mat.Formatted(z, mat.Squeeze()))
	// fmt.Printf("\ny matrix: \n%6.3f\n\n",
	// 	mat.Formatted(y, mat.Squeeze()))
	// fmt.Printf("\nx matrix: \n%6.3f\n\n",
	// 	mat.Formatted(x, mat.Squeeze()))

	return x
}

func plotPerimeterFlux(zj, aT []complex128, qj, lj []float64, zc complex128, p, r float64, nfaces, n, b int) {
	// print perimeter fluxes for testing
	sCtrl, qn, qi, hi, si := make([]float64, b), make([]float64, b), make([]float64, b), make([]float64, b), make([]float64, b)
	pPrev, pNext, j, sc, ijx := 0., lj[0], 0, p/float64(b), make([]int, b)
	for i := 0; i < b; i++ {
		sCtrl[i] = (float64(i) + 0.5) * sc
		if sCtrl[i] > pNext {
			pPrev = pNext
			j++
			pNext += lj[j]
		}
		ijx[i] = j
		qi[i] = qj[j] / lj[j]
		jj := (j + 1) % nfaces
		zGlobal := complex((sCtrl[i]-pPrev)/lj[j], 0.)*(zj[jj]-zj[j]) + zj[j] // eq. 3.6b
		zLocal := (zGlobal - zc) / complex(r, 0.)
		h, s := func(zl complex128) (float64, float64) {
			var omega complex128 // complex potential function (eq. 2.9a)
			for i := 0; i < n; i++ {
				fj := complex(float64(i), 0.)
				omega += aT[i] * cmplx.Pow(zl, fj)
			}
			return real(omega), imag(omega)
		}(zLocal)
		hi[i] = h
		si[i] = s
	}
	er := 0.
	for i := 0; i < b; i++ {
		ip := (i + 1) % b
		im := i - 1
		if im < 0 {
			im += b
		}
		qn[i] = (si[ip] - si[im]) / 2. / sc // qn=dPsi/ds
		er += math.Abs(qi[i] - qn[i])
	}
	er /= float64(b) * mmaths.SliceMax(qi)

	fmt.Printf("average flow error: %6.4f\n\n", er)

	m := make(map[string][]float64)
	m["potential"] = hi
	m["stream function"] = si
	m["normal flux"] = qn
	m["specified normal flux"] = qi

	plt.Line("perimeterFlux.png", sCtrl, m)
}

func saveControlPoints(zCtrl []complex128, zc complex128, r float64, m int) {
	txtw := mmio.NewTXTwriter("controlpoints.bln")
	for i := 0; i < m; i++ {
		zCtrlGlobal := zCtrl[i]*complex(r, 0.) + zc
		txtw.WriteLine("1")
		txtw.WriteLine(fmt.Sprintf("%v %v", real(zCtrlGlobal), imag(zCtrlGlobal)))
	}
	txtw.Close()
}

func ExportComplexPotentialField(zj, aT []complex128, r float64, d, n int) {
	yn, yx, xn, xx := extents(zj)
	var zc complex128 // cell centroid
	txtw := mmio.NewTXTwriter("cell.bln")
	txtw.WriteLine(fmt.Sprint(len(zj) + 1))
	for j := range zj {
		txtw.WriteLine(fmt.Sprintf("%v %v", real(zj[j]), imag(zj[j])))
		zc += zj[j]
	}
	txtw.WriteLine(fmt.Sprintf("%v %v", real(zj[0]), imag(zj[0])))
	txtw.Close()
	zc /= complex(float64(len(zj)), 0.)
	csvw := mmio.NewCSVwriter("hs.csv")
	csvw.WriteHead("x,y,h,s")
	for i := 0; i < d; i++ {
		fy := float64(i) / float64(d-1)
		fy *= yx - yn
		fy += yn
		for j := 0; j < d; j++ {
			fx := float64(j) / float64(d-1)
			fx *= xx - xn
			fx += xn
			if mmaths.PnPolyC(zj, complex(fx, fy)) {
				continue
			}
			h, s := func(zl complex128) (float64, float64) {
				var omega complex128 // complex potential function (eq. 2.9a)
				for i := 0; i < n; i++ {
					fj := complex(float64(i), 0.)
					omega += aT[i] * cmplx.Pow(zl, fj)
				}
				return real(omega), imag(omega)
			}((complex(fx, fy) - zc) / complex(r, 0.)) // local coordinate
			csvw.WriteLine(fx, fy, h, s)
		}
	}
	csvw.Close()
}

func extents(zj []complex128) (float64, float64, float64, float64) {
	yn, yx, xn, xx := math.MaxFloat64, -math.MaxFloat64, math.MaxFloat64, -math.MaxFloat64
	for _, v := range zj {
		yn = math.Min(yn, imag(v))
		yx = math.Max(yx, imag(v))
		xn = math.Min(xn, real(v))
		xx = math.Max(xx, real(v))
	}
	return yn, yx, xn, xx
}
