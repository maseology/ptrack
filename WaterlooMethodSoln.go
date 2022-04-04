package ptrack

import (
	"fmt"
	"math"
	"math/cmplx"

	mmplt "github.com/maseology/mmPlot"
	"github.com/maseology/mmio"
	"gonum.org/v1/gonum/mat"
)

func svdSolve(a *mat.Dense, b *mat.VecDense) *mat.VecDense {
	// following https://www.youtube.com/watch?v=oTCLm-WnX9Y
	// svdSolve(mat.NewDense(3, 2, []float64{1., 0., 0., 2., 0., 1.}), mat.NewVecDense(3, []float64{0., 1., 0.}))
	// Solve x in Ax=b
	ar, ac := a.Dims()

	var svd mat.SVD
	if !svd.Factorize(a, mat.SVDFull) {
		panic("SVD solver error")
	}
	u, v := &mat.Dense{}, &mat.Dense{}
	svd.UTo(u)
	svd.VTo(v)
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

	return x
}

func (w *WatMethSoln) cmplxPotFlow(zl complex128) complex128 {
	var omega complex128 // complex flow field function (eq. 2.9a and 3.9)
	for i := 0; i < w.n; i++ {
		fj := complex(float64(i), 0.)
		omega += w.aT[i] * cmplx.Pow(zl, fj)
	}
	return omega
}

func (w *WatMethSoln) cmplxPotVert(zl complex128) complex128 {
	h := -w.qv / 2. // vertical flux function (eq. 3.9)
	// h *= real(cmplx.Pow(w.r*zl, 2+0i))
	h *= math.Pow(real(w.r*zl), 2.)
	return complex(h, 0.)
}

func (w *WatMethSoln) cmplxPotWell(zl complex128, wID int) complex128 {
	h := w.qw // well function (eq. 3.9)
	h *= math.Log(real(w.r) * cmplx.Abs(zl-w.zwl[wID]))
	return complex(h, 0.)
}

func (w *WatMethSoln) cmplxVelFlow(zl complex128) complex128 {
	var omega complex128 // complex flow field function (eq. 3.16)
	for i := 1; i < w.n; i++ {
		fj := complex(float64(i-1), 0.)
		omega += complex(float64(i), 0.) * w.aT[i] * cmplx.Pow(zl, fj)
	}
	// fmt.Printf("  omega: %15.5e     zl: %15.5e\n", omega, zl)
	return -omega / w.r
}

func (w *WatMethSoln) cmplxVelVert(zl complex128) complex128 {
	v := w.qv * real(w.r*zl) // vertical flux function (eq. 3.16)
	return complex(v, 0.)
}

func (w *WatMethSoln) cmplxVelWell(zl complex128, wID int) complex128 {
	q := complex(-w.qw, 0.) // well function (eq. 3.16)
	q /= w.r * (zl - w.zwl[wID])
	return q
}

func (w *WatMethSoln) plotPerimeterFlux(zj []complex128, qj, lj []float64, p float64, b int) {
	// print perimeter fluxes for testing (such as Figure 3.13 in Ramadhan, 2015)
	sCtrl, qn, qi, qjj := make([]float64, b), make([]float64, b), make([]float64, b), make([]float64, b)
	zobs, hi, si := make([]complex128, b), make([]float64, b), make([]float64, b)
	pPrev, pNext, j, sc, ijx := 0., lj[0], 0, p/float64(b), make([]int, b)
	for i := 0; i < b; i++ {
		sCtrl[i] = (float64(i) + 0.5) * sc // eq. 3.6a (modified)
		if sCtrl[i] > pNext {
			pPrev = pNext
			j++
			pNext += lj[j]
		}
		ijx[i] = j
		qjj[i] = qj[j]
		qi[i] = qj[j] / lj[j] // normalized cell flows eq. 3.8
		jnext := (j + 1) % w.nf
		zGlobal := complex((sCtrl[i]-pPrev)/lj[j], 0.)*(zj[jnext]-zj[j]) + zj[j] // eq. 3.6b
		zLocal := (zGlobal - w.zc) / w.r                                         // eq. 3.6c
		zobs[i] = zGlobal
		o := w.cmplxPotFlow(zLocal)
		// cpot[i] = o
		hi[i] = real(o) // potential
		si[i] = imag(o) // stream line
	}

	er := 0.
	for i := 0; i < b; i++ {
		ip := (i + 1) % b
		im := i - 1
		if im < 0 {
			im += b
		}
		// ds := 0.5 * (sCtrl[ip] - sCtrl[im]) //:= sc
		// fmt.Println(ijx[i], sCtrl[i], ip, i, im, ds, sc, si[ip], si[i], si[im])

		// qn[i] = imag((cpot[ip] - cpot[im]) / complex(2*sc, 0.)) // qn=dPsi/ds
		qn[i] = (si[ip] - si[im]) / 2. / sc // qn=dPsi/ds
		er += math.Abs(qi[i] - qn[i])
		qn[i] *= lj[ijx[i]] // de-normalize
	}
	er /= float64(b) * sliceMax(qi)

	fmt.Printf("average flow error: %6.4f ", er)

	m := make(map[string][]float64)
	m["potential"] = hi
	m["stream function"] = si
	m["normal flux"] = qn
	m["specified normal flux"] = qjj

	mmplt.Line("perimeterFlux.png", sCtrl, m, 8, 3)

	txtw, _ := mmio.NewTXTwriter("surfer-perimeterFlux-obs.bln")
	for i := 0; i < b; i++ {
		txtw.WriteLine("1")
		txtw.WriteLine(fmt.Sprintf("%v %v", real(zobs[i]), imag(zobs[i])))
	}
	txtw.Close()
}

func sliceMax(s []float64) float64 {
	x := -math.MaxFloat64
	for _, v := range s {
		x = math.Max(x, v)
	}
	return x
}
