package main

import (
	"fmt"
	"math/cmplx"

	"github.com/maseology/mmio"
	pt "github.com/maseology/ptrack"
)

func main() {
	m := 80
	n := 30

	//////////////////////////////////////////////////////////////
	// EXAMPLES from Muhammad Ramadhan, 2015 //
	// // ex. 4.1
	// zj := []complex128{-.5 - .5i, -.5 + .5i, 1 + .5i, 1 - .5i} // vertex coordinate (clockwise order)
	// ztop := 0.5
	// zbot := -0.5
	// zw := 0. + 0i
	// qj := []float64{13., -13., -2., 2.}
	// Qtop := 0.
	// Qbot := 0.
	// Qwell := 0.
	// prt := pt.Particle{I: 0, X: -0.25, Y: -0.5, Z: 0., T: 0.}

	// // ex. 4.2
	// zj := []complex128{-.5 - 1.0i, -.5 + .5i, 0.5 + .5i, 0.5 - 1.0i} // vertex coordinate (clockwise order)
	// ztop := 0.5
	// zbot := -0.5
	// zw := 0. + 0i
	// qj := []float64{-3., 0., 0., 3.}
	// Qtop := 0.
	// Qbot := 0.
	// Qwell := 0.
	// prt := pt.Particle{I: 0, X: 0.48, Y: -1.0, Z: 0., T: 0.}

	// ex. 4.3
	// zj := []complex128{-.5 - 0.5i, -.5 + .5i, 0.5 + .5i, 0.5 - .5i} // vertex coordinate (clockwise order)
	// ztop := 1.0
	// zbot := 0.0
	// zw := 0. + 0i
	// qj := []float64{5., -2., -3., 4.}
	// Qtop := -3.
	// Qbot := -7.
	// Qwell := 0.
	// prt := pt.Particle{I: 0, X: -0.5, Y: -0.35, Z: 1., T: 0.}

	// // ex. 4.3 (2D for testing)
	// zj := []complex128{-.5 - 0.5i, -.5 + .5i, 0.5 + .5i, 0.5 - .5i} // vertex coordinate (clockwise order)
	// ztop := 1.0
	// zbot := 0.0
	// zw := 0. + 0i
	// qj := []float64{4., -3., -4., 3.}
	// Qtop := 0.
	// Qbot := 0.
	// Qwell := 0.
	// prt := pt.Particle{I: 0, X: -0.5, Y: -0.35, Z: 1., T: 0.}

	// fig. 3.2
	zj := []complex128{253.0483 + 53.0072i, 251.4244 + 49.7577i, 247.8321 + 50.2979i, 247.2358 + 53.8814i, 250.4596 + 55.5558i} // vertex coordinate (clockwise order)
	ztop := 0.5
	zbot := -0.5
	zw := 250. + 52.5i                         // complex well coordinate
	qj := []float64{-25., -25., 25., 15., 20.} // time-averaged flow rates through cell sides (Q_i) [L3/T]
	// qj := []float64{7., -7., 4., 3., -7.} // time-averaged flow rates through cell sides (Q_i) [L3/T]
	Qvert := -30. //10. //
	Qtop := 0.
	Qbot := -Qvert // Qvert = Qtop-Qbot
	Qwell := -40.  //0.    //  negative for pumping, positive for injecting (Q_well)
	prt := pt.Particle{I: 0, X: 248.4, Y: 54.1, Z: 0., T: 0.}

	//////////////////////////////////////////////////////////////

	var q pt.Prism
	q.New(zj, ztop, zbot, ztop-zbot, 0., 0.3)

	var wm pt.WatMethSoln
	if Qwell == 0. {
		zw = cmplx.NaN()
	}
	wm.New(&q, qj, zw, Qtop, Qbot, Qwell, m, n)
	wm.ExportComplexPotentialField(&q, 50)

	pl := pt.RungeKutta{Dt: 0.0001, Ds: 0.0001, Adaptive: false}
	// pl := pt.EulerSpace{Ds: 0.001}
	// pl := pt.EulerTime{Dt: 0.0001}

	ec, pathline := pt.TrackToExit(&prt, &q, &wm, &pl, zw)
	// for i, v := range pathline {
	// 	fmt.Println(i, v)
	// }
	printPathline(pathline)
	plast := pathline[len(pathline)-1]
	printExitCode(ec, len(zj))
	fmt.Printf(" particle exit point: %6.4f %6.4f %6.4f %6.4f", plast[0], plast[1], plast[2], plast[3])
}

func printPathline(p [][]float64) {
	txtw := mmio.NewTXTwriter("pathline.bln")
	txtw.WriteLine(fmt.Sprint(len(p)))
	for _, v := range p {
		txtw.WriteLine(fmt.Sprintf("%v %v %v %v", v[0], v[1], v[2], v[3]))
	}
	txtw.Close()
}

func printExitCode(ec, nf int) {
	switch {
	case ec == -9999:
		fmt.Println(" error: particle has not appeared to exit prism")
	case ec < 0:
		fmt.Printf(" particle has exited at well %d\n", -(ec + 1))
	case ec == nf:
		fmt.Println(" particle has exited at bottom face")
	case ec == nf+1:
		fmt.Println(" particle has exited at top face")
	default:
		fmt.Printf(" particle has exited at face %d\n", ec)
	}
}
