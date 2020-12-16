package ptrack

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"io/ioutil"
	"log"
	"sort"
	"time"
)

// ExportVTKpathlines saves particle tracking results as a *.vtk file for visualization.
func ExportVTKpathlines(filepath string, pl [][][]float64) {
	// write to data buffer
	buf, endi, np := new(bytes.Buffer), binary.BigEndian, 0
	for _, a := range pl {
		np += len(a)
	}

	binary.Write(buf, endi, []byte("# vtk DataFile Version 3.0\n"))
	binary.Write(buf, endi, []byte(fmt.Sprintf("Pathline: %d vertices, %s\n", np, time.Now().Format("2006-01-02 15:04:05"))))
	binary.Write(buf, endi, []byte("BINARY\n"))
	binary.Write(buf, endi, []byte("DATASET UNSTRUCTURED_GRID\n"))

	binary.Write(buf, endi, []byte(fmt.Sprintf("POINTS %d float\n", np)))
	for _, a := range pl {
		for _, aa := range a {
			binary.Write(buf, endi, float32(aa[0]))
			binary.Write(buf, endi, float32(aa[1]))
			binary.Write(buf, endi, float32(aa[2]))
		}
	}

	binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELLS %d %d\n", len(pl), np+len(pl))))
	ii := 0
	for _, a := range pl {
		binary.Write(buf, endi, int32(len(a)))
		for i := range a {
			binary.Write(buf, endi, int32(ii+i))
		}
		ii += len(a)
	}

	binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELL_TYPES %d\n", len(pl))))
	for i := 0; i < len(pl); i++ {
		binary.Write(buf, endi, int32(4)) // VTK_POLY_LINE
	}

	// // // cell index
	// // binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELL_DATA %d\n", np)))
	// // binary.Write(buf, endi, []byte(fmt.Sprintf("SCALARS time int\n")))
	// // binary.Write(buf, endi, []byte(fmt.Sprintf("LOOKUP_TABLE default\n")))
	// // for k, a := range pl {
	// // 	for i := 0; i < len(a); i++ {
	// // 		binary.Write(buf, endi, k)
	// // 	}
	// // }
	// // cell index
	// binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELL_DATA %d\n", len(pl))))
	// binary.Write(buf, endi, []byte(fmt.Sprintf("SCALARS time int\n")))
	// binary.Write(buf, endi, []byte(fmt.Sprintf("LOOKUP_TABLE default\n")))
	// for i := 0; i < len(pl); i++ {
	// 	binary.Write(buf, endi, i)
	// }

	// write to file
	if err := ioutil.WriteFile(filepath, buf.Bytes(), 0644); err != nil {
		log.Fatalf("ioutil.WriteFile failed: %v", err)
	}
}

// // ExportVTKpathline saves particle tracking results as a *.vtk file for visualization.
// func ExportVTKpathline(filepath string, pl [][]float64) {
// 	// write to data buffer
// 	buf, endi, np := new(bytes.Buffer), binary.BigEndian, len(pl)

// 	binary.Write(buf, endi, []byte("# vtk DataFile Version 3.0\n"))
// 	binary.Write(buf, endi, []byte(fmt.Sprintf("Pathline: %d vertices, %s\n", np, time.Now().Format("2006-01-02 15:04:05"))))
// 	binary.Write(buf, endi, []byte("BINARY\n"))
// 	binary.Write(buf, endi, []byte("DATASET UNSTRUCTURED_GRID\n"))

// 	binary.Write(buf, endi, []byte(fmt.Sprintf("POINTS %d float\n", np)))
// 	for i := 0; i < np; i++ {
// 		binary.Write(buf, endi, float32(pl[i][0]))
// 		binary.Write(buf, endi, float32(pl[i][1]))
// 		binary.Write(buf, endi, float32(pl[i][2]))
// 	}

// 	binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELLS %d %d\n", 1, np+1)))
// 	binary.Write(buf, endi, int32(np))
// 	for i := 0; i < np; i++ {
// 		binary.Write(buf, endi, int32(i))
// 	}

// 	binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELL_TYPES %d\n", 1)))
// 	binary.Write(buf, endi, int32(4)) // VTK_POLY_LINE

// 	// // cell index
// 	// binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELL_DATA %d\n", np)))
// 	// binary.Write(buf, endi, []byte(fmt.Sprintf("SCALARS time float\n")))
// 	// binary.Write(buf, endi, []byte(fmt.Sprintf("LOOKUP_TABLE default\n")))
// 	// for i := 0; i < np; i++ {
// 	// 	binary.Write(buf, endi, float32(pl[i][3]))
// 	// }

// 	// write to file
// 	if err := ioutil.WriteFile(filepath, buf.Bytes(), 0644); err != nil {
// 		log.Fatalf("ioutil.WriteFile failed: %v", err)
// 	}
// }

// ExportVTK saves particle tracking results as a *.vtk file for visualization.
func (d *Domain) ExportVTK(filepath string) {
	// collect cell ids
	cids := make([]int, 0)
	for i := range d.prsms {
		cids = append(cids, i)

		var wm WatMethSoln
		ql, qb, qt := d.flx[i].LatBotTop()
		wm.New(d.prsms[i], ql, d.zc[i], -qt, qb, d.flx[i].qw, d.m, d.n)
		d.VF[i] = &wm
	}
	sort.Ints(cids)
	nprsm := len(cids)

	// collect vertices
	v, vxr, nvert := func() (map[int][]float64, map[int][]int, int) {
		v, vxr, cnt := make(map[int][]float64), make(map[int][]int), 0
		for _, i := range cids {
			p, s1 := d.prsms[i], make([]int, 0)
			for _, c := range p.Z {
				v[cnt] = []float64{real(c), imag(c), p.Top}
				s1 = append(s1, cnt)
				cnt++
			}
			for _, c := range p.Z {
				v[cnt] = []float64{real(c), imag(c), p.Bot}
				s1 = append(s1, cnt)
				cnt++
			}
			vxr[i] = vtkReorder(s1)
		}
		nvert := cnt
		// // collect centroids
		// for _, i := range cids {
		// 	p := d.prsms[i]
		// 	cxy := p.CentroidXY()
		// 	v[cnt] = []float64{real(cxy), imag(cxy), (p.Top + p.Bot) / 2.}
		// 	cnt++
		// }
		return v, vxr, nvert
	}()

	// write to data buffer
	buf, endi := new(bytes.Buffer), binary.BigEndian

	binary.Write(buf, endi, []byte("# vtk DataFile Version 3.0\n"))
	binary.Write(buf, endi, []byte(fmt.Sprintf("Unstructured prism domain: %d Prisms, %d vertices, %s\n", nprsm, nvert, time.Now().Format("2006-01-02 15:04:05"))))
	binary.Write(buf, endi, []byte("BINARY\n"))
	binary.Write(buf, endi, []byte("DATASET UNSTRUCTURED_GRID\n"))

	binary.Write(buf, endi, []byte(fmt.Sprintf("POINTS %d float\n", nvert)))
	for i := 0; i < nvert; i++ {
		binary.Write(buf, endi, float32(v[i][0]))
		binary.Write(buf, endi, float32(v[i][1]))
		binary.Write(buf, endi, float32(v[i][2]))
	}

	binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELLS %d %d\n", nprsm, nprsm+nvert)))
	for _, i := range cids {
		binary.Write(buf, endi, int32(len(d.prsms[i].Z)*2))
		for _, nid := range vxr[i] {
			binary.Write(buf, endi, int32(nid))
		}
	}
	// binary.Write(buf, endi, []byte(fmt.Sprintf("POINTS %d float\n", nprsm+nvert)))
	// for i := 0; i < nprsm+nvert; i++ {
	// 	binary.Write(buf, endi, float32(v[i][0]))
	// 	binary.Write(buf, endi, float32(v[i][1]))
	// 	binary.Write(buf, endi, float32(v[i][2]))
	// }

	// binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELLS %d %d\n", 2*nprsm, 3*nprsm+nvert)))
	// for _, i := range cids {
	// 	binary.Write(buf, endi, int32(len(d.prsms[i].Z)*2))
	// 	for _, nid := range vxr[i] {
	// 		binary.Write(buf, endi, int32(nid))
	// 	}
	// }
	// // centroids
	// for i := range cids {
	// 	binary.Write(buf, endi, int32(1))
	// 	binary.Write(buf, endi, int32(nvert+i))
	// }

	// binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELL_TYPES %d\n", 2*nprsm)))
	binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELL_TYPES %d\n", nprsm)))
	for _, i := range cids {
		switch len(d.prsms[i].Z) {
		case 0, 1, 2:
			log.Fatalf("ExportVTK error: invalid prism shape")
		case 3:
			binary.Write(buf, endi, int32(13)) // VTK_WEDGE
		case 4:
			binary.Write(buf, endi, int32(12)) // VTK_HEXAHEDRON
		case 5:
			binary.Write(buf, endi, int32(15)) // VTK_PENTAGONAL_PRISM
		case 6:
			binary.Write(buf, endi, int32(16)) // VTK_HEXAGONAL_PRISM
		default:
			log.Fatalf("ExportVTK todo: >6 sided polyhedron")
		}
	}
	// // centroids
	// for i := 0; i < nprsm; i++ {
	// 	binary.Write(buf, endi, int32(1)) // VTK_VERTEX
	// }

	// cell index
	// binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELL_DATA %d\n", 2*nprsm)))
	binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELL_DATA %d\n", nprsm)))
	binary.Write(buf, endi, []byte(fmt.Sprintf("SCALARS cellID int32\n")))
	binary.Write(buf, endi, []byte(fmt.Sprintf("LOOKUP_TABLE default\n")))
	for _, i := range cids {
		binary.Write(buf, endi, int32(i))
	}
	// for _, i := range cids {
	// 	binary.Write(buf, endi, int32(i))
	// }

	// centroid flux
	binary.Write(buf, endi, []byte(fmt.Sprintf("\nVECTORS cflux float\n")))
	// binary.Write(buf, endi, []byte(fmt.Sprintf("LOOKUP_TABLE default\n")))
	for _, i := range cids {
		// if vf, ok := d.VF[i]; ok {
		// 	// vv := v[k+nvert]
		// 	q := d.prsms[i]
		// 	cxy := q.CentroidXY()
		// 	p := Particle{I: 0, X: real(cxy), Y: imag(cxy), Z: (q.Top + q.Bot) / 2., T: 0.}
		// 	vx, vy, vz := vf.PointVelocity(&p, q, 0.)
		// 	fmt.Printf("  %.3f  %.3f  %.3f\n", vx, vy, vz)
		// 	binary.Write(buf, endi, float32(vx))
		// 	binary.Write(buf, endi, float32(vy))
		// 	binary.Write(buf, endi, float32(vz))
		// } else {
		// 	binary.Write(buf, endi, float32(0.))
		// 	binary.Write(buf, endi, float32(0.))
		// 	binary.Write(buf, endi, float32(0.))
		// }
		q := d.prsms[i]
		cxy := q.CentroidXY()
		p := Particle{I: 0, X: real(cxy), Y: imag(cxy), Z: (q.Top + q.Bot) / 2., T: 0.}
		fmt.Println(p)
		vx, vy, vz := d.VF[i].PointVelocity(&p, q, 0.)
		fmt.Printf("  %.3f  %.3f  %.3f\n", vx, vy, vz)
		binary.Write(buf, endi, float32(vx))
		binary.Write(buf, endi, float32(vy))
		binary.Write(buf, endi, float32(vz))
	}

	// write to file
	if err := ioutil.WriteFile(filepath, buf.Bytes(), 0644); err != nil {
		log.Fatalf("ioutil.WriteFile failed: %v", err)
	}
}

func vtkReorder(s []int) []int {
	l := len(s) / 2
	f := func(s []int) []int {
		for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
			s[i], s[j] = s[j], s[i]
		}
		return s
	}
	return append(f(s[:l]), f(s[l:])...)
}
