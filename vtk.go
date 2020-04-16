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

// ExportVTKpathline saves particle tracking results as a *.vtk file for visualization.
func ExportVTKpathline(filepath string, pl [][]float64) {
	// write to data buffer
	buf, endi, np := new(bytes.Buffer), binary.BigEndian, len(pl)

	binary.Write(buf, endi, []byte("# vtk DataFile Version 3.0\n"))
	binary.Write(buf, endi, []byte(fmt.Sprintf("Pathline: %d vertices, %s\n", np, time.Now().Format("2006-01-02 15:04:05"))))
	binary.Write(buf, endi, []byte("BINARY\n"))
	binary.Write(buf, endi, []byte("DATASET UNSTRUCTURED_GRID\n"))

	binary.Write(buf, endi, []byte(fmt.Sprintf("POINTS %d float\n", np)))
	for i := 0; i < np; i++ {
		binary.Write(buf, endi, float32(pl[i][0]))
		binary.Write(buf, endi, float32(pl[i][1]))
		binary.Write(buf, endi, float32(pl[i][2]))
	}

	binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELLS %d %d\n", 1, np+1)))
	binary.Write(buf, endi, int32(np))
	for i := 0; i < np; i++ {
		binary.Write(buf, endi, int32(i))
	}

	binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELL_TYPES %d\n", 1)))
	binary.Write(buf, endi, int32(4)) // VTK_POLY_LINE

	// // cell index
	// binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELL_DATA %d\n", np)))
	// binary.Write(buf, endi, []byte(fmt.Sprintf("SCALARS time float\n")))
	// binary.Write(buf, endi, []byte(fmt.Sprintf("LOOKUP_TABLE default\n")))
	// for i := 0; i < np; i++ {
	// 	binary.Write(buf, endi, float32(pl[i][3]))
	// }

	// write to file
	if err := ioutil.WriteFile(filepath, buf.Bytes(), 0644); err != nil {
		log.Fatalf("ioutil.WriteFile failed: %v", err)
	}
}

// ExportVTK saves particle tracking results as a *.vtk file for visualization.
func ExportVTK(filepath string, d *Domain) {
	// collect cell ids
	cids := make([]int, 0)
	for k := range d.prsms {
		cids = append(cids, k)
	}
	sort.Ints(cids)

	// collect vertices
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

	// write to data buffer
	buf, endi := new(bytes.Buffer), binary.BigEndian

	binary.Write(buf, endi, []byte("# vtk DataFile Version 3.0\n"))
	binary.Write(buf, endi, []byte(fmt.Sprintf("Unstructured prism domain: %d Prisms, %d vertices, %s\n", len(d.prsms), len(v), time.Now().Format("2006-01-02 15:04:05"))))
	binary.Write(buf, endi, []byte("BINARY\n"))
	binary.Write(buf, endi, []byte("DATASET UNSTRUCTURED_GRID\n"))

	binary.Write(buf, endi, []byte(fmt.Sprintf("POINTS %d float\n", len(v))))
	for i := 0; i < cnt; i++ {
		binary.Write(buf, endi, float32(v[i][0]))
		binary.Write(buf, endi, float32(v[i][1]))
		binary.Write(buf, endi, float32(v[i][2]))
	}

	binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELLS %d %d\n", len(d.prsms), len(d.prsms)+cnt)))
	for _, i := range cids {
		binary.Write(buf, endi, int32(len(d.prsms[i].Z)*2))
		for _, nid := range vxr[i] {
			binary.Write(buf, endi, int32(nid))
		}
	}

	binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELL_TYPES %d\n", len(d.prsms))))
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

	// cell index
	binary.Write(buf, endi, []byte(fmt.Sprintf("\nCELL_DATA %d\n", len(d.prsms))))
	binary.Write(buf, endi, []byte(fmt.Sprintf("SCALARS cellID int32\n")))
	binary.Write(buf, endi, []byte(fmt.Sprintf("LOOKUP_TABLE default\n")))
	for _, i := range cids {
		binary.Write(buf, endi, int32(i))
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