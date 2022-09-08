package ptrack

import (
	"encoding/gob"
	"os"
)

func ExportPathlinesGob(fp string, pl [][][]float64, xr []int) error {
	mp := make(map[int][][]float64, len(pl))
	for i, pid := range xr {
		mp[pid] = pl[i]
	}
	f, err := os.Create(fp)
	if err != nil {
		return err
	}
	enc := gob.NewEncoder(f)
	err = enc.Encode(mp)
	if err != nil {
		return err
	}
	f.Close()
	return nil
}

func LoadPathlinesGOB(fp string) (map[int][][]float64, error) {
	var d map[int][][]float64
	f, err := os.Open(fp)
	if err != nil {
		return nil, err
	}
	enc := gob.NewDecoder(f)
	err = enc.Decode(&d)
	if err != nil {
		return nil, err
	}
	f.Close()
	return d, nil
}
