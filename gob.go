package ptrack

import (
	"encoding/gob"
	"os"
)

func ExportPathlinesGob(fp string, pl [][]Particle, xr []int) error {
	mp := make(map[int][]Particle, len(pl))
	for i, pid := range xr {
		mp[pid] = pl[i]
	}
	f, err := os.Create(fp)
	if err != nil {
		return err
	}
	enc := gob.NewEncoder(f)
	err = enc.Encode(mp) // NOTE variable "toobig" altered from `const tooBig = (1 << 30) << (^uint(0) >> 62)`; to `const tooBig = (1 << 31) << (^uint(0) >> 62)` in decoder.go
	if err != nil {
		return err
	}
	f.Close()
	return nil
}

func LoadPathlinesGOB(fp string) (map[int][]Particle, error) {
	var d map[int][]Particle
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
