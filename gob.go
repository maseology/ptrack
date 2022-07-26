package ptrack

import (
	"encoding/gob"
	"os"
)

func ExportPathlinesGob(fp string, pl [][][]float64) error {
	f, err := os.Create(fp)
	if err != nil {
		return err
	}
	enc := gob.NewEncoder(f)
	err = enc.Encode(pl)
	if err != nil {
		return err
	}
	f.Close()
	return nil
}

func LoadPathlinesGOB(fp string) ([][][]float64, error) {
	var d [][][]float64
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
