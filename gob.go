package ptrack

import (
	"encoding/gob"
	"fmt"
	"os"

	"github.com/maseology/mmio"
)

func ExportPathlinesGob(fp string, pl [][]Particle) error {
	f, err := os.Create(fp)
	if err != nil {
		return err
	}
	enc := gob.NewEncoder(f)
	err = enc.Encode(pl) // NOTE variable "toobig" altered from `const tooBig = (1 << 30) << (^uint(0) >> 62)`; to `const tooBig = (1 << 31) << (^uint(0) >> 62)` in decoder.go
	if err != nil {
		return err
	}
	f.Close()
	return nil
}

func LoadPathlinesGOB(fp string) ([][]Particle, int64, error) {
	if s, ok := mmio.FileExists(fp); ok {
		var d [][]Particle
		f, err := os.Open(fp)
		if err != nil {
			return nil, -1, err
		}
		enc := gob.NewDecoder(f)
		err = enc.Decode(&d)
		if err != nil {
			return nil, -1, err
		}
		f.Close()
		return d, s, nil
	}
	return nil, -1, fmt.Errorf("file %s cannot be found", fp)
}
