package ptrack

import "github.com/maseology/mmio"

const (
	reallyBig       = 9.999e100
	mingtzero       = 1e-8
	defaultPorosity = .3
	tol             = 1e-10 // "contains" tolerance
	wellTol         = .1
	rmax            = 1.5 // threshold of Contains()
)

func big(i int) string {
	return mmio.Thousands(int64(i))
}
