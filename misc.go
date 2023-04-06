package ptrack

import "strings"

func dirPrfx(fp string) string {
	return fp[:strings.Index(fp, "o.")]
}
