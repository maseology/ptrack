package ptrack

// PrismFlux contains the flux data normal to prism faces
type PrismFlux struct {
	q  []float64 // left-up-right-down-bottom-top
	qw float64   // source/sink (at cell centroid??)
}

// LatBotTop returns prism fluxes organized by lateral, top and bottom fluxes
func (pf *PrismFlux) LatBotTop() ([]float64, float64, float64) {
	nf := len(pf.q)
	return pf.q[:nf-2], pf.q[nf-2], pf.q[nf-1]
}
