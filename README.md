# ptrack.go

A generalized particle tracking and flux-field analysis suite, designed for unstructured model structures and computational scalability.

### beta-version 0.1 includes:

- the Waterloo Method: a semi-analytic groundwater particle tracking algorithm (Muhammad Ramadhan, 2015)

### work in progress:

- output to *.vtk
- FlowSource volumetric flow tracking (Foley and Black, 2017)
- model output conversion routines
- additional processing:
    - Monte Carlo uncertainty analysis
    - particle endpoint cluster analysis

*more details to come..*

## References

Foley and Black, 2017. Efficiently delineating volumetric capture areas and flow pathways using directed acyclic graphs and MODFLOW-description of the algorithms within FlowSource

Muhammad Ramadhan, 2015. A Semi-Analytic Particle Tracking Algorithm for Arbitrary Unstructured Grids. MASc thesis. University of Waterloo.