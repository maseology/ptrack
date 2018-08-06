# ptrack.go

A generalized particle tracking and flux-field analysis suite, designed for unstructured model structures and computational scalability.

### beta-version 0.1 includes:

- the Waterloo Method: a semi-analytic groundwater particle tracking algorithm (Muhammad Ramadhan, 2015)
- Euler, Runge-Kutta pathline integration schemes with adaptive time-stepping

### work in progress:

- currently able to compute pathlines for one prism at a time; intending to develop a standard unstructured grid format for large-scale particle tracking
- numerical model output conversion routines (i.e., convert to the "standard unstructured grid format")
- output results to *.vtk for 3D visualizations and animations
- FlowSource volumetric flow tracking (Foley and Black, 2017)
- additional processing:
    - Monte Carlo uncertainty analysis
    - particle endpoint cluster analysis

*more details to come..*

## References

Foley and Black, 2017. Efficiently delineating volumetric capture areas and flow pathways using directed acyclic graphs and MODFLOW-description of the algorithms within FlowSource.

Muhammad Ramadhan, 2015. A Semi-Analytic Particle Tracking Algorithm for Arbitrary Unstructured Grids. MASc thesis. University of Waterloo.