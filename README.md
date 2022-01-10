# ptrack.go

A generalized particle tracking and flux-field analysis suite, designed for unstructured model structures and computational scalability.

### current version includes:
- the Pollock (1989) method (only works for rectilinear model grids)
- the Waterloo method: a semi-analytic groundwater particle tracking algorithm (Muhammad Ramadhan, 2015)
- Euler and Runge-Kutta pathline integration schemes with adaptive time-stepping
- reads [MODFLOW6](https://www.usgs.gov/software/modflow-6-usgs-modular-hydrologic-model) output files
- output results to _*.vtk_ for 3D visualizations and animations

### work in progress:

- FlowSource volumetric flow tracking (Foley and Black, 2017)
- additional processing:
    - Monte Carlo uncertainty analysis
    - particle endpoint cluster analysis

*more details to come..*

## References

Foley and Black, 2017. Efficiently delineating volumetric capture areas and flow pathways using directed acyclic graphs and MODFLOW-description of the algorithms within FlowSource.

Pollock, D.W., 1989, Documentation of a computer program to compute and display pathlines using results from the U.S. Geological Survey modular three-dimensional finite-difference ground-water flow model: U.S. Geological Survey Open-File Report 89–381.

Muhammad Ramadhan, 2015. A Semi-Analytic Particle Tracking Algorithm for Arbitrary Unstructured Grids. MASc thesis. University of Waterloo.