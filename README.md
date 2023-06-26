# StaticsPy
Library of functions for structural beam statics applications.

## Basic Workflow
1. Create a beam
2. Add load data using static methods
3. Build and combine load effects
4. Display and export data


### Beam Class
Create a new beam instance:

`
beam = Beam(length)
`

### Simple_Point
Static class that provides `Beam` classes with simply supported, concentrated load shear and moment information. Add individual loads on a beam object but passing it as the first argument:

`
Simple_Point(beam, P, a)
`

Follow the `Beam` object argument with the load magnitude `P` and the distance `a` from the left support.

Positive `P` acts downward.

## Example Outputs
`Beam.Show_All()` writes .png files to the ~/Documents/Statics folder. Example shear and moment output images are below. Note the combination load effect is a bold blue line, while load effects from point loads are shown as pink, and distributed loads as lightblue.
<img src="https://github.com/benstanfish/StaticsPy/assets/34006582/4eb75649-dfd3-412c-98c5-b76d7f0c706a" height="400" width="600" alt="Shear Diagram">
<img src="https://github.com/benstanfish/StaticsPy/assets/34006582/f9ff2f1c-db0b-4f49-b74e-d1faf85ad5f1" height="400" width="600" alt="Moment Diagram">
