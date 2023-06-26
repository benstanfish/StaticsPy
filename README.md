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

Follow the Beam object argument with the load magnitude P and the distance "a" from the left support.