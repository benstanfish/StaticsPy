#StaticsPy
Library of functions for structural beam statics applications.

##Basic Workflow
Create a *Beam* object of specified length. Using the static methods of the various load types, e.g. *Simple_UDF* to add load parameters.

A beam can take several different loads, so when you're ready to run the loads use the beam object **Build_Loads()** then **Combine_Loads()** methods.

You can theb Show_Shear() or Show_Moment() to see the load effect diagrams.
