# julia-hackathon-hyderabad2025
Notebooks and exercises for the 2025 Julia Hackathon for IRIS-HEP at Hyderabad university


## Hackathon Projects
We are planning for two projects in the Hackathon sessions. The goal of both projects is to **produce new code that can be added to the [JuliaHEP](https://github.com/orgs/JuliaHEP/repositories) project**. In both cases the student will need to integrate the new code into an existing open source project providing all what is required for this contribution be accepted by the maintainers (tests, documentation, pull requests, etc.).

Depending on the level of knowledge of the Julia language, the student would be encourage to following the [Hands-on-Julia-for-particle-physicists](https://github.com/JuliaHEP/Hands-on-Julia-for-particle-physicists) tutorial. The tutorial has also been added to the Hackathon binder image. 

### [1] Data fitting related
 The idea of this project is to get access to some HEP experiment open data (or eventually link this project with the data produced by ePIC simulation), create a statistical and bayesian model for which we need to develop some new distribution that could be ultimately added to the [DistributionsHEP.jl](https://github.com/JuliaHEP/DistributionsHEP.jl) package. Specific files and information can be found in the `fitting` folder.

### [2] Detector simulation related
The idea of this project to develop a new detector simulation example using the [Geant4.jl](https://github.com/JuliaHEP/Geant4.jl) package. We can basically take an existing example from the Geant4 set of examples (basic or advanced) and convert it to Julia, which then can be added into the [G4Examples.jl](https://github.com/JuliaHEP/G4Examples.jl) package. 

Two examples in the `basic` the set of G4Examples are not yet implemented and are simple enough to be developed in the Hackathon:
- **B4:** Simplified calorimeter with layers of two materials. To demonstrate several possible ways of data scoring, the example is provided in four variants: B4a, B4b, B4c, B4d, with the same geometry but different ways to obtain the simulation data using scoring, user actions or sensitive detectors. 
- **B5:** A double-arm spectrometer with wire chambers, hodoscopes and calorimeters with a local constant magnetic field 

Other examples are also possible from the set of `extended` and `advanced` examples. Or even complete new ones based on the detector the student is currently working.

People interested in this project should start by following the [Geant4.jl-tutorial](https://github.com/peremato/Geant4.jl-tutorial), which introduces the Geant4.jl Julia package. The tutorial has been added to the Hackathon binder image.    
