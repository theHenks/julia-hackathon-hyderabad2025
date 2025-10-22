# Getting Started Instructions

We are collecting here the steps to developing a new Geant4 example in julia.

1. Fork the repository https://github.com/JuliaHEP/G4Examples.jl
2. Clone you forked repository
3. Create a folder in `basic/Bx` (if you decided to implement B4 or B5)
4. Start by defining the geometry and materials of the detector. Look at the C++ implementation in https://github.com/Geant4/geant4/tree/master/examples
    - You can get inspiration from existing examples in G4Examples.jl
    - Build the geometry step by step using a notebook and display the partial construction with `draw(logcial_volume)`
6. Define the primary particles and physics list
7. Define the magnetic field
8. You can now define an `G4JLApplication` and run without collecting any information. To increase the tracking verbosity use the command ```ui`/tracking/verbose 1` ```
9. Depending whether the example is using `user actions`, or `sensitive volumes`, or `scoring` start implementing them.




