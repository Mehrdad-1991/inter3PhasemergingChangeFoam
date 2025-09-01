# inter3PhasemergingChangeFoam
OpenFoam Solver

About the Solver (August 2025)
inter3PhasemergingChangeFoam is a customized three-phase solver developed for simulating incompressible, phase-changing flows involving three immiscible fluids. It is built upon a dedicated library structure under phaseChangeThreePhaseMixtures, which includes:

incompressibleThreePhaseMixture

ThreePhaseInterfaceProperties

phaseChangeThreePhaseMixtures

This solver is intended for advanced simulations such as cavitation, droplet breakup, and other complex multiphase phase-change phenomena.

🚀 New Feature: Dynamic Mesh Support
Unlike the original version of inter3PhasemergingChangeFoam, this enhanced version now supports dynamic mesh capabilities. This addition allows the solver to handle mesh motion and deformation, enabling more realistic modeling of physical processes such as moving boundaries, piston motions, or deformable geometries—something not possible in the original static-mesh-based implementation.

🧪 Development & Institutional Context
This updated solver was developed by M.Sc. Mehrdad Kazemi under the supervision of Prof. Dr.-Ing. habil. Nikolai Kornev at:

Universität Rostock
Fakultät für Maschinenbau und Schiffstechnik
Lehrstuhl für Modellierung und Simulation
Albert-Einstein-Str. 2
18059 Rostock, Deutschland

The solver and its underlying methodology have been submitted to the OpenFOAM® Journal, accompanied by a detailed research paper.

📚 Background & Lineage
This solver extends the work originally presented in:

"Implementation of a VoF solver with phase change for the simulation of internal cavitation and droplet breakup in injectors"
OpenFOAM® Journal, Vol. 4, 2024
DOI: 10.51560/ofj.v4.111 Authored by: Bjørn Christian Dueholm, Jesper de Claville Christiansen, Benny Endelt, Nikolaj Kristensen, Jakob Hærvig

Their solver was based on the standard interPhaseChangeFoam provided with OpenFOAM.

🛠️ Compilation Instructions

Source OpenFOAM v2406.

Run the Allwmake script. It compiles both the custom library and the solver using wmake.

Use Allwclean to clean all build artifacts.

⚠️ Compatibility Notice:
The solver has been tested on standard desktop PCs and HPC clusters running Linux Mint. Stability and compatibility are confirmed in these environments.

🚀 Running the Solver
After successful compilation, the solver can be launched via:
>>inter3PhasemergingChangeFoam

Happy Foaming!
