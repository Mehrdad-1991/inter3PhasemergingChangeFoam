# TestCase
Test Case
🧪 Development & Institutional Context
This updated test case was developed by Mehrdad Kazemi under the supervision of Prof. Dr.-Ing. habil. Nikolai Kornev at:

Universität Rostock
Fakultät für Maschinenbau und Schiffstechnik
Lehrstuhl für Modellierung und Simulation
Mehrdad Kazemi
Albert-Einstein-Str. 2
18059 Rostock, Deutschland

Run

    Set options in Allrun:

    PARALLEL_RUN=true   # or false
    nSubDomains=10      # number of cores


    Execute:

    ./Allrun


    Serial mode: runs solver directly.

    Parallel mode: copies mesh, runs decomposePar, solver with MPI, and reconstructPar.

Clean

    To remove results and reset the case:

    ./Allclean


    Deletes 0/, processor*/, postProcessing/, and restores mesh from mesh/.

Notes

    Make sure mesh/polyMesh/ exists before running.

    Adjust solver and settings in system/ and constant/ as needed.
