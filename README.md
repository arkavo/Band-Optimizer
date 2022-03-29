# Band-Optimizer
## Introduction
This is a set of scripts made to specifically plot band structures of substances(right now monoatomic support). Using Quantum-Espresso, we first optimize the .UPF files of the atoms given to obtain a set of optimized parameters. Using those parameters, we then create a band structure file and plot it directly using matplotlib. All the scripts are written in python for unix systems.

## Parameters
We will be varying 3 parameters in this calculation:
K points- The no of unit cells included in the calculation
Ecutwfc- The kinetic energy cutoff for the permissible wavefunctions
Lattice parameters- The dimensions of the said unit cell.

## Method
First, we will optimize K points from a suitable range 4-32 seems fine for the time being.
Note that there is obviously no "minimum" for this. A larger set of Unit cells will obviously give us more accurate and hence, lesser values of energy. What we instead seek is a k point that is sufficiently accurate and yet maintains some amount of speed.

Next, we optimize Ecutwfc in the same looping method with values from 10-120 as default.
This is the Kinetic energy upper cutoff for the wavefunction, a larger value will give us lower values of energy. An example diagram is given for ecut 200-400 as:![Tin ecut 200 400](https://user-images.githubusercontent.com/16555024/160578936-c1aa88aa-53a1-4234-8bc2-15a1fe22b943.png)
Our immediate goal is to choose a point that lies on the flats. That will give us a suitably accurate point for that accuracy domain.

Lastly, the lattice parameter optimization. This is the dimensions of the unit cell, there is no optimization behind this, the lattice constant with the lowest energy values will give us the most stable configuration. Since, this is the most sensitive parameter, we need to gauge this with some accuracy but keep speed.
For this, we divide this into Levels. Each "Level" will denote a decimal point, we will go short intervals on these levels till we get the desired accuracy.
For example, a value of 13.3658 will be obtained in 5 levels with a maximum of 5x10=50 calculations.
Level 0 gives us 13
Level 1 gives us 13.4
Level 2 gives us 13.36
Level 3 gives us 13.366
Level 4 gives us 13.3658
This is a good method to get speed.
