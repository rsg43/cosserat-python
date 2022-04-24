# cosserat-python
Python implemenation of the Cosserat rod model, based on work by Gazzola et al., translated from MATLAB code used for work carried out for my DPhil thesis. Includes a CosseratRod class to create filament objects, with parameters set to match the properties of actin filaments. This object stores the current state of the filament, and inlcudes a symplectic integrator method to allow numerical simulation. 

Only package dependency for main class is numpy, with remaining functions and operators included in class.

Completed Functionality updates:
- External forces
- Clamped filament segments

Planned Functionality updates:
- Visualisation code
- Dynamically clamped segments
- Cross-linked filaments
- Volume exclusion
