# cosserat-python
Python implemenation of the Cosserat rod model
Includes a CosseratRod class to create filament objects, with parameters set to match the properties of actin filaments. This object stores the current state of the filament, and inlcudes a symplectic integrator method to allow numerical simulation. 

Only package dependency is numpy, with remaining functions and operators included in class.

Planned Functionality updates:
- External forces
- Clamped filament segments
- Dynamically clamped segments
- Cross-linked filaments
