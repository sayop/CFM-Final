============
 Conclusions
============

In this computer project, a simplified passive scalar transport equation is introduced and added to the numerical simulation solution method. Typical diffusion flame simulation requires to set more than 3 balance equations in addition to the continuity and momentum equations. By employing the passive scalar concept, we can reduce the number of equations very effectively and take an advantage of less computational effort. Addition of the new transport equation is very simple because it only needs to be added to the state vector and flux vectors. Brief observations from this computer project can be summarized as follows:

- Introduction of passive scalar concept demands less computational time and resources.
- The simplified assumption of the thermodynamic and chemical kinetics enables us to reduce the number of equation set.
- In this project, a counter flow diffusion flame is simulated and gives a good practice for the diffusion flame exploration.
- The mixture fraction goes from 0 to 1 (from pure oxidizer to pure fuel) and the stoichiometric condition is made at lean mixture fraction position.
- Increasing Reynolds of inlet stream makes the flame thickness thinner and flow features unstable.
