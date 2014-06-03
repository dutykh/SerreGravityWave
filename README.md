Serre's equations spectral solver
================

![Head-on collision](/pics/SerreHeadOn.jpg)

This Matlab script is a pseudo-spectral solver for the Serre-Green-Naghdi equations which model the propagation of long gravity waves. Here, for the sake of simplicity, we restrict our attention to the case of the flat bottom.

The numerical scheme is described in the following publication:

* **D. Dutykh**, D. Clamond, P. Milewski & D. Mitsotakis. [*Finite volume and pseudo-spectral schemes for the fully nonlinear 1D Serre equations*](http://hal.archives-ouvertes.fr/hal-00587994/), European Journal of Applied Mathematics, **24**(5), 761-787, 2013

**Remark**  
: In this script we implemented a good time discretization scheme (Dormand-Prince 5(4)), but with a constant time step (for simplicity). However, the authors strongly recommend to use the adaptive time stepping procedures such as the PI-control strategy.

Any comments are welcome!

## Authors

* [Denys Dutykh](http://www.denys-dutykh.com/) ([CNRS](http://www.cnrs.fr/insmi/) - [LAMA](http://lama.univ-savoie.fr/index.php), [University of Savoie](http://www.univ-savoie.fr/))

* [Didier Clamond](http://math.unice.fr/~didierc/) ([University of Nice Sophia-Antipolis](http://unice.fr/), [LJAD](http://math.unice.fr/))