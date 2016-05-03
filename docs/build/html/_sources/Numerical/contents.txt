==========================
 Numerical solution method
==========================


In this project, arbitrary compressibility method (ACM) is employed to solve the set of balance equations. The ACM is supposed to solve the continuity and momentum equation in a way of compressible solution method by introducing pseudo time integration in continuity equation. Unlike the way of typical density based solution, this method is to solve the pressure that is manipulated with arbitrary compressibility in the continuity equation insteady of density itself. So the algorithm is not proper for unsteady solution because of its pseudo time integration. However, it is very efficient to achieve steady solutin with a very less computational effort because it does not come up with higher value of characteristic velocity in the convective terms. Therefore the time step for this method does not have to reach up to the level of time scale of wave propagation. 

In this particular problem set, we need to add one additional balance equation with dependent variable of mixture fraction, :math:`z`. This can be easily achieve by simply adding the variable into the state variable vector. Like the same way of original ACM approach, the flux vector only needs to contain additional element corresponding to the convection and diffusion terms of passive scalar transport equation.

For this reason, the solution method is pretty much similar with the previous homework problem. The brief introduction of ACM approach is summarized below and onlly additional description of modified state vector and flux vectors are listed after that.

- Continuity (incompressible)

  .. math::

     \frac{\partial u_{i}}{\partial x_{i}} = 0 

- Momentum equation:

  .. math::

     \frac{\partial u_{i}}{\partial t} + \frac{\partial u_{i}u_{j}}{\partial x_{j}} = -\frac{1}{\rho}\frac{\partial p}{\partial x_{i}} + \nu \frac{\partial}{\partial x_{j}}\left ( \frac{\partial u_{i}}{\partial x_{j}} \right )


- Non-dimensionalization of the Navier-Stokes equations

  In some cases, it is beneficial to non-dimensionalize the given transport equation because it eases the analysis of problem of interest, and also may reduce the number of parameters. The non-dimensionalized form of the Navier-Stokes equation can be achieved by first normalizing the primitive variables as followings:

  .. math::

     \tilde{u_{i}} = \frac{u_{i}}{U_{\text{ref}}},\;\;  \tilde{x_{i}} = \frac{x_{i}}{L_{\text{ref}}},\;\; \tilde{\rho}=\frac{\rho}{\rho_{\text{ref}}},\;\;\tilde{P} = \frac{P}{\rho_{\text{ref}}\, U^{2}_{\text{ref}}},\;\; \tilde{t}=\frac{t}{L/U_{\text{ref}}}

  For the final form of non-dimensionalized Navier-Stokes equation, tilda, :math:`\tilde{}`, will be dropped out for brevity and a new non-dimensional physical parameter :math:`Re` that represents the flow intertia against the fluid viscosity is introduced. Now we got:

  .. math::

     \frac{\partial u_{i}}{\partial t} + \frac{\partial u_{i}u_{j}}{\partial x_{j}} = -\frac{1}{\rho}\frac{\partial p}{\partial x_{i}} + \frac{1}{\text{Re}} \frac{\partial}{\partial x_{j}}\left ( \frac{\partial u_{i}}{\partial x_{j}} \right )

  where the Reynolds number is defined as:

  .. math::

     \text{Re} = \frac{U_{\text{ref}}L_{\text{ref}}}{\nu}


- Artificial Compressiblity Method (ACM)

  In the artificial compressibility method (ACM), the continuity equation is modified adding an unsteady term with ariticial compressiblity :math:`\beta`. To have this new form of continuity equation, an artificial equation of state that relates pressure, :math:`P`, to artificial density :math:`\tilde{\rho}` is emploeyd as following form:

  .. math::

     P = \frac{\tilde{\rho}}{\beta}


  Finally, the modified continuity equation can then be recast as:

  .. math::

     \frac{\partial P}{\partial t} + \frac{1}{\beta} \frac{\partial u_{i}}{\partial x_{i}} = 0


- Addition of passive scalar transport equation

  In this project, the mixture fraction :math:`z` field needs to be added to the equation set. Looking at given transport equation of mixture fraction, it looks pretty much similar to the above equation set in terms of having convection and diffusion terms in common. As opposed to the derived equation in the previous section, we need to drop density quantity out of the equation because of incompressible flow assumption. As a result, the mixture fraction transport equation is rewritten in the form:

  .. math::

     \frac{\partial z}{\partial t} + \frac{\partial u_{j}z }{\partial x_{j}} = \frac{\partial }{\partial x_{j}}\left ( \text{Re} D \frac{\partial z}{\partial x_{j}} \right )


- Vector form of transport equations

  Rewriting the previously drived non-dimensionalized continuity and momentum equation in vector form generates a simple format that eases implementation of the numerical method. The above transport equation can be newly formed as shown below:

  .. math::

     \frac{\partial \vec{U}}{\partial t} + \frac{\partial \vec{E}}{\partial x} + \frac{\partial \vec{F}}{\partial y} = \frac{1}{\text{Re}} \left ( \frac{\partial^{2}}{\partial x^{2}} + \frac{\partial^{2}}{\partial y^{2}} \right ) \vec{U}

  where the each of vector elements are summarized below:

  .. math::

     \vec{U} = \begin{bmatrix}P\\ u\\ v \\ z \end{bmatrix}, \;\; \vec{E} = \begin{bmatrix} \frac{u}{\beta}\\ uu + P\\ uv \\ uz \end{bmatrix}, \;\; \vec{F} = \begin{bmatrix} \frac{v}{\beta}\\ uv\\ vv + P \\ vz \end{bmatrix}

  Now this is good to go further for descritization because the given task is to solve explicit form of discretization equation. Even though the derived form of transport equation is not linearized, each of vectors above are easily discretized in terms of their elements that are combinations of each primitive variables. Thus, in this project, actual discretization has been doen form the driven transport equation above.


- Finding time step algorithm

  In order to find time step that may stabilize the numerical solution, we need to know system convecting velocity as we pick the coefficient of spatial derivative terms in Burger's and Euler equations as the convection velocity. The driven system of equation is not a single partial different equation but a set of three different partial different equation. To find the convection speed of numerical information in the time and space domains, we need to first linearize the given system of equations and find the Eigen values. The linearization can be obatained by following process. The driven system of PDE should be reformulated in linearized set of equations:

  .. math::

     \frac{\partial \vec{U}}{\partial t}  + \left [ A \right ] \frac{\partial \vec{U}}{\partial x} + \left [ B \right ] \frac{\partial \vec{U}}{\partial y} = \frac{1}{\text{Re}} \left ( \frac{\partial^{2}}{\partial x^{2}} + \frac{\partial^{2}}{\partial y^{2}} \right ) \vec{U}

  Now we have found two coefficient matrices of convection terms and the spatial derivatives is now taken with respect to :math:`\vec{U}` only. Despite the vector form, the PDE form is a identical with Burger's equation. The coefficient matrices are below listed:

  .. math::

     \left [ A \right ] = \begin{bmatrix} 0 & \frac{1}{\beta} & 0 \\ 1 & 2u & 0\\ 0 & v & u \end{bmatrix}, \;\; \left [ B \right ] = \begin{bmatrix} 0 & 0 & \frac{1}{\beta} \\ 0 & v & u\\ 1 & 0 & 2v \end{bmatrix}

  The resolved Eigen values of :math:`\left [ A \right ]` and :math:`\left [ B \right ]` matrices are :math:`u, u+a, u-a` and :math:`v, v+a, v-a`, respectively. Taking :math:`\left [ A \right ]` for example, the maximum convection velocity that transmit the numerical information can then be :math:`\left | u  \right | + a`. Therefore, the Courant number for this case can also be determined by:

  .. math::

     C = \frac{|u_{x}| + a}{\Delta x} + \frac{|u_{y}| + a}{\Delta y}


