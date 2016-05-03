=============
 Formulations
=============

Passive scalar transport analsys should start with the case where reaction involves only fuel, oxidizer, and products:

.. math::

   \nu_{F}F + \nu_{O}O \rightleftharpoons \nu_{P}P

The mass fraction :math:`Y_{k}` of each species follows a species conservation equation that can be expressed by:

.. math::

   \frac{\partial \rho Y_{k}}{\partial t} + \frac{\partial}{\partial x_{i}}\left ( \rho u_{i} Y_{k} \right ) = \frac{\partial}{\partial x_{i}}\left ( \rho D \frac{\partial Y_{k}}{\partial x_{i}} \right ) + \dot{\omega}_{k}

Species production rates :math:`\dot{\omega}_{k}` are all related to each species because their production and consumption are correlated in a single step reaction. Note that this approach works out for the single step reaction assumption. Let :math:`Q` is the single-step reaction rate such that following relation can be introduced:

.. math::

   \dot{\omega}_{k} = W_{k} \nu_{k} Q

It enables to relates the oxidizer reaction rate to the fuel reaction rate given by:

.. math::

   \dot{\omega}_{O} = s \dot{\omega}_{F} \;\;\; \text{with} \;\;\; s = \frac{\nu_{O} W_{O}}{\nu_{F} W_{F}}

where :math:`s` is the mass stoichimetric ratio. And the reaction rate for temperature is also obviously linked to the fuel reaction rate:

.. math::

   \dot{\omega}_{T} = -Q \dot{\omega}_{F}

Using the relation stated above, three conservation equation for all three species in the single step reaction can be formulated:

.. math::

   \frac{\partial \rho Y_{F}}{\partial t} + \frac{\partial}{\partial x_{i}} \left ( \rho u_{i} Y_{F} \right ) = \frac{\partial }{\partial x_{i}} \left ( \rho D \frac{\partial Y_{F}}{\partial x_{i}} \right ) + \dot{\omega}_{F}

   \frac{\partial \rho Y_{O}}{\partial t} + \frac{\partial}{\partial x_{i}} \left ( \rho u_{i} Y_{O} \right ) = \frac{\partial }{\partial x_{i}} \left ( \rho D \frac{\partial Y_{O}}{\partial x_{i}} \right ) + s\dot{\omega}_{F}

   \frac{\partial \rho T}{\partial t} + \frac{\partial}{\partial x_{i}} \left ( \rho u_{i} T \right ) = \frac{\partial }{\partial x_{i}} \left ( \frac{\lambda}{C_{p}} \frac{\partial T}{\partial x_{i}} \right ) - \frac{Q}{C_{p}} \dot{\omega}_{F}

As of now, we are still having three conservation equations in addition to continuity, momentum equations to solve the diffusion flame physics. Our goal is to reduce the number of balance equation by employing a new variable. To obtain this goal, we can simplify the given equation by combining above three quations two by two, assuming unity Lewis numbers (:math:`Le=\lambda / \rho C_{p} D = 1`). This can achieved by linking the dependent variables, :math:`Y_{F}`, :math:`Y_{O}`, and :math:`T`:

.. math::

   Z_{1} = sY_{F} - Y_{O}\;\; ; \;\;\;\; Z_{2} = \frac{C_{p}T}{Q} + Y_{F} \;\; ; \;\;\;\; Z_{3} = s \frac{C_{p}T}{Q} + Y_{O}

These three quantities follow the same balance equation dropping the source terms:

.. math::

   \frac{\partial \rho Z}{\partial t} + \frac{\partial}{\partial x_{i}}\left ( \rho u_{i} Z \right ) = \frac{\partial}{\partial x_{i}} \left ( \rho D \frac{\partial Z}{\partial x_{i}} \right )

Here :math:`Z` is now introduced as a passive (or called conserved) scalar and changes only due to diffusion and convection. Without source term, this quantity is supposed to vary monotonically with properly setup of boundary condition. Now the three variables :math:`Z_{1}`, :math:`Z_{2}` and :math:`Z_{3}` follow the same balance equation stated above but have different boundary condition. It means these quantities stil have to be resolved in their own transport equation. To finally reduce the set of equations, the quantities should be normalized such a way that they end up with same boundary condition as well as same balance equation. The normalized :math:`z_{j}` variables can be achieved by defining:

.. math::

   Z_{1} = sY_{F} - Y_{O}\;\; ; \;\;\;\; Z_{2} = \frac{C_{p}T}{Q} + Y_{F} \;\; ; \;\;\;\; Z_{3} = s \frac{C_{p}T}{Q} + Y_{O}

These three quantities follow the same balance equation dropping the source terms:

.. math::

   \frac{\partial \rho Z}{\partial t} + \frac{\partial}{\partial x_{i}}\left ( \rho u_{i} Z \right ) = \frac{\partial}{\partial x_{i}} \left ( \rho D \frac{\partial Z}{\partial x_{i}} \right )

Here :math:`Z` is now introduced as a passive (or called conserved) scalar and changes only due to diffusion and convection. Without source term, this quantity is supposed to vary monotonically with properly setup of boundary condition. Now the three variables :math:`Z_{1}`, :math:`Z_{2}` and :math:`Z_{3}` follow the same balance equation stated above but have different boundary condition. It means these quantities stil have to be resolved in their own transport equation. To finally reduce the set of equations, the quantities should be normalized such a way that they end up with same boundary condition as well as same balance equation. The normalized :math:`z_{j}` variables can be achieved by defining:

.. math::

   z_{j} = \frac{Z_{j}-Z_{j}^{O}}{Z_{j}^{F}-Z_{j}^{O}} \;\;\;\; \text{for}\;\; j = 1,2,3

Now all reduced variables :math:`z_{j}` follow the same transport equation expressed by:

.. math::

   \frac{\partial \rho z_{j}}{\partial t} + \frac{\partial }{\partial x_{i}}\left ( \rho u_{i} z_{j} \right ) = \frac{\partial}{\partial x_{i}}\left ( \rho D \frac{\partial z_{j}}{\partial x_{i}} \right )

and have the same boundary conditions: math:`z_{j}` = 1 in the fuel stream and :math:`z_{j}` = 0 in the oxidizer stream. And this :math:`z` is now called the mixture fraction and measures the local fuel/oxidizer ratio. Again, this mixture fraction can also be converted to desired thermodynamic properties by formulating the following relations:

.. math::

   z = \frac{s Y_{F}-Y_{O}+Y_{O}^{0}}{sY_{F}^{0}+Y_{O}^{0}} = \frac{\frac{C_{p}}{Q}(T - T_{O}^{0}) + Y_{F}}{\frac{C_{p}}{Q} (T_{F}^{0} - T_{O}^{0} + Y_{F}^{0})} = \frac{\frac{sC_{p}}{Q} (T - T_{O}^{0}) + Y_{O} - Y_{O}^{0}}{\frac{sC_{p}}{Q}(T_{F}^{0} - T_{O}^{0}) - Y_{O}^{0}}


Now we came up with the finalized form of reduced balance equation without source term such that is is easily to solve it numerically and analytically sometime. However, it is not sufficient to resolve the reacting field. The given form of equation only provides the mixing properties with no combustion. In order to simulate the reacting field, we need an additional assumption: infinitely fast chemistry.

In this assumption, the chemical kinetics goes much faster than all flow process such as mixing. Since our approach is based on the single step reaction, there should be no chance that fuel and oxidizer co-exist. It meands that there will be a pure 'fuel' side with :math:`Y_{F} = 0` and a pure 'oxidizer' side with :math:`Y_{O} = 0`. Based on this assumption, the solution can classified into two separate zone:

- Fuel side (:math:`z > z_{st}`):

  .. math::

     Y_{F}(z) = zY_{F}^{0} + (z-1)\frac{Y_{O}^{0}}{s} = Y_{F}^{0}\frac{z-z_{st}}{1-z_{st}}

     Y_{O}(z) = 0

     T(z) = zT_{F}^{0} + (1-z)T_{O}^{0} + \frac{QY_{F}^{0}}{C_{p}}z_{st}\frac{1-z}{1-z_{st}}

- Oxidizer side (:math:`z < z_{st}`):

  .. math::

     Y_{F}(z) = 0

     Y_{O}(z) = Y_{O}^{0} \left ( 1-\frac{z}{z_{st}} \right )

     T(z) = zT_{F}^{0} + (1-z)T_{O}^{0} + \frac{QY_{F}^{0}}{C_{p}}z

where the flame position in the :math:`z`-space, :math:`z_{st}` is determined by expressing that the flame is located where both :math:`Y_{F}` and :math:`Y_{O}` are zero. At this location, :math:`z` is equal to its stoichiometric value given by:

.. math::

   z_{st} = \frac{1}{1+ \frac{sY_{F}^{0}}{Y_{O}^{0}}}
