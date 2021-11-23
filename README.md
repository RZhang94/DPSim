# DPSim
RK4 based DPS simulators for MEng Capstone Project.

The MATLAB files are;
•RDSpatTimeEvoRK4wBSCont.m: A simulation of the diffusion reaction heat equationthat  utilizes  the  backstepping  approach  as  the  stabilizing  control  from  Section.3.1.Implementation of control effort has not been incorporated into each loop of the RK4simulation, this can be seen as the controller only only updating the control input onceper time increment.
•LinInelRK4v6LyapBasedBoundContrwTranslatablex.m:  A  simulation  of  the  heavyrope system with the the actuator end’s position or velocity being directly controllable.The controller implemented stabilizes the rope by moving the actuator end’s velocitydirectly, from Section.3.2.  The movement in the z direction has been locked.
•CascadeLinInelRK4v5Cascaded.m:  A simulation of the quadrotor-heavy rope sys-tem with the cascaded controller from Section.3.3.
