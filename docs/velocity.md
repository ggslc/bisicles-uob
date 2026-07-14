# BISICLES ice velocity

The ice velocity calculation depends at minimum on the [ice sheet
geometry](geometry.md) and [the in-ice and basal
stresses](stress.md). It may also depend on the
[thermodynamics](thermodynamics.md) and other factors that influence
the in-ice and basal stresses.

If you have built the doxygen code documentation, the C++ class
hierarchy underlying the ice velocity calculation is described at
[classIceVelocitySolver.md](../code/doc/doxygen/html/classIceVelocitySolver.html)

An ice velocity model is specified at run time by setting

      amr.velocity_solver_type = < integer solver_type, default 1 >

Additional options include

      amr.velocity_solver_interval = < integer n, default 1 > recompute the velocity every n timesteps.
      amr.regrid_interval = < integer m, default 1 > recompute the mesh and velocity every m timesteps.

Options for  < solver_type  > are

0.  [Picard Solver](#picard)  [deprecated, because the JFNKSolver
    provides the same functionality ]
1.  [JFNK Solver](#jfnk) (the default) solves the vertically integrated
    (SSA, L1L2) stress equations using a matrix-free Newton method.
2.  [Known Velocity](#known) specifies a known velocity
3.  [Petsc Nonlinear Solver](#snes)  [EXPERIMENTAL ] solves the
    vertically integrated (SSA, L1L2) stress equations using a PETSc
    SNES.
4.  [FASMGAMR](#fas)  [EXPERIMENTAL ] A full-approximation storage
    solver for the vertically integrated (SSA, L1L2) stress equations.
5.  [Python](#python) specifies velocity through the python interface
6.  [Inverse Vertically Integrated Velocity Solver](#inversevi)
    Estimates the basal friction and stiffness coefficients in order to
    provide a velocity field that fits observations.

## [JFNK Solver](#jfnk)

BISICLES JFNK (Jacobian-free Newton-Krylov) solver solves the nonlinear
vertically integrated (SSA, L1L2) stress equations using a matrix-free
Newton method. It is the default for conventional (forward) simulations,
and is also used by the [Inverse Vertically Integrated Velocity
Solver](#inversevi). BISICLES spends around 95% of its time solving the
stress equations, and it is also the calculation most likely to fail.
That said, it often fails or perform poorly when the input conditions -
basal traction, ice thickness, boundary conditions - are incorrect in
some way. For example, a region of floating ice unconnected to either
grounded ice, a rocky wall, or a domain boundary with Dirichlet
conditions does not have a unique solution - it is an ill-posed problem.
So if you have problems - especially on time step 0 - check the input
data first, especially the [geometry](geometry.md) and the [basal
traction](stress.md#basal).

However, problems that are in theory well-posed but ill-conditioned -
for example the Drygalski ice tongue - can see their iterations diverge,
and in those cases paying attention to the JFNKSolver options can be
useful. The most useful course in such cases is to try the [PETSc linear
solvers](#jfnk-petsc)

### [Option: PETSc linear solvers](#jfnk-petsc)

By default, the JFNK solver uses the Chombo geometric multigrid (GMG)
linear solvers to compute the Newton step. If you find that convergence
of the velocity solver is stalling (or even diverging), you might want
to try a PETSc algebraic multigrid (AMG) solver. To switch to the PETSc
solvers,

1.  You will need a .petsrc file to tell the PETSc solvers what to do.
    PETSc has an enormous range of options, but the important thing is
    to use a suitable preconditioner because the PETSc default does not
    work well in many cases. We usually find the
    [HYPRE](http://computation.llnl.gov/casc/linear_solvers/sls_hypre.md)
    BoomerAMG solver. to work well. A minimal .petscrc file would be
    something like

        #.petsrc
        #use the hypre / boomeramg preconditioner
        -pc_type hypre
        -pc_hypre_type boomeramg
        #gmres is the safe (and default) option for a KSP method. bcgs is sometimes quicker
        -ksp_type gmres
        #this makes the output less confusing, but it is optional
        -ksp_norm_type unpreconditioned
        #sensible values? depends what you are doing
        -ksp_max_it 20
        -ksp_rtol 1.0e-4
        #include this to see log the PETSc solver progress in pout.X etc
        -ksp_monitor
            

2.  Edit your BISICLES input file. There are two ways to use PETSC:

3.  As a bottom solver in the Chombo GMG solver. (if you don 't know
    what that means, don 't worry  -- suffice it to say that in many
    cases, this option occupies the middle ground between the Chombo
    native GMG solvers and a full AMR AMG solver.). Assuming you are
    using the JFNKSolver for your velocity solve, add the following line
    to your inputs file:

        JFNKSolver.bottom_solver_type = 1

4.  If that doesn 't help, then you 're ready for the full AMR PETSc
    solver. Assuming once again that you 're using the JFNKSolver, then
    change

         
          JFNKSolver.solverType = 0
          JFNKSolver.normType = 0

    to

         
          JFNKSolver.solverType = 4
          JNKKSolver.normType = 2

## [Known Velocity](##known)

## [Petsc Nonlinear Solver](#snes)

## [FASMGAMR](#fas)

## [Python](#pythin)

## [Inverse Vertically Integrated Velocity Solver](#inversevi)
