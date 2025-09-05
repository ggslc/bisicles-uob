main.poutBaseName = pout.ctrl_ideal_shelf_inv_btrc_tikhC@TIKHC_x
#domain details
main.domain_size = 128.0e+3 128.0e+3 1.0e+3
amr.num_cells = 256 256 10 
amr.is_periodic = 0 0 0

#bulk constitutive relation and parameters
main.constitutiveRelation = GlensLaw
main.rateFactor = constRate
constRate.epsSqr0 = 1.0e-12
constRate.A = 1.0e-17

#geometry,temperarture & basal friction input data
geometry.problem_type = Python
PythonIBC.module = ideal_shelf
PythonIBC.thicknessFunction = thickness_basic
PythonIBC.topographyFunction = topography	
PythonIBC.bc_lo = 1 1 # bc type at x = 0, y = 0
PythonIBC.bc_hi = 0 1 # bc type at x = L , y = W

#geometry.beta_type = Python
#PythonBasalFriction.module = ideal_shelf
#PythonBasalFriction.function = beta
#muCoefficient.type = Python
#PythonMuCoefficient.module = ideal_shelf
#PythonMuCoefficient.function = mucoef_rifted

temperature.type = constant
temperature.value = 268

amr.sigma = 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0

#basal friction relation parameters
main.basalFrictionRelation = powerLaw
BasalFrictionPowerLaw.m = 0.3333
BasalFrictionPowerLaw.includeEffectivePressure = false

#physical constants
constants.ice_density=918.0
constants.sea_water_density=1028.0

surfaceFlux.type = zeroFlux

basalFlux.type = maskedFlux
basalFlux.grounded.type = zeroFlux
basalFlux.floating.type = zeroFlux

#velocity solve type -- 0 = Picard, 1 = JFNK
amr.velocity_solver_type = 1

#initial guess 
amr.do_initial_velocity_guess = 1
amr.do_initial_velocity_solve = 1
amr.initial_velocity_guess_type = 1 # linear PDE
amr.initial_velocity_guess_const_mu = 1.0e+6
amr.initial_velocity_guess_solver_type = 1 
amr.initial_velocity_guess_const_vel = 1.0 0.0
#JFNK parameters
JFNKSolver.vtopSafety = 0.95
JFNKSolver.minPicardIterations = 8
JFNKSolver.maxIter = 32
JFNKSolver.absTol = 1.0e-10
JFNKSolver.relTol = 1.0e-10
JFNKSolver.solverType = 0 #Relax solver
JFNKSolver.solverType = 4 #PETSc
JFNKSolver.RelaxRelTol = 2.0e-4
JFNKSolver.RelaxHang = 0.975
JFNKSolver.maxRelaxIter = 20
JFNKSolver.normType = 0
JFNKSolver.verbosity = 5
JFNKSolver.vtopRelaxTol = 0.005
JFNKSolver.vtopRelaxMinIter = 8
JFNKSolver.numMGSmooth = 8
JFNKSolver.numMGIter = 1
JFNKSolver.h = 0.025
JFNKSolver.switchRate = 1.8
#JFNKSolver.writeResiduals = true

#time stepping
main.maxTime = 0.0  # maximum time to run to
main.maxStep = 0   # maximum number of steps to runx
amr.temporal_accuracy = 1
amr.cfl = 0.125 
amr.initial_cfl = 0.125 
amr.max_dt_grow_factor = 1.1 
amr.time_step_ticks = 0

#plotting options
amr.plot_interval = 1
amr.plot_prefix = plot.ideal_shelf_inv_btrc_tikhC@TIKHC.

amr.write_solver_rhs = 1
amr.write_viscous_tensor = true
amr.write_thickness_sources = true
amr.write_mask = true

#check points
amr.check_interval = 1000
amr.check_prefix = chk.ideal_shelf_inv_btrc_tikhC@TIKHC.

amr.check_overwrite = 0
amr.verbosity = 5

#AMR mesh options
amr.maxLevel = 10 # finest level allowed in simulation
amr.ref_ratio = 2 2 2 2 2 2 2 2 2 2 2
amr.regrid_interval = 1 # number of timesteps between regridding
amr.blockFactor = 4      # block factor used in grid generation
amr.fill_ratio = 0.85    # how efficient are the grids
amr.nestingRadius = 4    # proper nesting radius required
amr.tags_grow = 4        # amount to buffer tags
amr.tagCap = -1 # might increase this, but uniform grids are probably OK

amr.interpolate_zb = 0 #go back to the IBC to regrid geometry
amr.max_box_size = 64 # largest box length allowed

#AMR tagging
amr.tag_on_grad_velocity = 0

amr.tag_on_laplacian_velocity = 1
amr.lap_vel_tagging_val = 1.0

amr.tag_grounding_line = 1
amr.grounding_line_tagging_min_vel = 0.1

#misc options
amr.verbosity = 5
amr.velocity_solve_interval = 1

amr.velocity_solver_type = 6 # Inverse problem !

#initial guess for C!
geometry.beta_type = LevelData
inputLevelData.frictionFile = synthetic_observations_500m.2d.hdf5
inputLevelData.frictionName = beta_init

control.velMisfitType = speed # match observed speed (not velocity)
    
#field of data expressing confidence in speed observations. 0 for no confidence
control.velCoef.type = LevelData
control.velCoef.fileFormat = synthetic_observations_500m.2d.hdf5
control.velCoef.name = uc_grounded
control.velCoef.n = 1

#x-component of the observed velocity (could be the speed if the y-component is zero)
control.xVel.type = LevelData
control.xVel.fileFormat = synthetic_observations_500m.2d.hdf5
control.xVel.name = umod_dirty
control.xVel.n = 1

#y-component of the observed velocity (could be the speed if the x-component is zero))
#control.yVel.type = LevelData
#control.yVel.fileFormat = synthetic_observations_500m.2d.hdf5
#control.yVel.name = yVel
#control.yVel.n = 1

#optmize C/beta, not muCoef
control.optimizeX0 = true
control.optimizeX1 = false


#some box constraints. 
control.boundMethod = Projection
control.lowerX0 = -5
control.upperX0 = 5.0
control.lowerX1 = -3.0
control.upperX1 = 0.0125
    
control.writeInnerSteps = false
control.CGmaxIter = 1024
control.CGtol = 1.0e-300 # keep at it
control.CGhang = 2.0

control.minLevelForOptimization = 0 # refine to level 2 based on the observed speed, then optmize


control.outerStepFileNameBase = ctrl.ideal_shelf_inv_btrc_tikhC@TIKHC.

#inner steps
control.CGlineSearchMaxIter = 8
control.CGlineSearchInitialStep = 3.0e-8
control.CGlineSearchMethod = 20 # new secant
control.CGlineSearchTol = 1.0e-2

#control.rankRedBlack = true

#control.CGsaveStateInterval = -1
#control.CGsaveStateFileNameBase = cg-state.ideal_shelf_quick.

#Tikhonov regularization coefficients
control.gradCsqRegularization = @TIKHC
control.gradMuCoefsqRegularization = 0.0 # not relevant, but
control.gradX0sqRegularization = 1.0e-20
control.gradX1sqRegularization = 1.0e-20
control.X0sqRegularization = 1.0e-20
control.X1sqRegularization = 1.0e-20
