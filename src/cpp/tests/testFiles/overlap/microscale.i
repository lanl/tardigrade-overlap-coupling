############################################################
# Micro-scale simulation                                   #
############################################################
# The test micro-scale simulation for the overlap coupling #
############################################################

[Mesh]
  file = microscale_large.e#_small.e
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      save_in = 'FInternal_1 FInternal_2 FInternal_3'
      [./all]
        strain = FINITE
        add_variables = true
      [../]
    [../]
  [../]
[]

[BCs]
  active = 'left_x bottom_y front_z'
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
    preset = true
  [../]
  [./bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
    preset = true
  [../]
  [./front_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'front'
    preset = true
    function = front_bc
  [../]
[]

[Functions]
  [./front_bc]
    type  = ParsedFunction
    value = 1e-3*t
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 68e3
    poissons_ratio = 0.32
  [../]
  [./stress]
    type = ComputeFiniteStrainElasticStress
  [../]
  [./density]
    type = Density
    density = 1e-9
  [../]
[]

[Kernels]
  [./inertia_x]
    type = InertialForce
    variable = disp_x
    save_in = FInertial_1
    use_displaced_mesh = true
    density = density
  [../]
  [./inertia_y]
    type = InertialForce
    variable = disp_y
    save_in = FInertial_2
    use_displaced_mesh = true
    density = density
  [../]
  [./inertia_z]
    type = InertialForce
    variable = disp_z
    save_in = FInertial_3
    use_displaced_mesh = true
    density = density
  [../]
[]

#Stress outputs

[AuxVariables]
  [./s11]
    order = FIRST
    family = MONOMIAL
  [../]
  
  [./s12]
    order = FIRST
    family = MONOMIAL
  [../]
  
  [./s13]
    order = FIRST
    family = MONOMIAL
  [../]
  
  [./s21]
    order = FIRST
    family = MONOMIAL
  [../]
  
  [./s22]
    order = FIRST
    family = MONOMIAL
  [../]
  
  [./s23]
    order = FIRST
    family = MONOMIAL
  [../]
  
  [./s31]
    order = FIRST
    family = MONOMIAL
  [../]
  
  [./s32]
    order = FIRST
    family = MONOMIAL
  [../]
  
  [./s33]
    order = FIRST
    family = MONOMIAL
  [../]

  [./density]
    order = FIRST
    family = LAGRANGE
  [../]

  [./volume]
    order = FIRST
    family = LAGRANGE
  [../]

  [./uDot_1]
    order = FIRST
    family = LAGRANGE
  [../]

  [./uDot_2]
    order = FIRST
    family = LAGRANGE
  [../]

  [./uDot_3]
    order = FIRST
    family = LAGRANGE
  [../]

  [./uDotDot_1]
    order = FIRST
    family = LAGRANGE
  [../]

  [./uDotDot_2]
    order = FIRST
    family = LAGRANGE
  [../]

  [./uDotDot_3]
    order = FIRST
    family = LAGRANGE
  [../]

  [./FInternal_1 ]
    order = FIRST
    family = LAGRANGE
  [../]

  [./FInternal_2 ]
    order = FIRST
    family = LAGRANGE
  [../]

  [./FInternal_3 ]
    order = FIRST
    family = LAGRANGE
  [../]

  [./FInertial_1 ]
    order = FIRST
    family = LAGRANGE
  [../]

  [./FInertial_2 ]
    order = FIRST
    family = LAGRANGE
  [../]

  [./FInertial_3 ]
    order = FIRST
    family = LAGRANGE
  [../]

  [./FCoupling_1 ]
    order = FIRST
    family = LAGRANGE
  [../]

  [./FCoupling_2 ]
    order = FIRST
    family = LAGRANGE
  [../]

  [./FCoupling_3 ]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  [./s11]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 0
    j = 0
    variable = s11
    execute_on = "initial timestep_end"
  [../]
  [./s12]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 0
    j = 1
    variable = s12
    execute_on = "initial timestep_end"
  [../]
  [./s13]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 0
    j = 2
    variable = s13
    execute_on = "initial timestep_end"
  [../]
  [./s21]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 1
    j = 0
    variable = s21
    execute_on = "initial timestep_end"
  [../]
  [./s22]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 1
    j = 1
    variable = s22
    execute_on = "initial timestep_end"
  [../]
  [./s23]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 1
    j = 2
    variable = s23
    execute_on = "initial timestep_end"
  [../]
  [./s31]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 2
    j = 0
    variable = s31
    execute_on = "initial timestep_end"
  [../]
  [./s32]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 2
    j = 1
    variable = s32
    execute_on = "initial timestep_end"
  [../]
  [./s33]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 2
    j = 2
    variable = s33
    execute_on = "initial timestep_end"
  [../]
  [./density]
    type = NodalVolumeAverages
    variable = density
    execute_on = "initial timestep_end"
    ElementIntegrateUserObject = density_nodal_average
  [../]
  [./volume]
    type = NodalVolumeAverages
    variable = volume
    execute_on = "initial timestep_end"
    ElementIntegrateUserObject = density_nodal_average
    compute_nodal_volume = true
  [../]

  [./uDot_1]
    type = DOFTimeDerivative
    derivative_order = 1
    coupled = disp_x
    variable = uDot_1
  [../]

  [./uDot_2]
    type = DOFTimeDerivative
    derivative_order = 1
    coupled = disp_y
    variable = uDot_2
  [../]

  [./uDot_3]
    type = DOFTimeDerivative
    derivative_order = 1
    coupled = disp_z
    variable = uDot_3
  [../]

  [./uDotDot_1]
    type = DOFTimeDerivative
    derivative_order = 2
    coupled = disp_x
    variable = uDotDot_1
  [../]

  [./uDotDot_2]
    type = DOFTimeDerivative
    derivative_order = 2
    coupled = disp_y
    variable = uDotDot_2
  [../]

  [./uDotDot_3]
    type = DOFTimeDerivative
    derivative_order = 2
    coupled = disp_z
    variable = uDotDot_3
  [../]

[]

[UserObjects]
  [./density_nodal_average]
    type = ElementIntegrateUserObject
    variable = disp_x
    density = density
    execute_on = "INITIAL TIMESTEP_END"
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    #type = FDP
    full = true
  [../]
[]

[UserObjects]
    [./micro_coupling]
        type = OverlapCoupling
        is_macroscale = False
        overlap_configuration_filename = "testConfig.yaml"
        execute_on = "TIMESTEP_BEGIN"
    [../]
[]

[NodalKernels]
  [./bc_coupled_x]
    type = CouplingForce
    overlap_coupling_object = micro_coupling
    component = 0
    is_macroscale = False
    variable = disp_x
    save_in = FCoupling_1
  []
  [./bc_coupled_y]
    type = CouplingForce
    overlap_coupling_object = micro_coupling
    component = 1
    is_macroscale = False
    variable = disp_y
    save_in = FCoupling_2
  []
  [./bc_coupled_z]
    type = CouplingForce
    overlap_coupling_object = micro_coupling
    component = 2
    is_macroscale = False
    variable = disp_z
    save_in = FCoupling_3
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt        = 1e-2
  solve_type = 'PJFNK'
#  solve_type = 'NEWTON'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  nl_max_its = 50
  l_max_its  = 50
  [./TimeIntegrator]
    type = NewmarkBeta
    beta = 0.25
    gamma = 0.5
  [../]
[]

[Outputs]
  exodus = true
  perf_graph = true
  execute_on = "INITIAL TIMESTEP_END"
  [./xdmf]
    type = Xdmf
    execute_on = "INITIAL NONLINEAR"
  [../]
[]
