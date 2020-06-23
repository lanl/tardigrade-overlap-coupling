############################################################
# Micro-scale simulation                                   #
############################################################
# The test micro-scale simulation for the overlap coupling #
############################################################

[Mesh]
  file = /media/nathan/projects/overlapTestMeshes/microscale.e
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Modules]
  [./TensorMechanics]
    [./Master]
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
    value = 0.1*t
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 68
    poissons_ratio = 0.32
  [../]
  [./stress]
    type = ComputeFiniteStrainElasticStress
  [../]
  [./density]
    type = Density
    density = 2000.
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

[Executioner]
#  type = Steady
  type = Transient
  num_steps = 1
  dt        = 0.1
  solve_type = 'PJFNK'
#  solve_type = 'NEWTON'
#  nl_rel_tol = 1e-8
#  nl_abs_tol = 1e-8
#  nl_max_its = 100
  #Terms for debugging
#  petsc_options = '-ksp_monitor_true_residual -ksp_compute_singularvalues' 
#  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  nl_max_its = 20
  l_max_its  = 5
#  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
#  petsc_options_value = 'hypre    boomeramg      100'
#  petsc_options_iname = '-ksp_gmres_restart'
#  petsc_options_value = '100'
#  petsc_options = '-snes_ksp_ew -ksp_monitor_true_residual -ksp_compute_singularvalues'# -pc_svd_monitor'
#  petsc_options = '-ksp_monitor_true_residual -ksp_compute_singularvalues'# -pc_svd_monitor'
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart -print_linear_residuals'# -ksp_view_mat'
  petsc_options_value = 'asm      lu           1               101                false                  '# binary'
[]

[Outputs]
  exodus = true
  perf_graph = true
  [./xdmf]
    type = Xdmf
  [../]
[]
