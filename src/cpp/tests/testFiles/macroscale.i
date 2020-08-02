#####################################
#       Macroscale simulation       #
#####################################
# The definition of the macro-scale #
# deformation.                      #
#####################################

[Mesh]
  file = /media/nathan/projects/overlapTestMeshes/macroscale.e
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
  [./phi_xx]
  [../]
  [./phi_yy]
  [../]
  [./phi_zz]
  [../]
  [./phi_yz]
  [../]
  [./phi_xz]
  [../]
  [./phi_xy]
  [../]
  [./phi_zy]
  [../]
  [./phi_zx]
  [../]
  [./phi_yx]
  [../]
[]

[Kernels]
  #Define the internal force balance equations
  [./force_1]
    type = InternalForce
    component = 0
    dof_num   = 0
    variable  = disp_x
    save_in   = FInternal_1

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./force_2]
    type = InternalForce
    component = 1
    dof_num   = 1
    variable  = disp_y
    save_in   = FInternal_2

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./force_3]
    type = InternalForce
    component = 2
    dof_num   = 2
    variable  = disp_z
    save_in   = FInternal_3

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  #Define the internal couple balance equations
  [./couple_11]
    type = InternalCouple
    component_i = 0
    component_j = 0
    dof_num     = 3
    variable    = phi_xx
    save_in   = CInternal_11

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_12]
    type = InternalCouple
    component_i = 0
    component_j = 1
    dof_num     = 4
    variable    = phi_xy
    save_in   = CInternal_12

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_13]
    type = InternalCouple
    component_i = 0
    component_j = 2
    dof_num     = 5
    variable    = phi_xz
    save_in   = CInternal_13

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_21]
    type = InternalCouple
    component_i = 1
    component_j = 0
    dof_num     = 6
    variable    = phi_yx
    save_in   = CInternal_21

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_22]
    type = InternalCouple
    component_i = 1
    component_j = 1
    dof_num     = 7
    variable    = phi_yy
    save_in   = CInternal_22

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_23]
    type = InternalCouple
    component_i = 1
    component_j = 2
    dof_num     = 8
    variable    = phi_yz
    save_in   = CInternal_23

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_31]
    type = InternalCouple
    component_i = 2
    component_j = 0
    dof_num     = 9
    variable    = phi_zx
    save_in   = CInternal_31

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_32]
    type = InternalCouple
    component_i = 2
    component_j = 1
    dof_num     = 10
    variable    = phi_zy
    save_in   = CInternal_32

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_33]
    type = InternalCouple
    component_i = 2
    component_j = 2
    dof_num     = 11
    variable    = phi_zz
    save_in   = CInternal_33

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
[]


[AuxVariables]
  [./pk2_11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pk2_22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pk2_33]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sigma_11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sigma_22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sigma_33]
    order = CONSTANT
    family = MONOMIAL
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
  [./phiDot_11]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDot_12]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDot_13]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDot_21]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDot_22]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDot_23]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDot_31]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDot_32]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDot_33]
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
  [./phiDotDot_11]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDotDot_12]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDotDot_13]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDotDot_21]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDotDot_22]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDotDot_23]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDotDot_31]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDotDot_32]
    order = FIRST
    family = LAGRANGE
  [../]
  [./phiDotDot_33]
    order = FIRST
    family = LAGRANGE
  [../]

  [./FInternal_1]
    order = FIRST
    family = LAGRANGE
  [../]

  [./FInternal_2]
    order = FIRST
    family = LAGRANGE
  [../]

  [./FInternal_3]
    order = FIRST
    family = LAGRANGE
  [../]

  [./CInternal_11]
    order = FIRST
    family = LAGRANGE
  [../]

  [./CInternal_12]
    order = FIRST
    family = LAGRANGE
  [../]

  [./CInternal_13]
    order = FIRST
    family = LAGRANGE
  [../]

  [./CInternal_21]
    order = FIRST
    family = LAGRANGE
  [../]

  [./CInternal_22]
    order = FIRST
    family = LAGRANGE
  [../]

  [./CInternal_23]
    order = FIRST
    family = LAGRANGE
  [../]

  [./CInternal_31]
    order = FIRST
    family = LAGRANGE
  [../]

  [./CInternal_32]
    order = FIRST
    family = LAGRANGE
  [../]

  [./CInternal_33]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  [./pk2_11]
    type = MaterialStdVectorAux
    property = PK2
    index = 0
    variable = pk2_11
  [../]
  [./pk2_22]
    type = MaterialStdVectorAux
    property = PK2
    index = 4
    variable = pk2_22
  [../]
  [./pk2_33]
    type = MaterialStdVectorAux
    property = PK2
    index = 8
    variable = pk2_33
  [../]
  [./sigma_11]
    type = MaterialStdVectorAux
    property = SIGMA
    index = 0
    variable = sigma_11
  [../]
  [./sigma_22]
    type = MaterialStdVectorAux
    property = SIGMA
    index = 4
    variable = sigma_22
  [../]
  [./sigma_33]
    type = MaterialStdVectorAux
    property = SIGMA
    index = 8
    variable = sigma_33
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
  [./phiDot_11]
    type = DOFTimeDerivative
    derivative_order = 1
    coupled = phi_xx
    variable = phiDot_11
  [../]
  [./phiDot_12]
    type = DOFTimeDerivative
    derivative_order = 1
    coupled = phi_xy
    variable = phiDot_12
  [../]
  [./phiDot_13]
    type = DOFTimeDerivative
    derivative_order = 1
    coupled = phi_xz
    variable = phiDot_13
  [../]
  [./phiDot_21]
    type = DOFTimeDerivative
    derivative_order = 1
    coupled = phi_yx
    variable = phiDot_21
  [../]
  [./phiDot_22]
    type = DOFTimeDerivative
    derivative_order = 1
    coupled = phi_yy
    variable = phiDot_22
  [../]
  [./phiDot_23]
    type = DOFTimeDerivative
    derivative_order = 1
    coupled = phi_yz
    variable = phiDot_23
  [../]
  [./phiDot_31]
    type = DOFTimeDerivative
    derivative_order = 1
    coupled = phi_zx
    variable = phiDot_31
  [../]
  [./phiDot_32]
    type = DOFTimeDerivative
    derivative_order = 1
    coupled = phi_zy
    variable = phiDot_32
  [../]
  [./phiDot_33]
    type = DOFTimeDerivative
    derivative_order = 1
    coupled = phi_zz
    variable = phiDot_33
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
  [./phiDotDot_11]
    type = DOFTimeDerivative
    derivative_order = 2
    coupled = phi_xx
    variable = phiDotDot_11
  [../]
  [./phiDotDot_12]
    type = DOFTimeDerivative
    derivative_order = 2
    coupled = phi_xy
    variable = phiDotDot_12
  [../]
  [./phiDotDot_13]
    type = DOFTimeDerivative
    derivative_order = 2
    coupled = phi_xz
    variable = phiDotDot_13
  [../]
  [./phiDotDot_21]
    type = DOFTimeDerivative
    derivative_order = 2
    coupled = phi_yx
    variable = phiDotDot_21
  [../]
  [./phiDotDot_22]
    type = DOFTimeDerivative
    derivative_order = 2
    coupled = phi_yy
    variable = phiDotDot_22
  [../]
  [./phiDotDot_23]
    type = DOFTimeDerivative
    derivative_order = 2
    coupled = phi_yz
    variable = phiDotDot_23
  [../]
  [./phiDotDot_31]
    type = DOFTimeDerivative
    derivative_order = 2
    coupled = phi_zx
    variable = phiDotDot_31
  [../]
  [./phiDotDot_32]
    type = DOFTimeDerivative
    derivative_order = 2
    coupled = phi_zy
    variable = phiDotDot_32
  [../]
  [./phiDotDot_33]
    type = DOFTimeDerivative
    derivative_order = 2
    coupled = phi_zz
    variable = phiDotDot_33
  [../]
[]

[BCs]
  active = 'left_x back_z bottom_y'
#  active = 'left_x back_z bottom_y bottom_x top_y top_x'
  [./left_x]
    type = DirichletBC
    #type = PresetBC
    variable = disp_x
    boundary = 'left'
    #boundary = 'left right bottom top front back'
    preset = true
    value = 0
  [../]
  [./back_z]
    type = FunctionDirichletBC
    #type = PresetBC
    variable = disp_z
    boundary = 'back'
    #boundary = 'left right bottom top front back'
    preset = true
    function = back_bc
  [../]
  [./bottom_y]
    type = DirichletBC
    #type = PresetBC
    variable = disp_y
    boundary = 'bottom'
    #boundary = 'left right bottom top front back'
    preset = true
    value = 0
  [../]
[]

[Functions]
  [./top_bc]
    type  = ParsedFunction
    value = 0.1*t
  [../]
[]

[Functions]
  [./back_bc]
    type  = ParsedFunction
    value = -0.1*t
  [../]
[]

[Materials]
  [./linear_elastic]
    type = MicromorphicMaterial
    material_fparameters = '2 29.48e3 25.48e3 5 1e3 0.4e3 -1.5e3 -1.4e3 -3e3 11 0 0 0 0 0 0 10e5 0 0 0 0 2 .4e3 -3e3' 
    model_name = "LinearElasticity"

    #Coupled variables
    u1     = 'disp_x'
    u2     = 'disp_y'
    u3     = 'disp_z'
    phi_11 = 'phi_xx'
    phi_22 = 'phi_yy'
    phi_33 = 'phi_zz'
    phi_23 = 'phi_yz'
    phi_13 = 'phi_xz'
    phi_12 = 'phi_xy'
    phi_32 = 'phi_zy'
    phi_31 = 'phi_zx'
    phi_21 = 'phi_yx'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
#    type = FDP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 10
  dt        = 0.1
  solve_type = 'PJFNK'
#  solve_type = 'NEWTON'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  nl_max_its = 100
  [./TimeIntegrator]
    type = NewmarkBeta
    beta = 0.25
    gamma = 0.5
  [../]
[]

[Outputs]
  exodus = true
  perf_graph = true
  [./xdmf]
    type = Xdmf
  [../]
[]
