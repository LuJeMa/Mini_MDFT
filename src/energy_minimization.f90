subroutine energy_minimization

  use precision_kinds, only : dp
  use module_grid, only: grid
  use module_solvent, only: solvent
  use module_input, only: getinput
  use module_energy_and_gradient, only: energy_and_gradient
  !use module_lbfgs_nocedal_mdft, only: lbfgsb

  implicit none

  real(dp) :: f ! functional to minimize
  real(dp) :: df (grid%no, grid%nx, grid%ny, grid%nz, solvent(1)%nspec )
  real(dp) :: x_allsolv (grid%no, grid%nx, grid%ny, grid%nz, solvent(1)%nspec )
  real :: time(1:10)
  character(80) :: minimizerName

  call energy_and_gradient( f, df )
  
end subroutine energy_minimization
