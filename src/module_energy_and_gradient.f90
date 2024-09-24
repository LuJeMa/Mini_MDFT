module module_energy_and_gradient
  use iso_c_binding, only: c_float
  use precision_kinds, only: dp
  implicit none
  private
  type f_type
    real(dp) :: id = 0._dp,&
                ext = 0._dp,&
                exc_cs = 0._dp,&
                exc_cdeltacd = 0._dp,&
                exc_cproj = 0._dp,&
                exc_ck_angular = 0._dp,&
                exc_fmt = 0._dp,&
                exc_wca = 0._dp,&
                exc_3b = 0._dp,&
                exc_b = 0._dp,&
                exc_dipolar = 0._dp,&
                exc_multipolar_without_coupling_to_density = 0._dp,&
                exc_multipolar_with_coupling_to_density = 0._dp,&
                exc_hydro = 0._dp,&
                exc_nn_cs_plus_nbar = 0._dp, &
                tot=0._dp,&
                pscheme_correction=-999._dp,&
                pbc_correction=-999._dp
    integer :: ieval=0
    logical :: apply_energy_corrections_due_to_charged_solute = .false.  ! modified by daniel
  end type
  type (f_type), public :: ff
  public :: energy_and_gradient

contains




subroutine energy_and_gradient (f, df)

    ! In this subroutine one calls the different parts of the total energy
    ! This subroutine is the one called by the minimization stuff
    ! for computing the total energy and associated gradient
    ! FF is the TOTAL ENERGY of the system, it is thus the functional of the density that is minimized by solver
    ! dF_new is the gradient of FF with respect to all coordinates. Remember it is of the kind dF_new ( number of variables over density (ie angles etc))

    use iso_c_binding, only: c_float
    use precision_kinds, only: dp
    use module_solvent, only: solvent
    use module_grid, only: grid
    !use module_energy_ideal_and_external, only: energy_ideal_and_external
    ! use module_energy_cs, only: energy_cs
    ! use module_energy_cdeltacd, only: energy_cdeltacd
    use module_energy_cproj_mrso, only: energy_cproj_mrso
   ! use module_energy_cproj_mrso_mixture, only: energy_cproj_mrso_mixture
    ! use module_energy_cproj_no_symetry, only: energy_cproj_no_symetry
    ! use module_energy_ck_angular, only: energy_ck_angular
    ! use module_energy_luc, only: energy_luc
    ! use module_energy_luc_fast, only: energy_luc_fast
    use module_input, only: getinput
    ! use module_energy_cproj_slow, only: energy_cproj_slow
    !use module_lbfgs_nocedal_mdft, only:lbfgsb
    use proc_ptrs, only : excess_energy
    implicit none

    real(dp), intent(out) :: f
    real(dp), intent(out), optional :: df (grid%no, grid%nx, grid%ny, grid%nz, solvent(1)%nspec)
    real(dp), parameter :: zerodp=0._dp
    real(c_float) :: t(10)
    real(dp) :: fold
    integer :: ns, s
   
    ! This is an horrible trick to not update the count of SCF cycles if we're in a line search.
    ! If we're not in lbfgs, everything is easy: update the count if one needs the gradient only since the gradient is not computed in the line search.
    ! If the minimizer is lbfgs, then we always compute the gradient, even in the line search. Then, what you want it to use the dedicated variable for lbfgs: isave(30)
    !
    
    if (.not. allocated(solvent)) then
      print*, "in energy_and_gradient, solvent()% is not allocated"
      error stop
    end if
    ns = solvent(1)%nspec

    fold=f
    f  = zerodp
    if(present(df)) df = zerodp
    s=1

    
    !
    ! with Luc's routine
    !
    if (solvent(s)%do%exc_cproj) then
        if(present(df)) then
          call cpu_time(t(4))
          call energy_cproj_mrso( ff%exc_cproj, df, print_timers=.true.)
          call cpu_time(t(5))
        !   print*, "ff%exc_cproj_mrso =", real(ff%exc_cproj), " in",t(6)-t(5),"sec"
        else
          call cpu_time(t(4))
          call energy_cproj_mrso( ff%exc_cproj, print_timers=.false.)
          call cpu_time(t(5))
        end if
        f = f + ff%exc_cproj
    end if
    
    ! if (solvent(s)%do%exc_cproj) then
    !     call cpu_time(t(7))
    !     ! call energy_cproj_mrso (ff%exc_cproj, df)
    !     ! call energy_cproj_no_symetry (ff%exc_cproj, df)
    !     ! call energy_luc ( ff%exc_cproj , df )
    !     call energy_luc_fast ( ff%exc_cproj, df)
    !     call cpu_time(t(8))
    !     print*, "ff%exc_cproj      =", ff%exc_cproj,   "in",t(8)-t(7),"sec"
    !     ! stop "energy_and_gradient after call to energy_luc"
    !     ! print*, "ff%exc_cproj - (ff%exc_cs+ff%exc_cdeltacd) =", ff%exc_cproj-(ff%exc_cs + ff%exc_cdeltacd)
    !     f = f + ff%exc_cproj
    !     ! stop
    ! end if

   


    ! if (solvent(s)%do%exc_ck_angular) then
    !     call cpu_time(t(9))
    !     call energy_ck_angular (ff%exc_ck_angular, df)
    !     call cpu_time(t(10))
    !     print*, "ff%exc_ck_angular =", ff%exc_ck_angular,"in",t(10)-t(9),"sec"
    !     f = f + ff%exc_ck_angular
    ! end if

    !
    ! adhoc corrections to the solvation free energy (Hunenberger, pressure etc.)
    ! if you use systems with non-zero net charge (like one ion in water)
    ! and you use lattice schemes (FFT based) methods for computing the external potential
    !
    



    ! if (solvent(s)%do%wca) call lennard_jones_perturbation_to_hard_spheres (ff%exc_wca, df)
    ! if (solvent(s)%do%exc_multipolar_without_coupling_to_density) &
    !         call energy_polarization_multi (ff%exc_multipolar_without_coupling_to_density, df)
    ! if (solvent(s)%do%exc_multipolar_with_coupling_to_density) &
    !         call energy_polarization_multi_with_nccoupling (ff%exc_multipolar_with_coupling_to_density, df)
    ! if (solvent(s)%do%exc_hydro) call energy_hydro (ff%exc_hydro, df)
    ! if (solvent(s)%do%exc_nn_cs_plus_nbar) call energy_nn_cs_plus_nbar (ff%exc_nn_cs_plus_nbar, df)
    ! if (solvent(s)%do%exc_3b) call energy_threebody_faster (ff%exc_3d, df)

    ff%tot = f

    block
        logical, save :: printheader = .true.
        real(dp) :: reldf, Texc, Ttot, Textid, pgtol
        if(printheader) then
          !  write(*,'(A5,12A14)') "#eval","Ftot","Fext","Fid","Fexc","Fb","Cpbc","Cpsch","relF","pgtol","Ttot","Text+id","Texc"
            write(*,'(A5,8A10)') "#eval","Ftot","Fext","Fid","Fexc","Fb","relF","pgtol","Ttot"
            printheader = .false.
        end if
        Texc = t(5)-t(4)
        !Textid = t(2)-t(1)
        Ttot = Texc+Textid
        if (present(df)) then
            pgtol = real(maxval(df))
        else
            pgtol = 0
        end if
        reldf = (fold-f)/maxval([abs(fold),abs(f),1._dp])
        ! write(*,"(I5,12F14.4)") ff%ieval, ff%tot, ff%ext, ff%id, ff%exc_cproj, ff%exc_b, ff%pbc_correction, ff%pscheme_correction, reldf, pgtol, Ttot, Textid, Texc
        write(*,"(I5,5F10.3,2F10.4,F10.2)") ff%ieval, ff%tot, ff%ext, ff%id, ff%exc_cproj, ff%exc_b, reldf, pgtol, Ttot
    end block

end subroutine energy_and_gradient

end module




