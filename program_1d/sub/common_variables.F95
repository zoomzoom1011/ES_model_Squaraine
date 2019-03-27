!*********************************************************************!
!    Module with all variables
!*********************************************************************!
module common_variables
    implicit none

    !constant
    real*8::pi = atan(1.d0)*4 !define later
        
    !system parameters
    integer:: es_zmax
    integer:: nmax = 2
    integer:: es_n = 1, es_z1 = 2, es_z2 = 3

    !geometry
    real*8 x_offset
    real*8 y_offset
    real*8 z_offset
    real*8 armlength
    real*8, allocatable:: coulomb_coupling(:)
    real*8, allocatable:: molecule_pos(:, :, :)
    real*8, allocatable:: tq(:, :)
    real*8:: es_lambda_z1(3)=(/0.d0,1.d0,0.d0/) !N, Z1, Z2 for arm z1
    real*8:: es_lambda_z2(3)=(/0.d0,0.d0,1.d0/) !N, Z1, Z2 for arm z2, shift if in the specified electronic state
    
    !vibration
    integer:: vibmax=1
    integer:: sys_vibmax = 2
    real*8:: hw = 1400.d0
    real*8:: s_z = 1.00d0
    real*8:: s_cat = 0.50d0
    real*8:: s_ani = 0.50d0
    
    !CT coupling
    real*8:: tz
    real*8:: nz
    real*8:: nct
    real*8:: tct
    
    !transition dipole moment
    real*8,allocatable:: ux(:)
    real*8,allocatable:: uy(:)
    real*8,allocatable:: uz(:)
    
    !task title
    character*256:: task_title = 'default'
    
    !coupling
!    real*8:: jo_x = 350
    real*8:: dielectric = 1.d0
!    logical :: nofrenkel = .false.
!    real*8 :: tdm_scaling_factor = 1.d0
    
    !monomer states
    integer:: mon_kount
!    real*8, allocatable:: mon_h(:, :)
!    real*8, allocatable:: mon_eval(:)
    
    !system states
    integer:: sys_kount
    real*8, allocatable:: sys_h(:, :)
    real*8, allocatable:: sys_eval(:)
    
    !partition function 
!    real*8,allocatable:: partition(:, :)
    
    !franck condon factors
!    real*8,allocatable:: fc_gz(:,:)
!    real*8,allocatable:: fc_gc(:,:)
!    real*8,allocatable:: fc_ga(:,:)
!    real*8,allocatable:: fc_cz(:,:)
!    real*8,allocatable:: fc_az(:,:)
    
    
    !multiparticle states
    logical:: ct_state = .true.
    
    !oscillator strength
    real*8,allocatable::ab_osc_x(:)
    real*8,allocatable::ab_osc_y(:)
    real*8,allocatable::ab_osc_z(:)    
    real*8,allocatable::ab_x(:)
    real*8,allocatable::ab_y(:)
    real*8,allocatable::ab_z(:)
    logical::calc_pl = .false.
    real*8,allocatable::pl_osc_x(:,:)
    real*8,allocatable::pl_osc_y(:,:)
    real*8,allocatable::pl_osc_z(:,:)
    real*8,allocatable::pl_x(:)
    real*8,allocatable::pl_y(:)
    real*8,allocatable::pl_z(:)
    real*8::abs_lw = 200.d0
    logical:: no_frenkel = .false.
    integer:: spec_step
    real*8 spec_start_ab
    real*8 spec_end_ab
    real*8 spec_start_pl
    real*8 spec_end_pl
                
    !lineshape
    logical::lorentzian = .false.
    
    !harmonic potential well shift
    real*8:: lambda_n
    real*8:: lambda_z1
    real*8:: lambda_z2
    real*8:: lambda_a
    real*8:: lambda_c
    
    !normalization for spectra
!    logical:: normalization = .t.

    !periodic 
    logical:: periodic = .false. 
    
    type basis
        integer::es_state
        integer::vib1
        integer::vib2
    end type basis
    
    type(basis), allocatable:: mon_state(:)
    type(basis), allocatable:: sys_state(:,:)
    
    !z_state parts
    integer:: leftdonor = 1
    integer:: acceptor = 2
    integer:: rightdonor = 3
    
end module