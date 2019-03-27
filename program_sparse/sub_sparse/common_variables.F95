!*********************************************************************!
!    Module with all variables
!*********************************************************************!
module common_variables
    implicit none
    
!include petsc and slepc headers
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include <slepc/finclude/slepcsys.h>
#include <slepc/finclude/slepceps.h>


    !system parameters
    character*256:: task_title = 'default'
    integer:: nmax = 2
    integer:: es_n = 1, es_z1 = 2, es_z2 = 3
    
    !geometry
    ! real*8:: cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma
    ! character*200:: crystal_fname
    ! character*200:: xyz_file(4)
    real*8:: x_offset
    real*8:: y_offset
    real*8:: z_offset
    real*8:: armlength
    real*8:: dielectric
    integer:: sublattice = 2
    integer:: kount_duplicate_stacks
    real*8, allocatable:: coulomb_coupling(:)
    real*8, allocatable:: mol1pos(:, :, :)
    ! real*8, allocatable:: tq(:, :)
    real*8:: to_cartesian(3,3)
    real*8:: agg_angle
    
    !vibration
    integer:: vibmax = 1
    integer:: sys_vibmax = 2
    real*8:: hw = 1400.d0
    ! real*8:: s_z = 1.00d0
    ! real*8:: s_cat = 0.50d0
    ! real*8:: s_ani = 0.50d0
    
    !CT coupling
    real*8:: tz
    real*8:: nz
    real*8:: nct
    real*8:: tct
    
    !indexes
    integer, allocatable:: nx_lattice(:,:,:)
    
    !transition dipole moment
    real*8,allocatable:: ux(:)
    real*8,allocatable:: uy(:)
    real*8,allocatable:: uz(:)
    
    !coupling
!    real*8:: jo_x = 350
!    real*8:: dielectric = 1.d0
!    logical :: nofrenkel = .false.
!    real*8 :: tdm_scaling_factor = 1.d0
    
    !monomer states
    integer:: mon_kount
!    real*8, allocatable:: mon_h(:, :)
!    real*8, allocatable:: mon_eval(:)
    
    !sparse matrix technique
    Mat         h_sparse
    EPS         eps_cxt_sys
    PetscScalar er, ei        
    PetscScalar ger, gei
    Vec         wr, wi
    Vec         gr, gi
!    Vec         ux, uy, uz
    PetscInt    Istart, Iend
    PetscErrorCode ierr
    PetscInt:: rank
    PetscInt:: rankmax
    PetscInt:: sys_kount
    PetscInt:: nev
    
    !system states
!    integer:: sys_kount
    real*8, allocatable:: ab_sys_eval(:)
    real*8, allocatable:: pl_sys_eval(:,:)
    
    !intermolecule charge transfer
    logical:: ct_state = .false.
    ! integer:: es_a = 4, es_c1 = 5, es_c2=6, es_z3=7
    integer:: leftdonor = 1, acceptor = 2, rightdonor = 3
    real*8  :: es_lambda_z1(3)=[0.d0,1.d0,0.d0] !N, Z1, Z2 for arm z1
    real*8  :: es_lambda_z2(3)=[0.d0,0.d0,1.d0] !N, Z1, Z2 for arm z2, shift if in the specified 
    
    !harmonic potential well shift
    real*8:: lambda_n
    real*8:: lambda_z1
    real*8:: lambda_z2
    ! real*8:: lambda_a
    ! real*8:: lambda_c1
    ! real*8:: lambda_c2
    ! real*8:: lambda_z3
    
    !absorption & pl spectra
    real*8, allocatable:: ab_osc_x(:)
    real*8, allocatable:: ab_osc_y(:)
    real*8, allocatable:: ab_osc_z(:)
    real*8, allocatable:: ab_x(:)
    real*8, allocatable:: ab_y(:)
    real*8, allocatable:: ab_z(:)
    real*8, allocatable:: pl_osc_x(:, :)
    real*8, allocatable:: pl_osc_y(:, :)
    real*8, allocatable:: pl_osc_z(:, :)
    real*8, allocatable:: pl_x(:)
    real*8, allocatable:: pl_y(:)
    real*8, allocatable:: pl_z(:)
    logical:: calc_pl = .false.
    integer:: spec_step = 5000
    real*8:: spec_start_ab,spec_end_ab
    real*8:: spec_start = 0.d0
    real*8:: spec_end = 50000.d0
    real*8:: spec_start_pl,spec_end_pl
    logical:: lorentzian = .false.
    real*8:: abs_lw = 200.d0
    logical:: no_frenkel = .false.
    logical:: nearest_neighbor = .true.
    
    real*8, allocatable:: molecule_pos(:, :, :)
    real*8, allocatable:: tq(:, :)
    ! Real*8, allocatable:: tdm(:,:)
    ! Real*8, allocatable:: abt(:)
    ! Integer:: tmax = 15                    ! max time for TA in fs
    ! Real*8, allocatable:: ab2p(:)
    ! Real*8, allocatable  :: tdmbar(:,:)
    
    !temperature
!    integer:: tnum=13
!    real*8, allocatable :: em_temp(:)
    
    !empty
    integer, parameter :: empty = -1
    
    !normalization for spectra
    logical:: normalization = .true.

    !periodic 
    logical:: periodic = .false. 
    
    !state index
    integer::es_zmax = 2
    
    type basis
        integer::es_state
        integer::vib1
        integer::vib2
    end type basis
    
    type(basis), allocatable:: mon_state(:)
    type(basis), allocatable:: sys_state(:,:)
    
    !physical constants
    !constants for au to cm-1 conversion
    real*8, parameter:: pi = datan(1.d0) * 4
    real*8, parameter:: ev = 8065.54d0
    Real*8, parameter:: hbar = 5.308837367d3          ! cm-1 * fs
end module