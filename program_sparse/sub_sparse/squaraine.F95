program squaraine
    use common_variables
    implicit none
    
#include<slepc/finclude/slepcdef.h>
#include<petsc/finclude/petscdef.h>
    
    real*8:: total_time, time_in_seconds
    PetscViewer :: HOUT

    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
    total_time = MPI_Wtime()
    call MPI_COMM_RANK( PETSC_COMM_WORLD, rank, ierr )
    call MPI_COMM_SIZE( PETSC_COMM_WORLD, rankmax, ierr )
    
    !read para file
    call read_para_file()

    !read the crystal structure
    ! call read_in_crystal()

    !read geometry from xyz file
    ! call read_geom()
    
    if(rank == 0) then
        call system('if [ ! -d '//trim(task_title)//' ]; then '//&
                    'mkdir '//trim(task_title)//'; fi' )
        if( rank == 0 ) print*, 'Entering ./'//trim(adjustl(task_title))
        call chdir( task_title )
    end if
    
    !index the basis states
    call build_molecule_geometry()
    call index_mon_state()
    call index_sys_state()

    !calc coupling from geometry
    ! call calc_cpl_from_geom()
    
    !output the kounter
    call kount_out()

    !set the coupling and transition dipole moment
    call calc_coupling()

    !building the hamiltonian
    !output the hamiltonian to sparse matrix
    call hamiltonian_out()

    ! call MatView( h_sparse, PETSC_VIEWER_STDOUT_WORLD, ierr )
    !open viewer to write hamiltonian to file
    !PRINT HAMILTONIAN 
    ! call PetscViewerASCIIOpen( PETSC_COMM_WORLD, trim(task_title)//'_H.dat', HOUT, ierr )
    ! call MatView( h_sparse, HOUT, ierr )
    
    !adjust request evec to make sure it is not asking for more than possible
    nev = min( nev, sys_kount )
    if (rank==0) print'(a,i4,a)', ' Now looking for ', nev,' eigenvectors.'
    
    !slepc diagonalize
    !Create eigensolver context
    call EPSCreate( PETSC_COMM_WORLD, eps_cxt_sys, ierr )
    !Set operators. It is a standard eigenvalue problem
    call EPSSetOperators( eps_cxt_sys, h_sparse, PETSC_NULL_OBJECT, ierr )
    call EPSSetProblemType( eps_cxt_sys, EPS_HEP, ierr )
    !Set solver parameters at runtime
    call EPSSetFromOptions( eps_cxt_sys, ierr )
    !seek smallest real parts
    call EPSSetWhichEigenpairs( eps_cxt_sys, EPS_SMALLEST_REAL, ierr )
    call EPSSetDimensions( eps_cxt_sys, nev, PETSC_DEFAULT_INTEGER,  &
                           PETSC_DEFAULT_INTEGER, ierr )
    time_in_seconds = MPI_Wtime() 
    call EPSSolve(eps_cxt_sys, ierr)
    time_in_seconds = MPI_Wtime() - time_in_seconds
    if (rank==0) print'(a,f10.4)', ' Diagonalization Time: ', time_in_seconds

    !output eigenvalues and eigenvectors
        
    !calculate the absorption spectrum
    call abs_osc()
    call abs_spec()
    call abs_out()
    call normalization_ab()

    !calculate the pl spectrum
    if (calc_pl) call pl_osc()
    if (calc_pl) call pl_spec()
    if (calc_pl) call pl_out()
    if (calc_pl) call normalization_pl()

    !write the para to files
    call para_out()
    
    call EPSDestroy( eps_cxt_sys, ierr )
    call MatDestroy( h_sparse, ierr )

    total_time = MPI_Wtime() - total_time
    if (rank==0) print'(a,f10.4)', ' Calculation finished. Total Time: ', total_time
    call SlepcFinalize( ierr )

end program
