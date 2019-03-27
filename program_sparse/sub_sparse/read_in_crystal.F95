!*********************************************************************!
!    read in the crystal structure
!*********************************************************************!
subroutine read_in_crystal()
    use common_variables
    implicit none

#include<petsc/finclude/petscdef.h>
    integer :: fno = 105
    character*100 :: fname, readline1
    real*8 cnv_to_real
    character*128:: cwd
    Integer:: istat, getcwd
    Integer:: nline, n
    
    
    fname = crystal_fname
    if ( rank==0 ) print*, 'Reading Crystal Data from '//trim(fname)
    if ( rank==0 ) print*, '**********************************'//&
            '**********************************'
    !Get the current working directory
    istat = getcwd( cwd )
    if ( istat .ne. 0 ) stop 'getcwd: error'

    !Adjust the geometry file name
    fname = trim(cwd)//trim(fname)

    nline = 0
    OPEN (fno, file = fname) 
    DO 
        READ (fno,*, END=101) 
        nline = nline + 1 
    END DO 
    101 CLOSE (fno)

    open( unit = fno, file = trim(fname),          &
          status = 'old', action = 'read' )
    do
        read( fno, '(a)', end = 103) readline1
        if ( (readline1(1:14)) == '_cell_length_a' ) exit
    end do
    
    cell_a = cnv_to_real(readline1)
    read( fno,'(a)') readline1
    cell_b = cnv_to_real(readline1)
    read( fno,'(a)') readline1
    cell_c = cnv_to_real(readline1)
    read( fno,'(a)') readline1
    cell_alpha = cnv_to_real(readline1)/180.d0*pi
    read( fno,'(a)') readline1
    cell_beta = cnv_to_real(readline1)/180.d0*pi
    read( fno,'(a)') readline1
    cell_gamma = cnv_to_real(readline1)/180.d0*pi
    if ( rank==0 ) print*, 'crystal:',cell_a,cell_b,cell_c
    if ( rank==0 ) print*, 'crystal:',cell_alpha,cell_beta,cell_gamma

    103 close(fno)
end subroutine
