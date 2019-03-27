!*********************************************************************!
!    calculating factorial
!    * i think, mathimaticalla numbers should be integers, but
!      here, i use real numbers for convenience.
!     this subroutine is from haj's program 'haj_polamer.f95'
!     i have made no changes
!*********************************************************************!
real*8 function fact( n )
    implicit none
    integer, intent(in):: n
    integer i

    fact = 1.0d0
    if( n .ge. 0 ) then
        do i = 2, n
            fact = fact * i
        end do
    end if
end function