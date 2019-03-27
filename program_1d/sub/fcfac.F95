!*********************************************************************!
! franc-condon calculations     
! this subroutine is from haj's program 'haj_polamer.f95'
! i have made no changes
!*********************************************************************!
subroutine fcfac(n,m,s,fc)
    implicit none
    integer n,m,k
    real*8 s,fc,f_m,f_n,f_k,f_nmk,f_mk,facin
    real*8 fact

    fc = 0.d0

    do k = 0,m
        ! if n-m+k is negative, factorial is not calculatable.
        if(n-m+k < 0) go to 100    

        f_mk  = fact(m-k)
        f_nmk = fact(n-m+k)
        f_k   = fact(k)
        facin = 1.d0/(1.d0*f_k*f_mk*f_nmk)

        fc = fc + facin*s**(k*0.5d0)*s**(1.0d0*(n-m+k)*0.5d0)*&
                  (-1)**(n-m+k)
100        continue
    enddo

    f_n = fact(n)
    f_m = fact(m)
    fc = fc*dsqrt(1.d0*f_m*f_n)*dexp(-s/2.d0)

    return
end subroutine