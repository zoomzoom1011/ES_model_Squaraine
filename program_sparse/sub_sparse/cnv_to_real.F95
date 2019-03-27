!*********************************************************************!
!    convert character to real
!*********************************************************************!
real*8 function cnv_to_real( numin )
    implicit none
    integer i,j
    character*100 numin, numout
    integer:: pos
    integer :: stringLen 
    integer :: last, actual
    
    !the cif file has numbers in the form ####(#)
    !function removes the parenthesis and returns the number
    
    pos = scan(numin, ' ')
    numin = numin(pos + 1: )

    stringLen = len (numin)
    last = 1
    actual = 1

    do while (actual < stringLen)
        if (numin(last:last) == ' ') then
            actual = actual + 1
            numin(last:last) = numin(actual:actual)
            numin(actual:actual) = ' '
        else
            last = last + 1
            if (actual < last) &
                actual = last
        endif
    end do

    j = 1
    do i = 1, stringLen
        if ( numin(i:i) /= '(' .and. numin(i:i) /= ')' ) then
            numout(j:j) = numin(i:i)
            j = j + 1
        end if
    end do
    
    read( numout, * ) cnv_to_real    
end function