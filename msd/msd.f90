subroutine msd(dat,m,n,res)                         !return [[t0,msd(t0)],[t1,msd(t1)]..tm,msd(tm)]
    implicit none
    real (kind=8), intent(in) :: dat(n,3)
    integer, intent(in)       :: m,n
    real (kind=8), intent(out):: res(0:m)
    integer                   :: i,j,t
    res(:)=0.0
    tloop : do t=1,m                                !iterate through different lag time
        iloop : do i=1,(n-t)                        !iterate through frames with lag of
            res(t)=res(t)+ sum( (dat(i+t,:)-dat(i,:))**2) !sum of square displacements
        end do iloop
    end do tloop
    res=res/(n-t)                                   !normalise data with frame numbers
end subroutine msd
