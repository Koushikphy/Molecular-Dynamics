subroutine hetero(atom1,atom2,side,dr,bins,r,rdf,m,n)
    implicit none
    integer, intent(in)       :: m,n,bins
    real (kind=8), intent(in) :: atom1(m,3),atom2(n,3), side,dr
    real (kind=8), intent(out):: r(bins), rdf(bins)
    real (kind=8 )            :: xx,yy,zz,rd,vs,con1,con2
    integer                   :: i,j,bi

    do i= 1,bins
        rdf(i)=0
    end do

    iloop : do i =1,m
        jloop : do j= 1,n
            xx = atom1(i,1)-atom2(j,1)
            yy = atom1(i,2)-atom2(j,2)
            zz = atom1(i,3)-atom2(j,3)

            xx=xx-side*anint(xx/side)
            yy=yy-side*anint(yy/side)
            zz=zz-side*anint(zz/side)

            rd = dsqrt(xx*xx + yy*yy + zz*zz)
            bi=ceiling(rd/dr)           
            if (bi<=bins) then
                rdf(bi)= rdf(bi)+1
            endif
        end do jloop
    end do iloop

    con1=4*acos(-1.0)*dr**3
    con2=con1/3.0

    do i = 1,bins
        r(i)=(i-.5)*dr
        vs=con2+con1*i*(i-1)
        rdf(i)=rdf(i)*side**3/(vs*m*n)
    end do


end subroutine hetero


subroutine homo(atom,side,dr,bins,r,rdf,n)
    implicit none
    integer, intent(in)       :: n, bins
    real (kind=8), intent(in) :: atom(n,3), side,dr
    real (kind=8), intent(out):: r(bins), rdf(bins)
    real (kind=8 )            :: xx,yy,zz,rd,vs,con1,con2
    integer                   :: i,j,bi

    do i= 1,bins
        rdf(i)=0
    end do

    iloop : do i =1,n
        jloop : do j= i+1,n
            xx = atom(i,1)-atom(j,1)
            yy = atom(i,2)-atom(j,2)
            zz = atom(i,3)-atom(j,3)

            xx=xx-side*anint(xx/side)
            yy=yy-side*anint(yy/side)
            zz=zz-side*anint(zz/side)

            rd = dsqrt(xx*xx + yy*yy + zz*zz)
            bi=ceiling(rd/dr)           
            if (bi<=bins) then
                rdf(bi)= rdf(bi)+2
            endif
        end do jloop
    end do iloop

    con1=4*acos(-1.0)*dr**3
    con2=con1/3.0

    do i = 1,bins
        r(i)=(i-.5)*dr
        vs=con2+con1*i*(i-1)
        rdf(i)=rdf(i)*side**3/(vs*n**2)
    end do

end subroutine homo