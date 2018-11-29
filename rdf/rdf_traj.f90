subroutine homo(atom,side,dr,bins,hist,n)           !rdf for same type element 
    implicit none
    integer, intent(in)       :: n,bins
    real (kind=8), intent(in) :: atom(n,3), side,dr
    real (kind=8), intent(out):: hist(bins)
    real (kind=8 )            :: x(3),rd
    integer                   :: i,j,bi
    hist(:)=0
    iloop : do i =1,n
        jloop : do j= i+1,n                         !skip counting same atom pair
            x = atom(i,:)-atom(j,:)                 !x is a 3D vector between position of two atom
            x=x-side*anint(x/side)                  !Applying PBC
            rd = norm2(x)                           !norm of the vector/distance of two points
            bi=ceiling(rd/dr)           
            if (bi<=bins) then
                hist(bi)= hist(bi)+2
            endif
        end do jloop
    end do iloop
end subroutine homo


subroutine hetero(atom1,atom2,side,dr,bins,hist,m,n)!rdf for different type 
    implicit none
    integer, intent(in)       :: m,n,bins
    real (kind=8), intent(in) :: atom1(m,3), atom2(n,3), side,dr
    real (kind=8), intent(out):: hist(bins)
    real (kind=8 )            :: x(3),rd
    integer                   :: i,j,bi
    hist(:)=0
    iloop : do i =1,m
        jloop : do j= 1,n
            x = atom1(i,:)-atom2(j,:)
            x=x-side*anint(x/side)                  !apply PBC
            rd = norm2(x)
            bi=ceiling(rd/dr)           
            if (bi<=bins) then
                hist(bi)= hist(bi)+1
            endif
        end do jloop
    end do iloop
end subroutine hetero


subroutine normalise(hist,side,dr,m,n,r,rdf,rcn,kb,frame,bins)
    implicit none
    integer, intent(in)       :: frame,bins,m,n
    real (kind=8), intent(in) :: hist(frame,bins), side,dr
    real (kind=8), intent(out):: r(bins), rdf(bins), rcn(bins),kb(bins)
    real (kind=8 )            :: vs,con1,con2,tmp,tmp2
    integer                   :: i,j

! volume of the shell
! =(4/3)*pi*(r_upper^3 - r_lower^3)
! =(4/3)*pi*((i)^3 -(i-1)^3)*dr^3
! =(4/3)*pi*(1 + 3*i*(i-1))*dr^3
! =(4/3)*pi*dr^3 + 4*pi*i*(i-1)*dr^3
    con1=4*acos(-1.0)*dr**3
    con2=con1/3.0
    tmp=0
    tmp2=0
    iloop: do i = 1,bins
        rdf(i)=sum(hist(:,i))
        r(i)=(i-.5)*dr
        vs=con2+con1*i*(i-1)
        tmp=tmp+rdf(i)/(frame*m)
        rcn(i)=tmp                                  !get the RCN data
        rdf(i)=rdf(i)*side**3/(vs*frame*m*n)        !normalise the data  
        tmp2=tmp2+(rdf(i)-1)*vs   
        kb(i)=tmp2                                  !get the KB integral data
    end do iloop
end subroutine normalise
