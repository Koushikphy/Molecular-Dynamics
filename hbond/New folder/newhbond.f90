subroutine hbond(tmp,side,h_pair) !returns # of h-bond each molecule has in a particular time
    implicit none
    real (kind=8), intent(in) :: tmp(768,3),side
    integer, intent(out)      :: h_pair(32640)
    real (kind=8)             :: o1(3),o2(3),h(3),h1(3),dis1,dis2,theta
    integer                   :: i,j,c,pair(256,256)
    pair(:,:)=0
    iloop : do i =1,256                             !iterate over donor
        o1=tmp(3*i-2,:)
        h=tmp(3*i-1,:)
        h1=tmp(3*i,:)
        jloop : do j= 1,256                         !iterate over acceptor
            if (i==j) then                          !Skip same molecule
                cycle
            end if
            o2=tmp(3*j-2,:)
            o2=o2-side*anint((o2-o1)/side)          !Wrap o2 applying PBC
            dis1=dsqrt(sum((o1-o2)**2))
            if (dis1<=3.5) then
                dis2=dsqrt(sum((h-o2)**2))          !check 1st hydrogen of donor
                if (dis2<=2.45) then
                    call angle(o2,o1,h,theta)
                    if (theta<=30) then
                        pair(i,j)=1
                        pair(j,i)=1
                    endif
                endif
                dis2=dsqrt(sum((h1-o2)**2))         !check 2nd hydrogen of donor
                if (dis2<=2.45) then
                    call angle(o2,o1,h1,theta)
                    if (theta<=30) then
                        pair(i,j)=1
                        pair(j,i)=1
                    endif
                endif
            endif
        end do jloop
    end do iloop

    c=1
    kloop : do i=1,256
        tloop : do j=i+1,256
            h_pair(c)=pair(i,j)
            c=c+1
        end do tloop
    end do kloop

end subroutine hbond

subroutine angle(x,y,z,theta)                       !returns angle between points x,y,z in degrees
    real (kind=8), intent(in) :: x(3),y(3),z(3)
    real (kind=8), intent(out):: theta
    real (kind=8)             :: a(3),b(3),v
    a=x-y
    b=z-y
    v=dot_product(a,b)
    theta=acos(v/(norm2(a)*norm2(b)))*(180.0/3.14159265)
end subroutine angle



subroutine pop(hframe,shb,chb,n)                       !returns angle between points x,y,z in degrees
    real (kind=8), intent(in) :: hframe(n,32640)
    integer ,      intent(in) :: n
    real (kind=8), intent(out):: shb(0:n-1),chb(0:n-1)
    real (kind=8)             :: avg_h
    integer                   :: i,c,d,t,tau
    avg_h=0
    kloop : do i=1,n 
        avg_h=avg_h+sum(hframe(i,:))
    end do kloop
    avg_h=avg_h/(32640.0*n)
    tloop : do t=0,n-1
        c=0
        d=0
        ploop : do i=1,32640
            tauloop : do tau=1,n-t
                if (hframe(tau,i)==1) then 
                    if (hframe(tau+t,i)==1) then 
                        d=d+1 
                        if (all(hframe(tau+1:tau+t-1,i)==1)) then 
                            c=c+1
                        endif 
                    endif
                endif
            end do tauloop
        end do ploop
        shb(t)=log(c/(avg_h*32640.0*(n-t)))
        chb(t)=log(d/(avg_h*32640.0*(n-t)))
        print *,t
    end do tloop
end subroutine pop

