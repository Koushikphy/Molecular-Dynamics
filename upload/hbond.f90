! ***********************************************
! File : hbond.f90

! ***********************************************

subroutine hbond(tmp,side,hdist) !returns # of h-bond each molecule has in a particular time
    implicit none
    real (kind=8), intent(in) :: tmp(768,3),side
    integer, intent(out)      :: hdist(256)
    real (kind=8)             :: o1(3),o2(3),h(3),h1(3),dis1,dis2,theta
    integer                   :: i,j
    hdist(:)=0                                      !hist(i) denotes # of H-bond that i-th molecule has
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
                        hdist(i)=hdist(i)+1         
                        hdist(j)=hdist(j)+1         !if donor is making a H-bond then acceptor too
                    endif
                endif
                dis2=dsqrt(sum((h1-o2)**2))         !check 2nd hydrogen of donor
                if (dis2<=2.45) then
                    call angle(o2,o1,h1,theta)
                    if (theta<=30) then
                        hdist(i)=hdist(i)+1
                        hdist(j)=hdist(j)+1
                    endif
                endif
            endif
        end do jloop
    end do iloop
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

