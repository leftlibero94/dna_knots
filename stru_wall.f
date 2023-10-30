 	implicit double precision (a-h,o-z)
	dimension x(1000),y(1000),z(1000)
c	integer variable names start with i,j,k,l,m,n else dbl

	N=288
	do i=1,N
	x(i)=0
	y(i)=0
	z(i)=0
	enddo
	
	pi=4.d0*datan(1.d0)
        r=N*(2.0**(1/6.0))/pi
	s=180
	ds=5.0/8.0
	radian=pi/180.d0
	write(*,'(3f10.3)') r, ds, radian
	do i=1,N
	theta=radian*s
	x(i)=r*cos(theta)
	y(i)=0
	z(i)=r*sin(theta)
	write(*,'(3f10.3)') s, theta, cos(theta)
	s=s-ds
	enddo
	
	open(3, file='wall_288.txt', status='unknown')
	do i=1,N
	  write(3,4)x(i),y(i),z(i),'  ', i
	enddo
4	format(3f10.3,A,I3)
	close(3)
	
	stop
	end
