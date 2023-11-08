	program brownian
	implicit double precision(a-h,o-z)

	parameter ( men=300 )
	common/pos/x0(men),y0(men),z0(men)
	common/vel/x1(men),y1(men),z1(men)
	common/for/fx(men),fy(men),fz(men)
	common/prop/sumv,xsum,isum
        common/sig/sigma(men,men),sigma0,cut1,cut2,eshift
	common/bon/bond,H1,H2,b(men-1)
        common/nat/xn(men),yn(men),zn(men)
	COMMON/SAMPLE/ETA(men,men),CR(men,men)
        COMMON/PARA/CK,AR0,ALPHA
        COMMON/MONO/N
	common/cle/kb
        common/cutoff/rcut
	character*100 fname1,fname2,fname3,residue*2,residuename*5

        dimension aaa(100),dis_b(men)


	common/force/fstr

        residue=" P"
        residuename="    C"
        nchain=1

	fstr=0.000
c	tr=0.125d0
	part=men
        N=men

        write(fname1,"(a,f5.3,a)")
     &'thermo_',fstr,'_300bp_dt0.01_topo_endfree.dat'
	open(14,file=fname1,status='unknown')


        write(fname3,'(a,f5.3,a)')
     &'dishb_',fstr,'_300bp_dt0.01_topo_endfree.dat'
	open(15,file=fname3,status='unknown')

        rcut=1.122d0

	! SET PARAMETERS

        CK=100.0d0 !! ELASTIC CONSTANT
        AR0=1.122d0  !! MINIMUM POSITION OF THE L-J POTENTIAL
        ALPHA=0.0d0

	iflg  = 0	! save to 'simdata.kat' when iflg=1

c        CALL GOSAMPLE
	open(20,file='garcoords_tight.txt',status='old')

        DO I=1,men
	read(20,*)xn(i),yn(i),zn(i)
        x0(i)=xn(i)
        y0(i)=yn(i)
	z0(i)=zn(i)
        ENDDO
	close(20)

	CALL FORENER(AVCF,ENER)

	! READ PARAMETERS FROM FILE
	open(1,file='test1.inp',status='old')
	read(1,*)H1
	read(1,*)H2
c	read(1,*)tr
	read(1,*)delta
	read(1,*)GAMMA
	read(1,*)maxkb
	read(1,*)nconf
	read(1,*)kflg
	read(1,*)istart
	read(1,*)kwrite
	read(1,*)dfold
	read(1,*)ieq
	read(1,*)icq
	close(1)

	open(2, file='seed.inp', status='old')
	read(2,*)iseed
	close(2)


	! SET PARAMETERS IN PREDICTOR-CORRECTOR METHOD

	DeltEuler=delta/GAMMA

	! SCALE FACTORS FOR VELOCITIES DURING EQUILIBRATION
	pi=3.14159265d0
	twopi=8.d0*datan(1.d0)
	delsq=delta*delta
	deltsq=0.5d0*delsq

	! LANGEVIN PARAMETERS
	xmass=1.d0


	aaa(1)=0.100d0
	aaa(2)=0.205d0
		aaa(3)=0.210d0
    	aaa(4)=0.215d0
        aaa(5)=0.220d0
        aaa(6)=0.225d0
        aaa(7)=0.230d0
        aaa(8)=0.235d0
        aaa(9)=0.240d0
        aaa(10)=0.245d0
		aaa(11)=0.250d0
		aaa(12)=0.255d0
		aaa(13)=0.260d0
	aaa(14)=0.265d0
	aaa(15)=0.270d0
	aaa(16)=0.275d0
	aaa(17)=0.280d0
	aaa(18)=0.285d0
	aaa(19)=0.290d0
	aaa(20)=0.295d0



c        iseed=674673

	dfr=rand(iseed)
	loop=1.0
c	open(3, file='status.txt', status='unknown')
	do ii=1,1
	tr=aaa(ii)
c        write(3,*)tr
        write(fname2,'(a,f5.3,a,f5.3,a)')
     &'traj_',fstr,'_',tr,'_300bp_dt0.005_topo_endfree.dat'

	open(13,file=fname2,status='unknown')

        write(13,*)nchain
        write(13,*)men


        do i=1,men
        write(13,"(I4,a2,a5)")i,residue !,residuename
        enddo

        write(13,"(I10)")men


        aheat=delsq*part*3*tr

	const2=2.d0*tr*gamma/(delta*xmass)
	const2=dsqrt(const2)

	! LOAD INITIAL VELOCITIES OF PARTICLES
	call Intvel3d(aheat,part)

        CALL brownianforce(twopi,const2)

	time=0.d0
	tmp=tr
	ekin=part*tr
	etot=epot+ekin
	call tempe(part,delta,tmpins)

c -----------------------------------------
c	ENTER MAIN LOOP OF SIMULATION
c -----------------------------------------
	kb=0
	sume=0.d0
	sume2=0.d0

        do ic=1,men/2
        dis_b(ic)=0.0d0
        enddo



533	continue
	kb=kb+1

        CALL FORENER(AVCF,epot)


 	call brownianforce(twopi,const2)


	 do 37 i=1,men
	 if(i.eq.1.or.i.eq.men) goto 37
           x0(i)=x0(i)+DeltEuler*fx(i)
           y0(i)=y0(i)+DeltEuler*fy(i)
           z0(i)=z0(i)+DeltEuler*fz(i)
37       continue

	if(kb.gt.ieq)then
	sume=sume+epot
	sume2=sume2+epot*epot
	endif

	call tempe(part,delta,tmpins)	! Instant. Temp.

c        Call CPU_time(t2)

c        if((kb.gt.ieq).and.(kb/1000)*1000.eq.kb)then
        if((kb/1000)*1000.eq.kb)then
        write(13,*)kb
        do ik=1,men
        write(13,"(3f13.3)")x0(ik),y0(ik),z0(ik)
        enddo
        if(kb.lt.maxkb)then
        write(13,*)'continue'
        endif
        endif


        if(kb.gt.ieq)then
        nb=0
        do ib1=1,men/2
        do ib2=((men/2)+1), men
c        do ib2=men/2+1, MEN
        if((ib1+ib2).eq.(men+1))then
        nb=nb+1
        dx=x0(ib1)-x0(ib2)
        dy=y0(ib1)-y0(ib2)
        dz=z0(ib1)-z0(ib2)
        r=sqrt(dx**2+dy**2+dz**2)
        dis_b(nb)=dis_b(nb)+r
        endif
        enddo
        enddo
        endif

	if(kb.le.maxkb)goto 533

c        Call System("bzip2 " // fname2)

	avge=sume/(maxkb-ieq)
	avge2=sume2/(maxkb-ieq)
	avgn=cont/(maxkb-ieq)
	fluc=avge2-avge**2

        write(14,"(4f15.5)")fstr,tr,avge,(fluc/(tr**2))
c        write(3,"(4f15.5)")fstr,tr,avge,(fluc/(tr**2))
c        Call CPU_time(t2)
c        write(*,*)kb,t2


          do ik1=1,men/2
          write(15,"(f6.3,I4,f10.3)")tr,ik1,dis_b(ik1)/(maxkb-ieq)
c          write(*,"(f6.3,I4,f10.3)")tr,ik1,dis_b(ik1)/(maxkb-ieq)
          enddo

	close(13)
	enddo

          close(14)
          close(15)
c	   close (3)

        ! END LOOP OF SIMULATION


	! END OF EXECUTION
 	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine intvel3d(aheat,part)
	implicit double precision(a-h,o-z)

	! ASSIGN INITIAL VELOCITIES TO ATOMS

	parameter ( men=300 )
	common/vel/x1(men),y1(men),z1(men)

	sumx=0.
	sumy=0.
	sumz=0.
	do 200 i=1,men
c	if(i.eq.1.or.i.eq.(MEN/2).or.i.eq.((MEN/2)+1).or.i.eq.MEN)goto 200
        xx=2.d0*(rand(0)-0.5d0)
        yy=2.d0*(rand(0)-0.5d0)
        zz=2.d0*(rand(0)-0.5d0)
        xyz=1.d0/ dsqrt(xx*xx+yy*yy+zz*zz)
        x1(i)=xx*xyz
        y1(i)=yy*xyz
        z1(i)=zz*xyz
        sumx=sumx+x1(i)
        sumy=sumy+y1(i)
        sumz=sumz+z1(i)
200	continue

	! SCALE VELOCITIES SO THAT TOTAL MOMENTUM = ZERO
	x=0.d0
	do 210 i=1,men
c	if(i.eq.1.or.i.eq.(men/2).or.i.eq.((men/2)+1).or.i.eq.men)goto 210
	  x1(i)=x1(i)-sumx/part
	  y1(i)=y1(i)-sumy/part
	  z1(i)=z1(i)-sumz/part
	  x=x+x1(i)*x1(i)+y1(i)*y1(i)+z1(i)*z1(i)
210	continue

	! SCALE VELOCITIES TO DESIRED TEMPERATURE
	heat= dsqrt(aheat/x)
	do 220 i=1,men
c	if(i.eq.1.or.i.eq.(men/2).or.i.eq.((men/2)+1).or.i.eq.men)goto 220
	  x1(i)=x1(i)*heat
	  y1(i)=y1(i)*heat
	  z1(i)=z1(i)*heat
220	continue
	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine brownianforce(twopi,const2)
        implicit double precision(a-h,o-z)
        parameter(men= 300)
        common/vel/x1(men),y1(men),z1(men)
        common/for/fx(men),fy(men),fz(men)
	common/cle/kb
        common/cutoff/rcut

        do 26 i=1,men
c	if(i.eq.1.or.i.eq.(MEN/2).or.i.eq.((MEN/2)+1).or.i.eq.MEN)goto 26
101        r1=rand(0)
           r2=rand(0)
           r3=rand(0)
           r4=rand(0)
           r5=rand(0)
           r6=rand(0)

           if(r1.eq.0.or.r2.eq.0.or.r3.eq.0.or.r4.eq.0.or.r5.eq.0.or.
     &r6.eq.0)goto 101
           gam1=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
           gam2=dsqrt(-2.d0*log(r3))*dcos(twopi*r4)
           gam3=dsqrt(-2.d0*log(r5))*dcos(twopi*r6)

           fx(i)=fx(i)+const2*gam1  !gamma*delta
           fy(i)=fy(i)+const2*gam2
           fz(i)=fz(i)+const2*gam3

26      continue

         return
         end


c-------------------------------------------------------------------
	subroutine tempe(part,delta,tmpins)
        implicit double precision(a-h,o-z)
        parameter(men=300)
	common/vel/x1(men),y1(men),z1(men)

	! Center of Mass
        vxcm=0.d0
        vycm=0.d0
	vzcm=0.d0
        do i=1,men
        vxcm=vxcm+x1(i)
        vycm=vycm+y1(i)
        vzcm=vzcm+z1(i)
        enddo
        anorm=1.d0/part
        vxcm=vxcm*anorm
        vycm=vycm*anorm
        vzcm=vzcm*anorm

	! Directional instantaneous temperatures
        sumvelx=0.d0
        sumvely=0.d0
        sumvelz=0.d0
        do i=1,men
	x=x1(i)-vxcm	!  center of mass subtracted
	y=y1(i)-vycm
	z=z1(i)-vzcm
        sumvelx=sumvelx + x*x
        sumvely=sumvely + y*y
        sumvelz=sumvelz + z*z
        enddo
        anorm = 1.d0/(part*delta*delta)
        tmpx=sumvelx*anorm
        tmpy=sumvely*anorm
        tmpz=sumvelz*anorm
        tmpins = (tmpx+tmpy+tmpz)/3

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        SUBROUTINE FORENER(AVCF,ENER)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)

	PARAMETER (men=300)
        DIMENSION CX(men),CY(men)
        DIMENSION CFX(20),CFY(20)
        common/pos/x0(men),y0(men),z0(men)
        common/for/fx(men),fy(men),fz(men)
        COMMON/SAMPLE/ETA(men,men),CR(men,men)
        COMMON/PARA/CK,AR0,ALPHA
        COMMON/MONO/N
        common/force/fstr
	common/cle/kb
        common/cutoff/rcut

	ENER=0.d0
        do i=1,men
        fx(i) = 0.d0
        fy(i) = 0.d0
        fz(i) = 0.d0
        enddo


        do 495 i=1,men-1
        xi=x0(i)
        yi=y0(i)
        zi=z0(i)
        do 496 j=i+1,men
           dx = xi-x0(j)
           dy = yi-y0(j)
           dz = zi-z0(j)
           rsq=dx*dx+dy*dy+dz*dz
           r = dsqrt(rsq)


          IF(J.EQ.(I+1))THEN !! HARMONIC PART
          FCE=-2*CK*(R-AR0)/R
          ENE=CK*(R-AR0)**2
          ELSE                !! L-J PART


c	  if((i.eq.1.and.j.eq.(men/2)).or.
c    &(i.eq.((men/2)+1).and.j.eq.men))then
c          FCE=-2*CK*(R-1.12d0)/R
c          ENE=CK*(R-1.12d0)**2
c	  goto 775
c	  endif


          if(r.gt.rcut)goto 496
c	  if(i.eq.18.and.j.eq.19)then
c	  write(*,*)i,j, CR(I,J), ETA(I,J)
c	  pause
c	  endif

          AA=1/RSQ**6
          BB=1/RSQ**3
          ENE=4.0d0*(AA-BB+0.25)
          FCE=48.0d0*AA-24.0d0*BB
          FCE=FCE/RSQ
          ENDIF

775        continue

c          if(kb.ge.5503.and.kb.le.5505)then
c          write(*,*)FCE
c          endif


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           repx=fce*dx
           repy=fce*dy
           repz=fce*dz

           fx(i) = fx(i) + repx
           fx(j) = fx(j) - repx
           fy(i) = fy(i) + repy
           fy(j) = fy(j) - repy
           fz(i) = fz(i) + repz
           fz(j) = fz(j) - repz
           ENER=ENER+ENE
496     continue
495     continue

c	wall force
	do 613 i=5,men-4
	  dz1=z0(1)-z0(i)
	  dz2=z0(i)-z0(men)
c	  write(*,*) i, x0(i), x0(1), x0(men)
	    if (dz1.lt.rcut) then
	      AA=1/dz1**12
              BB=1/dz1**6
              ENE=4.0d0*(AA-BB)
              ENER=ENER+ENE
              FCE=48.0d0*AA-24.0d0*BB
              fz(i)=fz(i)-FCE
c 	      write(*,*) FCE, ENE, i, dy1             
            else if (dz2.lt.rcut) then
              AA=1/dz2**12
              BB=1/dz2**6
              ENE=4.0d0*(AA-BB)
              ENER=ENER+ENE
              FCE=48.0d0*AA-24.0d0*BB
              fz(i)=fz(i)+FCE
c             write(*,*) FCE, ENE, i, dy2
            endif
613 	continue



C        external force
          if(y0(1).gt.y0(men))then
          fy(1)=fy(1)+fstr
          fy(men)=fy(men)-fstr
          else
          fy(1)=fy(1)-fstr
          fy(men)=fy(men)+fstr
          endif


C     CALCULATE AVERAGE FORCE F
            AVCF=0.D0
            DO 782 I=1,N
            AVCF=AVCF+(FX(I)**2+FY(I)**2+FZ(I)**2)
782         CONTINUE
            AVCF=DSQRT(AVCF)/N
         RETURN
         END


ccc-------------------------------------------------------------
        SUBROUTINE GOSAMPLE
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER(men=300)
        COMMON/SAMPLE/ETA(men, men),CR(men, men)
        COMMON/PARA/CK,AR0,ALPHA
        COMMON/MONO/N

        DO I=1,N-1
	DO J=I+1,N

        IF(J.EQ.(I+1).and.j.ne.((men/2)+1))THEN
        ETA(I,J)=0.D0
        CR(I,J)=0.D0
        ELSE


        if((i+j).eq.(men+1))then
        ETA(I,J)=1.D0 !!! ALPHA=0.
        CR(I,J)=1.D0
        else
        ETA(I,J)=ALPHA !!! ALPHA=0.
        CR(I,J)=1.D0
        endif

c	pause

        ENDIF
        ENDDO
        ENDDO


C   PUT NONZERO NATIVE CONTACT ENERGIES FOR TARGET (CELL 63)




	 DO I=1,N-1
         DO J=I+1,N
         ETA(J,I)=ETA(I,J)
         CR(J,I)=CR(I,J)    !!!! SYMMETRY
         ENDDO
         ENDDO

        RETURN
        END
