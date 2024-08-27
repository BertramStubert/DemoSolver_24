
c     fortran code of a Demo-2D-Solver (euler method)

c     Author: Bertram Stubert
c     Time of Project-Start: 21.02.2024

c     Goal of this code is a demonstration for interested people to get
c     an idea of a 1D-aero-code functionality as simple as possible

c     Application deals with the calculation of a
c     one-block-configuration building a channel with curved wall.


c     !!! There is no warranty due to the code !!!
c     !!! This code is written for demonstration purposes only !!!
c     !!! Do not deliver any results of this code for any projects !!!
      
      include 'param2D.h'
      include 'common2D.h'

      open (unit = 10, file = 'Demo2D_input.dat')
      open (unit = 11, file = 'Demo2D_output.dat')
      open (unit = 12, file = 'Demo2D_konti.dat')
      open (unit = 13, file = 'Demo2D_evmax.dat')
      open (unit = 14, file = 'Demo2D_coord.dat')


c     konwn indices
      im = id
      jm = jd

c     Mainprogram
c     ***********

c----------------------------------------------------------------------
c     phase 1 - prpare all data for iteration in time

c     initialization variables and arrays
c     ***********************************
      call init

c     parameter input with namelist l1
c     ********************************
      call readst

c     read coordinates
c     ****************
      call read_coord

c     grid geometry
c     *************
      call geo

c     initialization of variables
c     ***************************
      call start

c----------------------------------------------------------------------
c     phase 2 - Lax method for update of all conservative variables

c     iteration over nmax timesteps
c     *****************************
      call update

c----------------------------------------------------------------------
c     phase 3 - write out all results

c     postprocessing
c     **************
      call post

c     write result cgns file
c     **********************
      call write_cgns
      call write_boundary


      call cpu_time(time)

      write(11,*) ' '
      write(11,'(a13,f10.2)') 'cpu-seconds: ', time
      write(11,'(a13,f10.2)') 'cpu-minutes: ', time/60.0
      write(11,*) ' '
      

      write(11,*) ' Demo2D ended'
      write(11,*) ' ************'


      close (10)
      close (11)
      close (12)
      close (13)
      close (14)


      stop
      end
      subroutine cfln
c     ===============

      include 'param2D.h'
      include 'common2D.h'

c     calculate local timesteps
c     *************************

c     Scalar product of (velocity + speed of sound) * area at center
c     of volume gives a term like delta(volume)/delta(timestep).
c     The inverse step is the required faktor for multiplying
c     with the fluxbalance

c     calculate speed of sound
      do  j = 1,jm
       do i = 1,im
        aq       = amax1(100.0,cpzucv*p(i,j)/u1(i,j))
        ar1(i,j) = sqrt(aq)
       enddo
      enddo


       do  j = 1,jm-1
        do i = 1,im-1

c        volume mean values for speed of sound and velocity

         rav = 0.25*(ar1(i  ,j  )
     &              +ar1(i  ,j+1)
     &              +ar1(i+1,j  )
     &              +ar1(i+1,j+1))

         vxv = 0.25*(vx (i  ,j  )
     &              +vx (i  ,j+1)
     &              +vx (i+1,j  )
     &              +vx (i+1,j+1))

         vyv = 0.25*(vy (i  ,j  )
     &              +vy (i  ,j+1)
     &              +vy (i+1,j  )
     &              +vy (i+1,j+1))




c     area components at center of single cell volume

         rix = 0.5*(aix(i+1,j  )+aix(i,j))
         riy = 0.5*(aiy(i+1,j  )+aiy(i,j))

         rjx = 0.5*(ajx(i  ,j+1)+ajx(i,j))
         rjy = 0.5*(ajy(i  ,j+1)+ajy(i,j))


c        scalar product at center of volume of
c        (speed of sound + velocity) times area-components

c        normal to ax-areas
         ctx = rav*sqrt(rix*rix+riy*riy)

c        scalar product x-area*velocity
         cai = rix*vxv + riy*vyv

c        delta(time)/delta(volume) for i-direction
         rdti = abs(cai) + ctx


c        normal to ay-areas
         cty = rav*sqrt(rjx*rjx+rjy*rjy)

c        scalar product x-area*velocity
         caj = rjx*vxv + rjy*vyv

c        delta(time)/delta(volume) for j-direction
         rdtj = abs(caj) + cty

c        search for maximum reverse as minimum timestep

         rla(i,j) = amax1(rdti,rdtj)

       enddo
      enddo


c     values at boundary elements

      do  i = 1,im-1
c      lower wall
       rla(i,1 )   = 2.0*rla(i,1)
c      upper wall
       rla(i,jm-1) = 2.0*rla(i,jm-1)
      enddo

c     inlet and outlet
      do j = 1,jm-1
       rla(1   ,j) = 2.0*rla(1   ,j)
       rla(im-1,j) = 2.0*rla(im-1,j)
      enddo


c     delta(volume)/delta(time) at nodes for multiplying flux-balances

c     sum up for control volumes
      call sfl(0,rla,step)

c     invert to correct step values delta(time)/delta(volume)

      do  j = 1,jm
       do i = 1,im
        step(i,j) = 1.0/step(i,j)
       enddo
      enddo


      return
      end
      subroutine flbal(m)
c     ===================

      include 'param2D.h'
      include 'common2D.h'

c     calculate flux through area type ai

      do  j = 1,jm-1
       do i = 1,im

         flai(i,j) = (ar1(i,j  )
     &               +ar1(i,j+1))
     &               *aix(i,j)*0.5
     &              +(ar2(i,j  )
     &               +ar2(i,j+1))
     &               *aiy(i,j)*0.5
       enddo
      enddo

c     calculate flux through area type aj

      do  j = 1,jm
       do i = 1,im-1

         flaj(i,j) = (ar1(i  ,j)
     &               +ar1(i+1,j))
     &               *ajx(i,j)*0.5
     &              +(ar2(i  ,j)
     &               +ar2(i+1,j))
     &               *ajy(i,j)*0.5
       enddo
      enddo


c     boundary condition at walls j=1, j=jm

      if(m.eq.1.or.m.eq.4) then
       do i = 1,im-1

         flaj(i,1 ) = 0.0
         flaj(i,jm) = 0.0

       enddo
      endif

      if(m.eq.2) then

c     boundary condition at walls j=1, j=jm

       do i = 1,im-1
         flaj(i,1)  = (p(i  ,1)
     &                +p(i+1,1))
     &                *ajx(i,1)*0.5

         flaj(i,jm) = (p(i  ,jm)
     &                +p(i+1,jm))
     &                *ajx(i,jm)*0.5
       enddo
      endif

      if(m.eq.3) then

c     boundary condition at walls j=1, j=jm

       do i = 1,im-1
         flaj(i,1)  = (p(i  ,1)
     &                +p(i+1,1))
     &                *ajy(i,1)*0.5

         flaj(i,jm) = (p(i  ,jm)
     &                +p(i+1,jm))
     &                *ajy(i,jm)*0.5
       enddo
      endif


c     fluxbalance for single volume

      do  j = 1,jm-1
       do i = 1,im-1

         fblc(i,j) = flai(i,j)-flai(i+1,j  )
     &              +flaj(i,j)-flaj(i  ,j+1)
       enddo
      enddo

c     sum up for flux balances of node control volumes
      call sfl(0,fblc,fblp)


      return
      end
      subroutine geo
c     ==============

      include 'param2D.h'
      include 'common2D.h'

c     calculation of all area components and volumes

c     area components ai

      do  j = 1,jm-1
       do i = 1,im
        aix(i,j) = y(i,j+1) - y(i,j)
       enddo
      enddo

      do  j = 1,jm-1
       do i = 1,im
        aiy(i,j) = -(x(i,j+1) - x(i,j))
       enddo
      enddo


c     area components aj

      do  j = 1,jm
       do i = 1,im-1
        ajx(i,j) = -(y(i+1,j) - y(i,j))
       enddo
      enddo

      do  j = 1,jm
       do i = 1,im-1
        ajy(i,j) = x(i+1,j) - x(i,j)
       enddo
      enddo


c     calculate single volumes

      do  j = 1,jm-1
       do i = 1,im-1

          call gfl(x  (i  ,j  ),x  (i+1,j  ),
     &             x  (i+1,j+1),x  (i  ,j+1),
     &             y  (i  ,j  ),y  (i+1,j  ),
     &             y  (i+1,j+1),y  (i  ,j+1),
     &             volc(i,j))
       enddo
      enddo

c     sum up for node control volumes
      call sfl(0,volc,volp)


c      call write_geo



      return
      end
      subroutine gfl(x1,x2,x3,x4,y1,y2,y3,y4,a)
c     =========================================

c     Calculation of cross prduct of two vectors

      a = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)


      return
      end
      subroutine init
c     ===============

      include 'param2D.h'
      include 'common2D.h'

c     initialization variables

      nmax   = 0
      nstep  = 0
      rfd    = 0.0
      ft     = 0.0
      rg     = 0.0
      ermax  = 0.0
      rfpin  = 0.0
      rfpout = 0.0
      cpzucv = 0.0
      cp     = 0.0
      rgzucv = 0.0
      rgzucp = 0.0
      cvzurg = 0.0
      pt1    = 0.0
      tt1    = 0.0
      ptvbl  = 0.0
      pout   = 0.0

c     initialization arrays

      do  j = 1,jm
       do i = 1,im
         x   (i,j) = 0.0
         y   (i,j) = 0.0
         u1  (i,j) = 0.0
         u2  (i,j) = 0.0
         u3  (i,j) = 0.0
         u4  (i,j) = 0.0
         vx  (i,j) = 0.0
         vy  (i,j) = 0.0
         p   (i,j) = 0.0
         qc  (i,j) = 0.0
         u1e (i,j) = 0.0
         u2e (i,j) = 0.0
         u3e (i,j) = 0.0
         u4e (i,j) = 0.0
         u1m (i,j) = 0.0
         u2m (i,j) = 0.0
         u3m (i,j) = 0.0
         u4m (i,j) = 0.0

         step(i,j) = 0.0
         flai(i,j) = 0.0
         flaj(i,j) = 0.0
         rla (i,j) = 0.0
         fblc(i,j) = 0.0
         fblp(i,j) = 0.0
         aix (i,j) = 0.0
         aiy (i,j) = 0.0
         ajx (i,j) = 0.0
         ajy (i,j) = 0.0
         ar1 (i,j) = 0.0
         ar2 (i,j) = 0.0

        enddo
      enddo

      return
      end
      subroutine inlet
c     ================

      include 'param2D.h'
      include 'common2D.h'


c     boundary condition at inflow

      do j = 1,jm

c       from isentropic relations
        p(1,j)    = (1.0-rfpin)*p(2,j)+rfpin*pinold(j)
        pinold(j) = p(1,j)
        ts1       = tt1*(p(1,j)/pt1)**rgzucp
        rma       = sqrt(2.0*cvzurg*(tt1/ts1-1.0))
        u1(1,j)   = p(1,j)/(rg*ts1)
c       axial inflow
        vx(1,j)   = rma*sqrt(cpzucv*rg*ts1)
        vy(1,j)   = 0.0
        u2(1,j)   = u1(1,j)*vx(1,j)
        u3(1,j)   = u1(1,j)*vy(1,j)
        u4(1,j)   = p(1,j)/rgzucv+0.5*u1(1,j)*(vx(1,j)**2+vy(1,j)**2)

      enddo

      return
      end
      subroutine lax(q,fbal,qe,qm)
c     ============================


      include 'param2D.h'
      include 'common2D.h'

      dimension q(id,jd),qe(id,jd),qm(id,jd)
      dimension fbal(id,jd)


c     volume averages single volumes

      do  j = 1,jm-1
       do i = 1,im-1

         qe(i,j) = 0.25*(q(i  ,j  )+q(i  ,j+1)
     &                  +q(i+1,j+1)+q(i+1,j  ))

       enddo
      enddo


c     sum up volume averages to nodes
      call sfl(1,qe,qm)


c     update conservative variable
      do  j = 1,jm
       do i = 1,im

         q(i,j)  = (1.0-rfd)*q(i,j)+rfd*qm(i,j)
     &                +ft*fbal(i,j)*step(i,j)

       enddo
      enddo



      return
      end
      subroutine outlet
c     =================

      include 'param2D.h'
      include 'common2D.h'


c     boundary condition at outflow

       do j = 1,jm

        p(im,j) = (1.0-rfpout)*p(im,j)+rfpout*pout

       enddo

      return
      end
      subroutine post
c     ===============

      include 'param2D.h'
      include 'common2D.h'

      dimension flowi(id)


c     calculate massflow for all i-planes
      do  j = 1,jm
       do i = 1,im

          ar1(i,j) = u1(i,j)*vx(i,j)
          ar2(i,j) = u1(i,j)*vy(i,j)

       enddo
      enddo

      call flbal(1)

      write(11,*) '   '
      write(11,*) 'massflow ratio at i-planes: '
      write(11,*) '   '

      do   i = 1,im
       flowi(i) = 0.0
       do j = 1,jm-1
        flowi(i) = flowi(i) + flai(i,j)
       enddo
       if(i.eq.1)  flowin = flowi(i)
       if(i.eq.im) flowot = flowi(i)
       write(11,*) i,flowi(i),flowi(i)/flowin
      enddo

      write(11,*) '   '
      write(11,'(A17,F12.5)') 'massflow inlet  = ',flowin
      write(11,'(A17,F12.5)') 'massflow outlet = ',flowot
      write(11,*) '   '

c     calculate total pressure ratio

c     total pressure and temperature at inlet
      asin   = 0.0
      ptinlt = 0.0
      ttinlt = 0.0
      ptout  = 0.0
      ttout  = 0.0

      do j = 1,jm-1

c       inlet
        ts        = p(1,j)/(u1(1,j)*rg)
        aq        = cpzucv*rg*ts
        rmaq      = (vx(1,j)**2+vy(1,j)**2)/aq
        tt        = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
        pt        = p(1,j)/((ts/tt)**cpzurg)
        ptinlt    = ptinlt+pt*flai(1,j)
        ttinlt    = ttinlt+tt*flai(1,j)

c       outlet
        ts        = p(im,j)/(u1(im,j)*rg)
        aq        = cpzucv*rg*ts
        rmaq      = (vx(im,j)**2+vy(im,j)**2)/aq
        tt        = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
        pt        = p(im,j)/((ts/tt)**cpzurg)
        ptout     = ptout+pt*flai(im,j)
        ttout     = ttout+tt*flai(im,j)

      enddo

      ptinlt = ptinlt/flowin
      ttinlt = ttinlt/flowin

      ptout  = ptout/flowot
      ttout  = ttout/flowot

      ptv    = ptout/ptinlt
      ttv    = ttout/ttinlt

      write(11,*) '   '
      write(11,*) 'mass averaged total pressure and temperature:'
      write(11,*) '   '

      write(11,'(A17,F12.5)') 'pt inlet   = ',ptinlt
      write(11,'(A17,F12.5)') 'tt inlet   = ',ttinlt
      write(11,'(A17,F12.5)') 'pt outlet  = ',ptout
      write(11,'(A17,F12.5)') 'tt outlet  = ',ttout
      write(11,'(A17,F12.5)') 'ptout/ptin = ',ptv
      write(11,'(A17,F12.5)') 'ttout/ttin = ',ttv

      write( 6,'(A17,F12.5)') 'massflow inlet  = ',flowin
      write( 6,'(A17,F12.5)') 'massflow outlet = ',flowot
      write( 6,'(A17,F12.5)') 'pt inlet   = ',ptinlt
      write( 6,'(A17,F12.5)') 'tt inlet   = ',ttinlt
      write( 6,'(A17,F12.5)') 'pt outlet  = ',ptout
      write( 6,'(A17,F12.5)') 'tt outlet  = ',ttout
      write( 6,'(A17,F12.5)') 'ptout/ptin = ',ptv
      write( 6,'(A17,F12.5)') 'ttout/ttin = ',ttv

      write(11,*) '   '
      write(11,*) '   '
      write(11,*) '   '


      return
      end
      subroutine readst
c     =================

      include 'param2D.h'
      include 'common2D.h'

      character*80 titel


      namelist /l1/ nmax,nrfd,
     &              rfd,ft,rg,cpzucv,ermax,
     &              rfpin,rfpout,pt1,tt1,pout,ptvbl,
     &              pstart,u1star,vxstar,vystar,vzstar


c     data input
c     **********

      read(10,'(a80)') titel

      read(10,l1)

      write(11,'(A80)') titel
      write(11,'(A9,I5)')    'nmax   = ', nmax
      write(11,'(A9,I5)')    'nrfd   = ', nrfd
      write(11,'(A9,F10.5)') 'rg     = ', rg
      write(11,'(A9,F10.5)') 'cpzucv = ', cpzucv
      write(11,'(A9,F10.3)') 'pstart = ', pstart
      write(11,'(A9,F10.5)') 'vxstar = ', vxstar
      write(11,'(A9,F10.5)') 'vystar = ', vystar
      write(11,'(A9,F10.5)') 'u1star = ', u1star
      write(11,'(A9,F10.3)') 'pt1    = ', pt1
      write(11,'(A9,F10.3)') 'tt1    = ', tt1
      write(11,'(A9,F10.3)') 'pout   = ', pout
      write(11,'(A9,F10.5)') 'ft     = ', ft
      write(11,'(A9,F10.5)') 'rfd    = ', rfd
      write(11,'(A9,F10.5)') 'rfpin  = ', rfpin
      write(11,'(A9,F10.5)') 'rfpout = ', rfpout
      write(11,'(A9,F10.5)') 'ermax  = ', ermax

      write(6,'(A80)') titel
      write(6,'(A9,I5)')    'nmax   = ', nmax
      write(6,'(A9,I5)')    'nrfd   = ', nrfd
      write(6,'(A9,F10.5)') 'rg     = ', rg
      write(6,'(A9,F10.5)') 'cpzucv = ', cpzucv
      write(6,'(A9,F10.3)') 'pstart = ', pstart
      write(6,'(A9,F10.5)') 'vxstar = ', vxstar
      write(6,'(A9,F10.5)') 'vystar = ', vystar
      write(6,'(A9,F10.5)') 'u1star = ', u1star
      write(6,'(A9,F10.3)') 'pt1    = ', pt1
      write(6,'(A9,F10.3)') 'tt1    = ', tt1
      write(6,'(A9,F10.3)') 'pout   = ', pout
      write(6,'(A9,F10.5)') 'ft     = ', ft
      write(6,'(A9,F10.5)') 'rfd    = ', rfd
      write(6,'(A9,F10.5)') 'rfpin  = ', rfpin
      write(6,'(A9,F10.5)') 'rfpout = ', rfpout
      write(6,'(A9,F10.5)') 'ermax  = ', ermax

      
      return
      end
      subroutine read_coord
c     =====================

      include 'cgnslib_f.h'
      include 'param2D.h'
      include 'common2D.h'

      do  j = 1,jm
       do i = 1,im
        read(14,*) x(i,j),y(i,j)
       enddo
      enddo

      return
      end
      subroutine sfl(m,q,qm)
c     ======================

      include 'param2D.h'
      include 'common2D.h'

c     sum up single balances to control volume

      dimension q(id,jd),qm(id,jd)

      if(m.eq.0) then
       f4 = 1.0
       f2 = 1.0
      else
       f4 = 0.25
       f2 = 0.50
      endif

c     regular field
       do  j = 2,jm-1
        do i = 2,im-1

         qm(i,j)  = (q(i-1,j-1)
     &              +q(i-1,j  )
     &              +q(i  ,j  )
     &              +q(i  ,j-1))*f4
       enddo
      enddo

c     values at boundaries

      do  i = 2,im-1
c       lower wall
        qm(i,1 ) = (q(i-1,1)
     &             +q(i  ,1))*f2
c       upper wall
        qm(i,jm) = (q(i-1,jm-1)
     &             +q(i  ,jm-1))*f2
      enddo

      do  j = 2,jm-1
c       inlet
        qm(1,j)  = (q(1,j-1)
     &             +q(1,j  ))*f2
c       outlet
        qm(im,j) = (q(im-1,j-1)
     &             +q(im-1,j  ))*f2
      enddo


c     inlet/wall
      qm(1 ,1 ) = q(1 ,1   )
      qm(1 ,jm) = q(1 ,jm-1)


c     outlet/wall
      qm(im,1 ) = q(im-1,1   )
      qm(im,jm) = q(im-1,jm-1)




      return
      end
      subroutine start
c     ================

c     Initialization of variables

      include 'param2D.h'
      include 'common2D.h'

c     some gas relations

      cpzurg = cpzucv/(cpzucv-1.0)
      cvzurg = 1.0/(cpzucv-1.0)
      rgzucp = 1.0/cpzurg
      rgzucv = 1.0/cvzurg

c     initialization of density, pressure, velocity components
        do  j = 1,jm
         do i = 1,im

          u1(i,j) = u1star
          p (i,j) = pstart
          vx(i,j) = vxstar
          vy(i,j) = vystar

        enddo
       enddo

      do  j = 1,jm
       do i = 1,im

         u2(i,j) = u1(i,j)*vx(i,j)
         u3(i,j) = u1(i,j)*vy(i,j)
         u4(i,j) = p(i,j)/rgzucv+0.5*u1(i,j)*(vx(i,j)**2+vy(i,j)**2)

       enddo
      enddo


      return
      end
      subroutine update
c     =================

c     Iteration of conservative variables to convergence

      include 'param2D.h'
      include 'common2D.h'

      dimension flowi(id)

      rfdmin = rfd
      nconv  = 0

      do n = 1,nmax

       nstep = n

       if(n.lt.nrfd) then
        rfd = 1.0 - real(n)/real(nrfd)*(1.0-rfdmin)
       endif

c      calculate new timesteps
       call cfln


c      update density
c      **************

c      calculate flux variables

       do  j = 1,jm
        do i = 1,im

          ar1(i,j) = u1(i,j)*vx(i,j)
          ar2(i,j) = u1(i,j)*vy(i,j)

        enddo
       enddo

c      calculate flux balances
       call flbal(1)

c      update u1
       call lax(u1,fblp,u1e,u1m)


c      calculate massflow at inlet and outlet
       do   i = 1,im,im-1
        flowi(i) = 0.0
        do j = 1,jm-1
          flowi(i) = flowi(i) + flai(i,j)
        enddo
        if(i.eq.1)  flowin = flowi(i)
        if(i.eq.im) flowot = flowi(i)
        flrate = flowot/flowin
       enddo
       write(6 ,*)  'n = ',nstep,'massflow out/in',flrate,' rfd = ',rfd
       write(11,*)  'n = ',nstep,'massflow out/in',flrate,' rfd = ',rfd

       write(12,*)  real(n),flrate


c      update x-momentum
c      *****************

c      calculate flux variables

       do  j = 1,jm
        do i = 1,im

          ar1(i,j) = u2(i,j)*vx(i,j)+p(i,j)
          ar2(i,j) = u2(i,j)*vy(i,j)

        enddo
       enddo

c      calculate flux balances
       call flbal(2)

c      update u2
       call lax(u2,fblp,u2e,u2m)



c      update y-momentum
c      *****************

c      calculate flux variables

       do  j = 1,jm
        do i = 1,im

          ar1(i,j) = u3(i,j)*vx(i,j)
          ar2(i,j) = u3(i,j)*vy(i,j)+p(i,j)

        enddo
       enddo

c      calculate flux balances
       call flbal(3)

c      update u3
       call lax(u3,fblp,u3e,u3m)


c      update total energy
c      *******************

c      calculate flux variables

       do  j = 1,jm
        do i = 1,im

          ar1(i,j) = (u4(i,j)+p(i,j))*vx(i,j)
          ar2(i,j) = (u4(i,j)+p(i,j))*vy(i,j)

        enddo
       enddo


c      calculate flux balances
       call flbal(4)

c      update u4
       call lax(u4,fblp,u4e,u4m)


c      new primary variables

       do  j = 1,jm
        do i = 1,im

          qc(i,j) = amax1(0.00001,abs(vx(i,j)))

          vx(i,j) = u2(i,j)/u1(i,j)
          vy(i,j) = u3(i,j)/u1(i,j)

          p (i,j) = (u4(i,j) - 0.5*(vx(i,j)**2+vy(i,j)**2)
     &              *u1(i,j))*rgzucv

        enddo
      enddo


c      inflow boundary
       call inlet

c      outflow boundary
       call outlet


c      convergence control
       iconv = 0
       iemax = 0
       jemax = 0
       evmax = 0.0

       do  j = 1,jm
        do i = 1,im
          erij = 1.0-sqrt(vx(i,j)**2+vy(i,j)**2)/qc(i,j)
          if(abs(erij).gt.1000.0) then
           write(6 ,*) 'vx:',vx(i,j),'vy:',vy(i,j),'qc:',qc(i,j)
           write(6 ,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(11,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(6 ,*) 'erij:',erij,' Demo3D diverges!'
           write(11,*) 'erij:',erij,' Demo3D diverges!'
           goto 999
          endif

          if(erij.gt.evmax) then
           evmax = erij
           iemax = i
           jemax = j
          endif
        enddo
       enddo

       write(6 , *) 'n = ',n,' at',iemax,jemax,'evmax = ',evmax
       write(11, *) 'n = ',n,' at',iemax,jemax,'evmax = ',evmax
       write(13, *) evmax

c      check convergence
       if(abs(flrate-1.0).lt.ermax.and.n.gt.1000) then
        nconv = nconv +1
        write(6 ,*) 'nconv = ',nconv
        write(11,*) 'nconv = ',nconv
        if(nconv.eq.10) goto 999
       endif

      enddo

      write(11,*) nstep, ' '
      write(11,*) nstep, ' steps calculated'
      write(11,*) nstep, ' '

  999 continue


      return
      end
      subroutine write_boundary
c     =========================

      include 'cgnslib_f.h'
      include 'param2D.h'

      character*32 basename,zonename,famname

      dimension isize(3,3),ipnts(3,3)

      call cg_open_f('Demo2D_result.cgns',CG_MODE_MODIFY,index_file,ier)
      write(6,*) 'BC open: index_file =',index_file
      if(ier.ne.CG_OK) call cg_error_exit_f

      index_base = 1
      call cg_base_read_f(index_file,index_base,basename,
     &                   icelldim,iphysdim,ier)
      write(6,*) 'base: index_base =',index_base
      if(ier.ne.CG_OK) call cg_error_exit_f

      write(6,*) 'base: basename =',basename
      write(6,*) 'base: Celldim =',icelldim
      write(6,*) 'base: Physdim =',iphysdim

      index_zone = 1
      call cg_zone_read_f(index_file,index_base,index_zone,zonename,
     &                    isize,ier)
      if(ier.ne.CG_OK) call cg_error_exit_f

      write(6,*) 'zone: index_zone =',index_zone
      write(6,*) 'zone: zonename =',zonename

c     index ranges
      ilo = 1
      jlo = 1
      klo = 1
      ihi = isize(1,1)
      jhi = isize(2,1)
      khi = isize(3,1)

      write(6,*) 'ihi = ',ihi
      write(6,*) 'jhi = ',jhi
      write(6,*) 'khi = ',khi

c     create Inlet panel

      famname = 'Inlet'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family Inlet'
      write(6,*) 'Inlet: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

c     index space for Inlet
      ipnts(1,1) = ilo
      ipnts(2,1) = jlo
      ipnts(3,1) = klo
      ipnts(1,2) = ilo
      ipnts(2,2) = jhi
      ipnts(3,2) = khi

      call cg_boco_write_f(index_file,index_base,index_zone,'Inlet1',
     &             BCTunnelInflow,PointRange,2,ipnts,index_bc,ier)

c     write Outlet panel

      famname = 'Outlet'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family Outlet'
      write(6,*) 'Outlet: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f


c     index space for Outlet
      ipnts(1,1) = ihi
      ipnts(2,1) = jlo
      ipnts(3,1) = klo
      ipnts(1,2) = ihi
      ipnts(2,2) = jhi
      ipnts(3,2) = khi

      call cg_boco_write_f(index_file,index_base,index_zone,'Outlet',
     &             BCTunnelOutflow,PointRange,2,ipnts,index_bc,ier)


c     write WallHub panel

      famname = 'WallHub'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family WallHub'
      write(6,*) 'WallHub: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

c     index space for WallHub
      ipnts(1,1) = ilo
      ipnts(2,1) = jlo
      ipnts(3,1) = klo
      ipnts(1,2) = ihi
      ipnts(2,2) = jlo
      ipnts(3,2) = khi

      call cg_boco_write_f(index_file,index_base,index_zone,'WallHub',
     &             BCWallViscous,PointRange,2,ipnts,index_bc,ier)

c     write WallTip panel

      famname = 'WallTip'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family WallTip'
      write(6,*) 'WallTip: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

c     index space for WallTip
      ipnts(1,1) = ilo
      ipnts(2,1) = jhi
      ipnts(3,1) = klo
      ipnts(1,2) = ihi
      ipnts(2,2) = jhi
      ipnts(3,2) = khi

      call cg_boco_write_f(index_file,index_base,index_zone,'WallTip',
     &             BCWallViscous,PointRange,2,ipnts,index_bc,ier)

c     write WallLeft panel

      famname = 'WallLeft'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family WallLeft'
      write(6,*) 'WallLeft: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

c     index space for WallLeft
      ipnts(1,1) = ilo
      ipnts(2,1) = jlo
      ipnts(3,1) = klo
      ipnts(1,2) = ihi
      ipnts(2,2) = jhi
      ipnts(3,2) = klo

      call cg_boco_write_f(index_file,index_base,index_zone,'WallLeft',
     &             BCWallViscous,PointRange,2,ipnts,index_bc,ier)


c     write WallRight panel

      famname = 'WallRight'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family WallRight'
      write(6,*) 'WallLeft: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

c     index space for WallRight
      ipnts(1,1) = ilo
      ipnts(2,1) = jlo
      ipnts(3,1) = khi
      ipnts(1,2) = ihi
      ipnts(2,2) = jhi
      ipnts(3,2) = khi

      call cg_boco_write_f(index_file,index_base,index_zone,'WallRight',
     &             BCWallViscous,PointRange,2,ipnts,index_bc,ier)




c close CGNS file
      call cg_close_f(index_file,ier)
      write(6,'('' Successfully wrote boundary panels to cgns'')')

      return
      end
      subroutine write_cgns
c     =====================

      include 'cgnslib_f.h'
      include 'param2D.h'
      include 'common2D.h'

      dimension x3d (id,jd,kd),y3d  (id,jd,kd),z3d (id,jd,kd)
      dimension vx3d(id,jd,kd),vy3d (id,jd,kd),vz3d(id,jd,kd)
      dimension p3d (id,jd,kd),u13d (id,jd,kd)

      dimension isize(3,3)

      character basename*32,zonename*32,solname*32

      km = 3

      do   k = 1,km
       do  j = 1,jm
        do i = 1,im

         x3d (i,j,k) = x(i,j)
         y3d (i,j,k) = y(i,j)
         z3d (i,j,k) = real(k)*0.001
         vx3d(i,j,k) = vx(i,j)
         vy3d(i,j,k) = vy(i,j)
         vz3d(i,j,k) = 0.0
         p3d (i,j,k) = p(i,j)
         u13d(i,j,k) = u1(i,j)

        enddo
       enddo
      enddo


c     WRITE X, Y, Z GRID POINTS TO CGNS FILE
c     open CGNS file for write

      call cg_open_f('Demo2D_result.cgns',CG_MODE_WRITE,index_file,ier)

      write(6,*) 'index_file =',index_file

      if (ier .ne. CG_OK) call cg_error_exit_f

c create base (user can give any name)

      basename = 'Base1'
      icelldim = 3
      iphysdim = 3

      call cg_base_write_f(index_file,basename,icelldim,iphysdim,
     &                     index_base,ier)

      write(6,*) 'index_base =',index_base
      write(6,*) 'icelldim   =',icelldim
      write(6,*) 'iphysdim   =',iphysdim


c define zone name (user can give any name)

      zonename = 'Block 1'

c vertex size
      isize(1,1) = im
      isize(2,1) = jm
      isize(3,1) = km

c cell size
      isize(1,2) = isize(1,1)-1
      isize(2,2) = isize(2,1)-1
      isize(3,2) = isize(3,1)-1

c boundary vertex size (always zero for structured grids)
      isize(1,3)=0
      isize(2,3)=0
      isize(3,3)=0


c create zone
      call cg_zone_write_f(index_file,index_base,zonename,isize,
     &                     Structured,index_zone,ier)

      write(6,*) 'index_zone =',index_zone

c write grid coordinates (user must use SIDS-standard names here)
      call cg_coord_write_f(index_file,index_base,index_zone,
     &                     RealSingle,'CoordinateX',x3d,index_coord,ier)
      write(6,*) 'index_coord X =',index_coord

      call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateY',y3d,index_coord,ier)
      write(6,*) 'index_coord Y =',index_coord

      call cg_coord_write_f(index_file,index_base,index_zone,
     &                     RealSingle,'CoordinateZ',z3d,index_coord,ier)
      write(6,*) 'index_coord Z =',index_coord


c    define flow solution node name (user can give any name)
      solname = 'FlowSolution'

c    create flow solution node
      call cg_sol_write_f(index_file,index_base,index_zone,solname,
     &                    Vertex,index_flow,ier)

c    write flow solution (user must use SIDS-standard names here)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Density',u13d,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Pressure',p3d,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityX',vx3d,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityY',vy3d,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityZ',vz3d,index_field,ier)



c close CGNS file
      call cg_close_f(index_file,ier)
      write(6,'('' Successfully wrote grid to file Demo2D_result.cgns
     & '')')



      return
      end
      subroutine write_geo
c     ====================

      include 'param2D.h'
      include 'common2D.h'


       do  j = 1,jm-1
        do i = 1,im-1

         ar1(i,j)   =  aix(i+1,j  )-aix(i,j)
     &               + ajx(i  ,j+1)-ajx(i,j)
     &               + aiy(i+1,j  )-aiy(i,j)
     &               + ajy(i  ,j+1)-ajy(i,j)
       enddo
      enddo

      write(11,*) 'area aix,aiy'
      do  i = 1,im
       do j = 1,jm-1
         write(11,*) i,j,aix(i,j),aiy(i,j)
       enddo
      enddo

      write(11,*) 'area ajx,ajy'
      do  i = 1,im-1
       do j = 1,jm
         write(11,*) i,j,ajx(i,j),ajy(i,j)
       enddo
      enddo

      write(11,*) 'volc'
      do  i = 1,im-1
       do j = 1,jm-1
         write(11,*) i,j,volc(i,j)
       enddo
      enddo

      write(11,*) 'area check'

      do  j = 1,jm-1
       do i = 1,im-1
         write(11,*) i,j,ar1(i,j)
       enddo
      enddo

      return
      end
      subroutine write_prime
c     ======================

      include 'param2D.h'
      include 'common2D.h'


      write(11,*) 'write prime','im=',im,'jm=',jm
c      do  j = 1,jm
          j = 61
       do i = 1,im,4
         write(11,*) i,j,vx(i,j),vy(i,j),p(i,j),u1(i,j)
       enddo
c      enddo

      write(11,*) 'write conservatives'

c      do  j = 1,jm
          j = 61
       do i = 1,im,4
         write(11,*) i,j,u1(i,j),u2(i,j),u3(i,j),u4(i,j)
       enddo
c      enddo

      return
      end
