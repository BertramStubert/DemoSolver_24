c     fortran code of a Demo-1D-Solver (euler method)

c     Author: Bertram Stubert
c     Time of Project-Start: 15.03.2023

c     Goal of this code is a demonstration for interested people to get
c     an idea of a 1D-aero-code functionality as simple as possible

c     !!! There is no warranty due to the code !!!
c     !!! This code is written for demonstration purposes only !!!
c     !!! Do not deliver any results of this code for any projects !!!
      
      
      include 'param1D.h'
      include 'common1D.h'

      open (unit = 10, file = 'Demo1D_input.dat')
      open (unit = 11, file = 'Demo1D_output.dat')
      open (unit = 13, file = 'Demo1D_evmax.dat')
      open (unit = 15, file = 'Demo1D_vx.plot')
      open (unit = 16, file = 'Demo1D_ps.plot')
      open (unit = 17, file = 'Demo1D_ma.plot')
      open (unit = 18, file = 'Demo1D_konti.dat')

c     konwn indices
      im = id


c     Mainprogram
c     ***********

c----------------------------------------------------------------------
c     phase 1 - prepare all data for iteration in time

c     initialization variables and arrays
c     ***********************************
      call init

c     parameter input with namelist l1
c     ********************************
      call readst

c     read areas aix
c     **************
      call read_area

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

      

      write(11,*) ' Demo3D ended'
      write(11,*) ' ************'


      call cpu_time(time)

      write(11,*) ' '
      write(11,'(a13,f10.2)') 'cpu-seconds: ', time
      write(11,*) ' '

      close (10)
      close (11)
      close (13)
      close (15)
      close (16)
      close (17)
      close (18)

      stop
      end
      subroutine cfln
c     ===============

      include 'param1D.h'
      include 'common1D.h'

c     calculate local timesteps
c     *************************

c     product of (velocity + speed of sound) * area at center
c     of volume gives a term like delta(volume)/delta(timestep).
c     The inverse value is the required faktor for multiplying
c     with the fluxbalance to get time changes for conservatives

c     calculate speed of sound
      do i = 1,im
        aq     = amax1(100.0,cpzucv*p(i)/u1(i))
        ar1(i) = sqrt(aq)
      enddo


       do i = 1,im-1

c        volume mean values for speed of sound and velocity

         rav = 0.5*(ar1(i+1)+ar1(i))

         vxv = 0.5*(vx (i+1)+vx(i))


c     area components at center of single cell volume

         rix = 0.5*(aix(i+1)+aix(i))

c        product at center of volume of
c        (speed of sound + velocity) times area-components

c        normal to ai-areas
         ctx = rav*rix

c        product x-area*velocity
         cai = rix*vxv

c        delta(time)/delta(volume) for i-direction
         rdti = abs(cai) + ctx

         rla(i) = rdti

      enddo


c     values at boundary elements

c     inlet and outlet
      rla(1   ) = 2.0*rla(1   )
      rla(im-1) = 2.0*rla(im-1)


c     delta(volume)/delta(time) at nodes for multiplying flux-balances

c     sum up for control volumes
      call sfl(0,rla,step)


c     invert to correct step values delta(time)/delta(volume)
      do i = 1,im
        step(i) = 1.0/step(i)
      enddo


      return
      end
      subroutine flbal(m)
c     ===================

      include 'param1D.h'
      include 'common1D.h'

c     calculate flux through area type ai

      do i = 1,im
         flai(i) =  ar1(i)*aix(i)
      enddo

c     fluxbalance for single volume

      if(m.eq.2) then

       do i = 1,im-1
         fblc(i) =  flai(i) - flai(i+1)
     &            + 0.5*(p(i+1)+p(i))*(aix(i+1)-aix(i))
       enddo

      else

       do i = 1,im-1
         fblc(i) = flai(i)-flai(i+1)
       enddo
      endif


c     sum up for flux balances of node control volumes
      call sfl(0,fblc,fblp)


      return
      end
      subroutine init
c     ===============

      include 'param1D.h'
      include 'common1D.h'

c     initialization variables

      nmax   = 0
      nstep  = 0
      istart = 0
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
      pout   = 0.0
      pinold = 0.0

c     initialization arrays

       do i = 1,im
         x   (i) = 0.0
         u1  (i) = 0.0
         u2  (i) = 0.0
         u3  (i) = 0.0
         vx  (i) = 0.0
         p   (i) = 0.0
         qc  (i) = 0.0
         u1e (i) = 0.0
         u2e (i) = 0.0
         u3e (i) = 0.0
         u1m (i) = 0.0
         u2m (i) = 0.0
         u3m (i) = 0.0

         step(i) = 0.0
         flai(i) = 0.0
         rla (i) = 0.0
         fblc(i) = 0.0
         fblp(i) = 0.0
         aix (i) = 0.0
         ar1 (i) = 0.0
       enddo

      return
      end
      subroutine inlet
c     ================

      include 'param1D.h'
      include 'common1D.h'


c     boundary condition at inflow

c       from isentropic relations
        p(1)    = (1.0-rfpin)*p(2)+rfpin*pinold
        pinold  = p(1)
        ts1     = tt1*(p(1)/pt1)**rgzucp
        rma     = sqrt(2.0*cvzurg*(tt1/ts1-1.0))
        u1(1)   = p(1)/(rg*ts1)
c       axial inflow
        vx(1)   = rma*sqrt(cpzucv*rg*ts1)
        u2(1)   = u1(1)*vx(1)
        u3(1)   = p(1)/rgzucv+0.5*u1(1)*vx(1)**2

      return
      end
      subroutine lax(q,fbal,qe,qm)
c     ============================


      include 'param1D.h'
      include 'common1D.h'

      dimension q(id),qe(id),qm(id)
      dimension fbal(id)


c     volume averages single volumes

      do i = 1,im-1

         qe(i) = 0.5*(q(i+1)+q(i))

      enddo


c     sum up volume averages to nodes
      call sfl(1,qe,qm)


c     update conservative variable
      do i = 1,im
         q(i)  = (1.0-rfd)*q(i)+rfd*qm(i) + ft*fbal(i)*step(i)
      enddo


      return
      end
      subroutine outlet
c     =================

      include 'param1D.h'
      include 'common1D.h'


c     boundary condition at outflow

      p(im) = (1.0-rfpout)*p(im)+rfpout*pout


      return
      end
      subroutine post
c     ===============

      include 'param1D.h'
      include 'common1D.h'

c     ------------------------------------
c      probable test output
c      write(11,*) ' '
c      write(11,*) 'post vx,vy,vz,p,u1'
c      call write_prime
c     ------------------------------------

c     calculate massflow for all i-planes
      do i = 1,im
          ar1(i) = u1(i)*vx(i)
      enddo

      call flbal(1)

      write(11,*) 'post flai(1),flai(im):',flai(1),flai(im)
      write(11,*) '   '
      write(11,*) 'massflow ratio at i-planes: '
      write(11,*) '   '
      do  i = 1,im
       write(11,*) i,flai(i),flai(i)/flai(1)
      enddo

      write(11,*) '   '
      write(11,'(A18,F12.5)') 'massflow inlet   = ',flai(1)
      write(11,'(A18,F12.5)') 'massflow outlet  = ',flai(im)
      write(11,'(A18,F12.5)') 'massratio out/in = ',flai(im)/flai(1)
      write(11,*) '   '

c     calculate total pressure ratio

c     total pressure and temperature at inlet

c     inlet
      ts        = p(1)/(u1(1)*rg)
      aq        = cpzucv*rg*ts
      rmaq      = vx(1)**2/aq
      tt        = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
      pt        = p(1)/((ts/tt)**cpzurg)
      ptinlt    = pt
      ttinlt    = tt

c     outlet
      ts        = p(im)/(u1(im)*rg)
      aq        = cpzucv*rg*ts
      rmaq      = vx(im)**2/aq
      tt        = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
      pt        = p(im)/((ts/tt)**cpzurg)
      ptout     = pt
      ttout     = tt

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

      write( 6,'(A17,F12.5)') 'massflow inlet  = ',flai(1)
      write( 6,'(A17,F12.5)') 'massflow outlet = ',flai(im)
      write( 6,'(A17,F12.5)') 'pt inlet   = ',ptinlt
      write( 6,'(A17,F12.5)') 'tt inlet   = ',ttinlt
      write( 6,'(A17,F12.5)') 'pt outlet  = ',ptout
      write( 6,'(A17,F12.5)') 'tt outlet  = ',ttout
      write( 6,'(A17,F12.5)') 'ptout/ptin = ',ptv
      write( 6,'(A17,F12.5)') 'ttout/ttin = ',ttv

      write(11,*) '   '
      write(11,*) '   '
      write(11,*) '   '

      do i = 1,im
       write(15,*) vx(i)
       write(16,*) p(i)
       write(17,*) vx(i)/sqrt(cpzucv*p(i)/u1(i))
      enddo


      return
      end
      subroutine readst
c     =================

      include 'param1D.h'
      include 'common1D.h'

      character*80 titel


      namelist /l1/ nmax,nrfd,
     &              rfd,ft,rg,cpzucv,ermax,
     &              rfpin,rfpout,pt1,tt1,pout,ptvbl,
     &              pstart,u1star,vxstar


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
      subroutine read_area
c     ====================

      include 'param1D.h'
      include 'common1D.h'

      do  i = 1,im
        aix(i) = 1.0
      enddo

      return
      end
      subroutine sfl(m,dudt,dqdt)
c     ===========================

      include 'param1D.h'
      include 'common1D.h'

c     sum up single balances to control volume

      dimension dudt(id),dqdt(id)

      if(m.eq.1) then
       f2 = 0.5
      else
       f2 = 1.0
      endif

c     regular field
      do i = 2,im-1
         dqdt(i)  = (dudt(i-1)
     &              +dudt(i  ))*f2
      enddo

c     values at boundaries

c     inlet
      dqdt(1 ) = dudt(1)
c     outlet
      dqdt(im) = dudt(im-1)


      return
      end
      subroutine start
c     ================

c     Initialization of variables

      include 'param1D.h'
      include 'common1D.h'


c     some gas relations

      cpzurg = cpzucv/(cpzucv-1.0)
      cvzurg = 1.0/(cpzucv-1.0)
      rgzucp = 1.0/cpzurg
      rgzucv = 1.0/cvzurg

c     initialization of density, pressure, velocity components
      do i = 1,im
          u1(i) = u1star
          p (i) = pstart
          vx(i) = vxstar
      enddo

      do i = 1,im
         u2(i) = u1(i)*vx(i)
         u3(i) = p(i)/rgzucv+0.5*u1(i)*vx(i)**2
      enddo

      pinold = p(1)


      return
      end
      subroutine update
c     =================

c     Iteration of conservative variables to convergence

      include 'param1D.h'
      include 'common1D.h'

      dimension flowi(id)

      rfdmin = rfd

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

       do i = 1,im
         ar1(i) = u1(i)*vx(i)
       enddo

c      calculate flux balances
       call flbal(1)

c      update u1
       call lax(u1,fblp,u1e,u1m)


c      claculate massflow at inlet and outlet

       flrate = flai(im)/flai(1)

       write(6 ,*)  'n = ',nstep,'massflow out/in',flrate,' rfd = ',rfd
       write(11,*)  'n = ',nstep,'massflow out/in',flrate,' rfd = ',rfd

c      output for visualization
       write(18,*)  real(n),flrate


c      update x-momentum
c      *****************

c      calculate flux variables

       do i = 1,im
          ar1(i) = u2(i)*vx(i)+p(i)
       enddo

c      calculate flux balances
       call flbal(2)

c      update u2
       call lax(u2,fblp,u2e,u2m)


c      update total energy
c      *******************

c      calculate flux variables

       do i = 1,im
          ar1(i) = (u3(i)+p(i))*vx(i)
       enddo


c      calculate flux balances
       call flbal(3)

c      update u4
       call lax(u3,fblp,u3e,u3m)


c      new primary variables

       do i = 1,im

          qc(i) = amax1(0.00001,abs(vx(i)))

          vx(i) = u2(i)/u1(i)

          p (i) = (u3(i)-0.5*u1(i)*vx(i)**2)*rgzucv
     &
       enddo


c      inflow boundary
       call inlet

c      outflow boundary
       call outlet


c      convergence control
       iconv = 0
       iemax = 0
       evmax = 0.0

       do i = 1,im
          eri = 1.0-vx(i)/qc(i)
          if(abs(eri).gt.1000.0) then
           write(6 ,*) 'vx:',vx(i),'qc:',qc(i)
           write(6 ,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(11,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(6 ,*) 'erij:',erij,' Demo3D diverges!'
           write(11,*) 'erij:',erij,' Demo3D diverges!'
           goto 999
          endif

          if(erij.gt.evmax) then
           evmax = eri
           iemax = i
          endif
       enddo

c      possible output
c       write(6 , *) 'n = ',n,' at',iemax,'evmax = ',evmax
c       write(11, *) 'n = ',n,' at',iemax,'evmax = ',evmax
       write(13, *) evmax

c      check convergence
       if(abs(flrate-1.0).lt.ermax.and.n.gt.1000) goto 999

      enddo

  999 continue

      write(11,*) ' '
      write(11,*) 'calculated steps: ',nstep
      write(11,*) ' '

      return
      end
      subroutine write_prime
c     ======================

      include 'param1D.h'
      include 'common1D.h'


      write(11,*) 'write prime vx,p,rho,ma'

      do i = 1,im

       a   = sqrt(cpzucv*p(i)/u1(i))
       write(11,*) i,vx(i),p(i),u1(i),vx(i)/a
     &
      enddo

c      write(11,*) 'write conservatives'

c      do i = 1,im
c         write(11,*) i,u2(i),u3(i)
c      enddo

      return
      end
