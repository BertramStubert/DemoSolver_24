
c     fortran code of a Demo-2D-Solver (euler method)

c     Author: Bertram Stubert
c     Time of Project-Start: 21.02.2024

c     Goal of this code is a demonstration for interested people to get
c     an idea of a 1D-aero-code functionality as simple as possible

c     Application deals with the calculation of a
c     3D-one-block-configuration building a channel with curved walls.

c     !!! There is no warranty due to the code !!!
c     !!! This code is written for demonstration purposes only !!!
c     !!! Do not deliver any results of this code for any projects !!!

      
      
      include 'param3D.h'
      include 'common3D.h'

      open (unit = 10, file = 'Demo3D_input.dat')
      open (unit = 11, file = 'Demo3D_output.dat')
      open (unit = 12, file = 'Demo3D_konti.dat')
      open (unit = 13, file = 'Demo3D_evmax.dat')
      open (unit = 14, file = 'Demo3D_stop.dat')


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

c     read cgns file
c     ***************
      call read_cgns

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
      

      write(11,*) ' Demo3D ended'
      write(11,*) ' ************'


      call cpu_time(time)

      write(11,*) ' '
      write(11,'(a13,f10.2)') 'cpu-hours  : ', time/3600.0
      write(11,'(a13,f10.2)') 'cpu-minutes: ', time/60.0
      write(11,'(a13,f10.2)') 'cpu-seconds: ', time
      write(11,*) ' '


      close (10)
      close (11)
      close (12)
      close (13)

      stop
      end
      subroutine cfln
c     ===============

      include 'param3D.h'
      include 'common3D.h'

c     calculate local timesteps
c     *************************

c     Scalar product of (velocity + speed of sound) * area at center
c     of volume gives a term like delta(volume)/delta(timestep).
c     The inverse step is the required faktor for multiplying
c     with the fluxbalance

c     calculate speed of sound
      do   k = 1,km
       do  j = 1,jm
        do i = 1,im
         aq         = amax1(100.0,cpzucv*p(i,j,k)/u1(i,j,k))
         ar1(i,j,k) = sqrt(aq)
        enddo
       enddo
      enddo

      do   k = 1,km-1
       do  j = 1,jm-1
        do i = 1,im-1

c        volume mean values for speed of sound and velocity

         rav = 0.125*(ar1(i  ,j  ,k  )+ar1(i  ,j+1,k  )
     &               +ar1(i  ,j+1,k+1)+ar1(i  ,j  ,k+1)
     &               +ar1(i+1,j  ,k  )+ar1(i+1,j+1,k  )
     &               +ar1(i+1,j+1,k+1)+ar1(i+1,j  ,k+1))

         vxv = 0.125*(vx (i  ,j  ,k  )+vx (i  ,j+1,k  )
     &               +vx (i  ,j+1,k+1)+vx (i  ,j  ,k+1)
     &               +vx (i+1,j  ,k  )+vx (i+1,j+1,k  )
     &               +vx (i+1,j+1,k+1)+vx (i+1,j  ,k+1))

         vyv = 0.125*(vy (i  ,j  ,k  )+vy (i  ,j+1,k  )
     &               +vy (i  ,j+1,k+1)+vy (i  ,j  ,k+1)
     &               +vy (i+1,j  ,k  )+vy (i+1,j+1,k  )
     &               +vy (i+1,j+1,k+1)+vy (i+1,j  ,k+1))

         vzv = 0.125*(vz (i  ,j  ,k  )+vz (i  ,j+1,k  )
     &               +vz (i  ,j+1,k+1)+vz (i  ,j  ,k+1)
     &               +vz (i+1,j  ,k  )+vz (i+1,j+1,k  )
     &               +vz (i+1,j+1,k+1)+vz (i+1,j  ,k+1))


c     area components at center of single cell volume

         rix = 0.5*(aix(i+1,j  ,k  )+aix(i,j,k))
         riy = 0.5*(aiy(i+1,j  ,k  )+aiy(i,j,k))
         riz = 0.5*(aiz(i+1,j  ,k  )+aiz(i,j,k))
         rjx = 0.5*(ajx(i  ,j+1,k  )+ajx(i,j,k))
         rjy = 0.5*(ajy(i  ,j+1,k  )+ajy(i,j,k))
         rjz = 0.5*(ajz(i  ,j+1,k  )+ajz(i,j,k))
         rkx = 0.5*(akx(i  ,j  ,k+1)+akx(i,j,k))
         rky = 0.5*(aky(i  ,j  ,k+1)+aky(i,j,k))
         rkz = 0.5*(akz(i  ,j  ,k+1)+akz(i,j,k))

c        scalar product at center of volume of
c        (speed of sound + velocity) times area-components

c        normal to ax-areas
         ctx = rav*sqrt(rix*rix+riy*riy+riz*riz)

c        scalar product x-area*velocity
         cai = rix*vxv + riy*vyv + riz*vzv

c        delta(time)/delta(volume) for i-direction
         rdti = abs(cai) + ctx


c        normal to ay-areas
         cty = rav*sqrt(rjx*rjx+rjy*rjy+rjz*rjz)

c        scalar product x-area*velocity
         caj = rjx*vxv + rjy*vyv + rjz*vzv

c        delta(time)/delta(volume) for j-direction
         rdtj = abs(caj) + cty


c        normal to az-areas
         ctz = rav*sqrt(rkx*rkx+rky*rky+rkz*rkz)

c        scalar product z-area*velocity
         cak = rkx*vxv + rky*vyv + rkz*vzv

c        delta(time)/delta(volume) for k-direction
         rdtk = abs(cak) + ctz

c        search for maximum reverse as minimum timestep

         rla(i,j,k) = amax1(rdti,rdtj,rdtk)

        enddo
       enddo
      enddo


c     values at boundary elements

      do   j = 1,jm-1
       do  i = 1,im-1
c       left wall
        rla(i,j,1)    = 2.0*rla(i,j,1)
c       right wall
        rla(i,j,km-1) = 2.0*rla(i,j,km-1)
       enddo
      enddo

      do   k = 1,km-1
       do  i = 1,im-1
c       lower wall
        rla(i,1 ,k)   = 2.0*rla(i,1,k)
c       upper wall
        rla(i,jm-1,k) = 2.0*rla(i,jm-1,k)
       enddo
      enddo

c     inlet and outlet
      do  k = 1,km-1
       do j = 1,im-1
        rla(1   ,j,k) = 2.0*rla(1   ,j,k)
        rla(im-1,j,k) = 2.0*rla(im-1,j,k)
       enddo
      enddo


c     delta(volume)/delta(time) at nodes for multiplying flux-balances

c     sum up for control volumes
      call sfl(0,rla,step)

c     invert to correct step values delta(time)/delta(volume)

      do   k = 1,km
       do  j = 1,jm
        do i = 1,im
         step(i,j,k) = 1.0/step(i,j,k)
        enddo
       enddo
      enddo


      return
      end
      subroutine flbal(m)
c     ===================

      include 'param3D.h'
      include 'common3D.h'

c     calculate flux through area type ai

      do   k = 1,km-1
       do  j = 1,jm-1
        do i = 1,im

         flai(i,j,k) = (ar1(i,j  ,k  )
     &                 +ar1(i,j+1,k  )
     &                 +ar1(i,j+1,k+1)
     &                 +ar1(i,j  ,k+1))
     &                 *aix(i,j,k)*0.25
     &                +(ar2(i,j  ,k  )
     &                 +ar2(i,j+1,k  )
     &                 +ar2(i,j+1,k+1)
     &                 +ar2(i,j  ,k+1))
     &                 *aiy(i,j,k)*0.25
     &                +(ar3(i,j  ,k  )
     &                 +ar3(i,j+1,k  )
     &                 +ar3(i,j+1,k+1)
     &                 +ar3(i,j  ,k+1))
     &                 *aiz(i,j,k)*0.25
        enddo
       enddo
      enddo

c     calculate flux through area type aj

      do   k = 1,km-1
       do  j = 1,jm
        do i = 1,im-1

         flaj(i,j,k) = (ar1(i  ,j,k  )
     &                 +ar1(i  ,j,k+1)
     &                 +ar1(i+1,j,k+1)
     &                 +ar1(i+1,j,k  ))
     &                 *ajx(i,j,k)*0.25
     &                +(ar2(i  ,j,k  )
     &                 +ar2(i  ,j,k+1)
     &                 +ar2(i+1,j,k+1)
     &                 +ar2(i+1,j,k  ))
     &                 *ajy(i,j,k)*0.25
     &                +(ar3(i  ,j,k  )
     &                 +ar3(i  ,j,k+1)
     &                 +ar3(i+1,j,k+1)
     &                 +ar3(i+1,j,k  ))
     &                 *ajz(i,j,k)*0.25
        enddo
       enddo
      enddo

c     calculate flux through area type ak

      do   k = 1,km
       do  j = 1,jm-1
        do i = 1,im-1

         flak(i,j,k) = (ar1(i  ,j  ,k)
     &                 +ar1(i+1,j  ,k)
     &                 +ar1(i+1,j+1,k)
     &                 +ar1(i  ,j+1,k))
     &                 *akx(i,j,k)*0.25
     &                +(ar2(i  ,j  ,k)
     &                 +ar2(i+1,j  ,k)
     &                 +ar2(i+1,j+1,k)
     &                 +ar2(i  ,j+1,k))
     &                 *aky(i,j,k)*0.25
     &                +(ar3(i  ,j  ,k)
     &                 +ar3(i+1,j  ,k)
     &                 +ar3(i+1,j+1,k)
     &                 +ar3(i  ,j+1,k))
     &                 *akz(i,j,k)*0.25
        enddo
       enddo
      enddo

c     boundary condition at walls j=1, j=jm

      if(m.eq.1.or.m.eq.5) then
       do  k = 1,km-1
        do i = 1,im-1

         flaj(i,1 ,k) = 0.0
         flaj(i,jm,k) = 0.0

        enddo
       enddo

c     boundary condition at walls k=1, k=km

       do  j = 1,jm-1
        do i = 1,im-1

         flak(i,j,1)  = 0.0
         flak(i,j,km) = 0.0

        enddo
       enddo
      endif

      if(m.eq.2) then

c     boundary condition at walls j=1, j=jm

       do  k = 1,km-1
        do i = 1,im-1
         flaj(i,1,k)  = (p(i  ,1 ,k  )
     &                  +p(i  ,1 ,k+1)
     &                  +p(i+1,1 ,k+1)
     &                  +p(i+1,1 ,k  ))
     &                  *ajx(i,1 ,k)*0.25

         flaj(i,jm,k) = (p(i  ,jm,k  )
     &                  +p(i  ,jm,k+1)
     &                  +p(i+1,jm,k+1)
     &                  +p(i+1,jm,k  ))
     &                  *ajx(i,jm,k)*0.25

        enddo
       enddo

c     boundary condition at walls k=1, k=km

       do  j = 1,jm-1
        do i = 1,im-1
         flak(i,j,1)  = (p(i  ,j  ,1)
     &                  +p(i+1,j  ,1)
     &                  +p(i+1,j+1,1)
     &                  +p(i  ,j+1,1))
     &                  *akx(i,j,1)*0.25

         flak(i,j,km) = (p(i  ,j  ,km)
     &                  +p(i+1,j  ,km)
     &                  +p(i+1,j+1,km)
     &                  +p(i  ,j+1,km))
     &                  *akx(i,j,km)*0.25
        enddo
       enddo
      endif

      if(m.eq.3) then

c     boundary condition at walls j=1, j=jm

       do  k = 1,km-1
        do i = 1,im-1
         flaj(i,1,k)  = (p(i  ,1 ,k  )
     &                  +p(i  ,1 ,k+1)
     &                  +p(i+1,1 ,k+1)
     &                  +p(i+1,1 ,k  ))
     &                  *ajy(i,1 ,k)*0.25

         flaj(i,jm,k) = (p(i  ,jm,k  )
     &                  +p(i  ,jm,k+1)
     &                  +p(i+1,jm,k+1)
     &                  +p(i+1,jm,k  ))
     &                  *ajy(i,jm,k)*0.25

        enddo
       enddo

c     boundary condition at walls k=1, k=km

       do  j = 1,jm-1
        do i = 1,im-1
         flak(i,j,1)  = (p(i  ,j  ,1)
     &                  +p(i+1,j  ,1)
     &                  +p(i+1,j+1,1)
     &                  +p(i  ,j+1,1))
     &                  *aky(i,j,1)*0.25

         flak(i,j,km) = (p(i  ,j  ,km)
     &                  +p(i+1,j  ,km)
     &                  +p(i+1,j+1,km)
     &                  +p(i  ,j+1,km))
     &                  *aky(i,j,km)*0.25
        enddo
       enddo
      endif

      if(m.eq.4) then

c     boundary condition at walls j=1, j=jm

       do  k = 1,km-1
        do i = 1,im-1
         flaj(i,1,k)  = (p(i  ,1 ,k  )
     &                  +p(i  ,1 ,k+1)
     &                  +p(i+1,1 ,k+1)
     &                  +p(i+1,1 ,k  ))
     &                  *ajz(i,1 ,k)*0.25

         flaj(i,jm,k) = (p(i  ,jm,k  )
     &                  +p(i  ,jm,k+1)
     &                  +p(i+1,jm,k+1)
     &                  +p(i+1,jm,k  ))
     &                  *ajz(i,jm,k)*0.25

        enddo
       enddo

c     boundary condition at walls k=1, k=km

       do  j = 1,jm-1
        do i = 1,im-1
         flak(i,j,1)  = (p(i  ,j  ,1)
     &                  +p(i+1,j  ,1)
     &                  +p(i+1,j+1,1)
     &                  +p(i  ,j+1,1))
     &                  *akz(i,j,1)*0.25

         flak(i,j,km) = (p(i  ,j  ,km)
     &                  +p(i+1,j  ,km)
     &                  +p(i+1,j+1,km)
     &                  +p(i  ,j+1,km))
     &                  *akz(i,j,km)*0.25
        enddo
       enddo
      endif

c     fluxbalance for single volume

      do   k = 1,km-1
       do  j = 1,jm-1
        do i = 1,im-1

         fblc(i,j,k) = flai(i,j,k)-flai(i+1,j  ,k  )
     &                +flaj(i,j,k)-flaj(i  ,j+1,k  )
     &                +flak(i,j,k)-flak(i  ,j  ,k+1)
        enddo
       enddo
      enddo

c     sum up for flux balances of node control volumes
      call sfl(0,fblc,fblp)


      return
      end
      subroutine geo
c     ==============

      include 'param3D.h'
      include 'common3D.h'

c     calculation of all area components and volumes

c     area components ai

      do   k = 1,km-1
       do  j = 1,jm-1
        do i = 1,im

          call gfl(y  (i,j  ,k  ),y  (i,j+1,k  ),
     &             y  (i,j+1,k+1),y  (i,j  ,k+1),
     &             z  (i,j  ,k  ),z  (i,j+1,k  ),
     &             z  (i,j+1,k+1),z  (i,j  ,k+1),
     &             aix(i,j,k))
        enddo
       enddo
      enddo

      do   k = 1,km-1
       do  j = 1,jm-1
        do i = 1,im

          call gfl(z  (i,j  ,k  ),z  (i,j+1,k  ),
     &             z  (i,j+1,k+1),z  (i,j  ,k+1),
     &             x  (i,j  ,k  ),x  (i,j+1,k  ),
     &             x  (i,j+1,k+1),x  (i,j  ,k+1),
     &             aiy(i,j,k))
        enddo
       enddo
      enddo


      do   k = 1,km-1
       do  j = 1,jm-1
        do i = 1,im

          call gfl(x  (i,j  ,k  ),x  (i,j+1,k  ),
     &             x  (i,j+1,k+1),x  (i,j  ,k+1),
     &             y  (i,j  ,k  ),y  (i,j+1,k  ),
     &             y  (i,j+1,k+1),y  (i,j  ,k+1),
     &             aiz(i,j,k))
        enddo
       enddo
      enddo

c     area components aj

      do   k = 1,km-1
       do  j = 1,jm
        do i = 1,im-1

          call gfl(y  (i  ,j,k  ),y  (i  ,j,k+1),
     &             y  (i+1,j,k+1),y  (i+1,j,k  ),
     &             z  (i  ,j,k  ),z  (i  ,j,k+1),
     &             z  (i+1,j,k+1),z  (i+1,j,k  ),
     &             ajx(i,j,k))
        enddo
       enddo
      enddo

      do   k = 1,km-1
       do  j = 1,jm
        do i = 1,im-1

          call gfl(z  (i  ,j,k  ),z  (i  ,j,k+1),
     &             z  (i+1,j,k+1),z  (i+1,j,k  ),
     &             x  (i  ,j,k  ),x  (i  ,j,k+1),
     &             x  (i+1,j,k+1),x  (i+1,j,k  ),
     &             ajy(i,j,k))
        enddo
       enddo
      enddo

      do   k = 1,km-1
       do  j = 1,jm
        do i = 1,im-1

          call gfl(x  (i  ,j,k  ),x  (i  ,j,k+1),
     &             x  (i+1,j,k+1),x  (i+1,j,k  ),
     &             y  (i  ,j,k  ),y  (i  ,j,k+1),
     &             y  (i+1,j,k+1),y  (i+1,j,k  ),
     &             ajz(i,j,k))
        enddo
       enddo
      enddo

c     area components ak

      do   k = 1,km
       do  j = 1,jm-1
        do i = 1,im-1

          call gfl(y  (i  ,j  ,k),y  (i+1,j  ,k),
     &             y  (i+1,j+1,k),y  (i  ,j+1,k),
     &             z  (i  ,j  ,k),z  (i+1,j  ,k),
     &             z  (i+1,j+1,k),z  (i  ,j+1,k),
     &             akx(i,j,k))
        enddo
       enddo
      enddo

      do   k = 1,km
       do  j = 1,jm-1
        do i = 1,im-1

          call gfl(z  (i  ,j  ,k),z  (i+1,j  ,k),
     &             z  (i+1,j+1,k),z  (i  ,j+1,k),
     &             x  (i  ,j  ,k),x  (i+1,j  ,k),
     &             x  (i+1,j+1,k),x  (i  ,j+1,k),
     &             aky(i,j,k))
        enddo
       enddo
      enddo

      do   k = 1,km
       do  j = 1,jm-1
        do i = 1,im-1

          call gfl(x  (i  ,j  ,k),x  (i+1,j  ,k),
     &             x  (i+1,j+1,k),x  (i  ,j+1,k),
     &             y  (i  ,j  ,k),y  (i+1,j  ,k),
     &             y  (i+1,j+1,k),y  (i  ,j+1,k),
     &             akz(i,j,k))
        enddo
       enddo
      enddo

c     calculate single volumes

      do   k = 1,km-1
       do  j = 1,jm-1
        do i = 1,im-1

         call voltet(x(i  ,j  ,k  ),x(i  ,j+1,k  ),
     &               x(i  ,j+1,k+1),x(i  ,j  ,k+1),
     &               x(i+1,j  ,k  ),x(i+1,j+1,k  ),
     &               x(i+1,j+1,k+1),x(i+1,j  ,k+1),
     &               y(i  ,j  ,k  ),y(i  ,j+1,k  ),
     &               y(i  ,j+1,k+1),y(i  ,j  ,k+1),
     &               y(i+1,j  ,k  ),y(i+1,j+1,k  ),
     &               y(i+1,j+1,k+1),y(i+1,j  ,k+1),
     &               z(i  ,j  ,k  ),z(i  ,j+1,k  ),
     &               z(i  ,j+1,k+1),z(i  ,j  ,k+1),
     &               z(i+1,j  ,k  ),z(i+1,j+1,k  ),
     &               z(i+1,j+1,k+1),z(i+1,j  ,k+1),
     &               volc(i,j,k))

        enddo
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

      include 'param3D.h'
      include 'common3D.h'

c     initialization variables

      im     = 0
      jm     = 0
      km     = 0
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
      ptvbl  = 0.0
      pout   = 0.0

c     initialization arrays

      do   k = 1,km
       do  j = 1,jm
        do i = 1,im
         x   (i,j,k) = 0.0
         y   (i,j,k) = 0.0
         z   (i,j,k) = 0.0
         u1  (i,j,k) = 0.0
         u2  (i,j,k) = 0.0
         u3  (i,j,k) = 0.0
         u4  (i,j,k) = 0.0
         u5  (i,j,k) = 0.0
         vx  (i,j,k) = 0.0
         vy  (i,j,k) = 0.0
         vz  (i,j,k) = 0.0
         p   (i,j,k) = 0.0
         qc  (i,j,k) = 0.0
         u1e (i,j,k) = 0.0
         u2e (i,j,k) = 0.0
         u3e (i,j,k) = 0.0
         u4e (i,j,k) = 0.0
         u5e (i,j,k) = 0.0
         u1m (i,j,k) = 0.0
         u2m (i,j,k) = 0.0
         u3m (i,j,k) = 0.0
         u4m (i,j,k) = 0.0
         u5m (i,j,k) = 0.0

         step(i,j,k) = 0.0
         flai(i,j,k) = 0.0
         flaj(i,j,k) = 0.0
         flak(i,j,k) = 0.0
         rla (i,j,k) = 0.0
         fblc(i,j,k) = 0.0
         fblp(i,j,k) = 0.0
         aix (i,j,k) = 0.0
         aiy (i,j,k) = 0.0
         aiz (i,j,k) = 0.0
         ajx (i,j,k) = 0.0
         ajy (i,j,k) = 0.0
         ajz (i,j,k) = 0.0
         akx (i,j,k) = 0.0
         aky (i,j,k) = 0.0
         akz (i,j,k) = 0.0
         ar1 (i,j,k) = 0.0
         ar2 (i,j,k) = 0.0
         ar3 (i,j,k) = 0.0

        enddo

        pinold(j,k) = 0.0

       enddo
      enddo

      return
      end
      subroutine inlet
c     ================

      include 'param3D.h'
      include 'common3D.h'

      dimension ain(jd,kd)

      asin   = 0.0
      pinsum = 0.0

      do  k = 1,km-1
       do j = 1,jm-1

        ain(j,k)  = sqrt(aix(1,j,k)**2+aiy(1,j,k)**2+aiz(1,j,k)**2)
        asin      = asin+ain(j,k)
        pinsum    = pinsum + p(2,j,k)*ain(j,k)
       enddo
      enddo

      pinsum = pinsum/asin

      write( 6,*) 'mean static pressure inlet: ',pinsum
      write(11,*) 'mean static pressure inlet: ',pinsum

c     boundary condition at inflow

      do  k = 1,km
       do j = 1,jm

c       from isentropic relations
        p(1,j,k)    = (1.0-rfpin)*p(2,j,k)+rfpin*pinold(j,k)
        pinold(j,k) = p(1,j,k)
        ts1         = tt1*(p(1,j,k)/ptin(j))**rgzucp
        rma         = sqrt(2.0*cvzurg*(tt1/ts1-1.0))
        u1(1,j,k)   = p(1,j,k)/(rg*ts1)
c       axial inflow
        vx(1,j,k)   = rma*sqrt(cpzucv*rg*ts1)
        vy(1,j,k)   = 0.0
        vz(1,j,k)   = 0.0
        u2(1,j,k)   = u1(1,j,k)*vx(1,j,k)
        u3(1,j,k)   = u1(1,j,k)*vy(1,j,k)
        u4(1,j,k)   = u1(1,j,k)*vz(1,j,k)
        u5(1,j,k)   = p(1,j,k)/rgzucv + u1(1,j,k)*
     &                0.5*(vx(1,j,k)**2+vy(1,j,k)**2+vz(1,j,k)**2)

       enddo
      enddo

      return
      end
      subroutine lax(q,fbal,qe,qm)
c     ============================


      include 'param3D.h'
      include 'common3D.h'

      dimension q(id,jd,kd),qe(id,jd,kd),qm(id,jd,kd)
      dimension fbal(id,jd,kd)


c     volume averages single volumes

      do   k = 1,km-1
       do  j = 1,jm-1
        do i = 1,im-1

         qe(i,j,k) = 0.125*(q(i  ,j  ,k  )+q(i  ,j+1,k  )
     &                     +q(i  ,j+1,k+1)+q(i  ,j  ,k+1)
     &                     +q(i+1,j  ,k  )+q(i+1,j+1,k  )
     &                     +q(i+1,j+1,k+1)+q(i+1,j  ,k+1))


        enddo
       enddo
      enddo


c     sum up volume averages to nodes
      call sfl(1,qe,qm)


c     update conservative variable
      do   k = 1,km
       do  j = 1,jm
        do i = 1,im

         q(i,j,k)  = (1.0-rfd)*q(i,j,k)+rfd*qm(i,j,k)
     &                  +ft*fbal(i,j,k)*step(i,j,k)

        enddo
       enddo
      enddo



      return
      end
      subroutine outlet
c     =================

      include 'param3D.h'
      include 'common3D.h'

      dimension aout(jd,kd)

c     boundary condition at outflow

       do  k = 1,km
        do j = 1,jm

        p(im,j,k) = (1.0-rfpout)*p(im,j,k)+rfpout*pout

        enddo
       enddo


      asout  = 0.0
      pimsum = 0.0

      do  k = 1,km-1
       do j = 1,jm-1

        aout(j,k)  = sqrt(aix(im,j,k)**2+aiy(im,j,k)**2+aiz(im,j,k)**2)
        asout      = asout+aout(j,k)

        pimsum = pimsum+0.25*(p(im,j  ,k  )
     &                       +p(im,j+1,k  )
     &                       +p(im,j+1,k+1)
     &                       +p(im,j  ,k+1))*aout(j,k)

       enddo
      enddo

      pimsum   = pimsum/asout

      write(6,*)  'mean static pressure outlet:',pimsum,' pout:',pout
      write(11,*) 'mean static pressure outlet:',pimsum,' pout:',pout


      return
      end
      subroutine post
c     ===============

      include 'param3D.h'
      include 'common3D.h'

      dimension flowi(id)

      write(11,*) ' '
      write(11,*) 'update vx,vy,vz,p,u1'
      call write_prime


c     calculate massflow for all i-planes
       do   k = 1,km
        do  j = 1,jm
         do i = 1,im

          ar1(i,j,k) = u1(i,j,k)*vx(i,j,k)
          ar2(i,j,k) = u1(i,j,k)*vy(i,j,k)
          ar3(i,j,k) = u1(i,j,k)*vz(i,j,k)

         enddo
        enddo
       enddo

      call flbal(1)

      write(11,*) '   '
      write(11,*) 'massflow ratio at i-planes: '
      write(11,*) '   '

      do   i = 1,im
       flowi(i) = 0.0
       do  k = 1,km-1
        do j = 1,jm-1
         flowi(i) = flowi(i) + flai(i,j,k)
        enddo
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

      do  k = 1,km-1
       do j = 1,jm-1

c       inlet
        ts        = p(1,j,k)/(u1(1,j,k)*rg)
        aq        = cpzucv*rg*ts
        rmaq      = (vx(1,j,k)**2+vy(1,j,k)**2+vz(1,j,k)**2)/aq
        tt        = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
        pt        = p(1,j,k)/((ts/tt)**cpzurg)
        ptinlt    = ptinlt+pt*flai(1,j,k)
        ttinlt    = ttinlt+tt*flai(1,j,k)

c       outlet
        ts        = p(im,j,k)/(u1(im,j,k)*rg)
        aq        = cpzucv*rg*ts
        rmaq      = (vx(im,j,k)**2+vy(im,j,k)**2+vz(im,j,k)**2)/aq
        tt        = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
        pt        = p(im,j,k)/((ts/tt)**cpzurg)
        ptout     = ptout+pt*flai(im,j,k)
        ttout     = ttout+tt*flai(im,j,k)

       enddo
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

      include 'param3D.h'     
      include 'common3D.h'

      character*80 titel


      namelist /l1/ nmax,nrfd,istart,
     &              rfd,ft,rg,cpzucv,ermax,
     &              rfpin,rfpout,pt1,tt1,pout,ptvbl,
     &              pstart,u1star,vxstar,vystar,vzstar


c     data input
c     **********

      read(10,'(a80)') titel

      read(10,l1)

      write(11,'(A80)') titel
      write(11,'(A9,I5)')    'istart = ', istart
      write(11,'(A9,I5)')    'nmax   = ', nmax
      write(11,'(A9,I5)')    'nrfd   = ', nrfd
      write(11,'(A9,F10.5)') 'rg     = ', rg
      write(11,'(A9,F10.5)') 'cpzucv = ', cpzucv
      write(11,'(A9,F10.3)') 'pstart = ', pstart
      write(11,'(A9,F10.5)') 'vxstar = ', vxstar
      write(11,'(A9,F10.5)') 'vystar = ', vystar
      write(11,'(A9,F10.5)') 'vzstar = ', vzstar
      write(11,'(A9,F10.5)') 'u1star = ', u1star
      write(11,'(A9,F10.3)') 'pt1    = ', pt1
      write(11,'(A9,F10.3)') 'ptvbl  = ', ptvbl
      write(11,'(A9,F10.3)') 'tt1    = ', tt1
      write(11,'(A9,F10.3)') 'pout   = ', pout
      write(11,'(A9,F10.5)') 'ft     = ', ft
      write(11,'(A9,F10.5)') 'rfd    = ', rfd
      write(11,'(A9,F10.5)') 'rfpin  = ', rfpin
      write(11,'(A9,F10.5)') 'rfpout = ', rfpout
      write(11,'(A9,F10.5)') 'ermax  = ', ermax

      write(6,'(A80)') titel
      write(6,'(A9,I5)')    'istart = ', istart
      write(6,'(A9,I5)')    'nmax   = ', nmax
      write(6,'(A9,I5)')    'nrfd   = ', nrfd
      write(6,'(A9,F10.5)') 'rg     = ', rg
      write(6,'(A9,F10.5)') 'cpzucv = ', cpzucv
      write(6,'(A9,F10.3)') 'pstart = ', pstart
      write(6,'(A9,F10.5)') 'vxstar = ', vxstar
      write(6,'(A9,F10.5)') 'vystar = ', vystar
      write(6,'(A9,F10.5)') 'vzstar = ', vzstar
      write(6,'(A9,F10.5)') 'u1star = ', u1star
      write(6,'(A9,F10.3)') 'pt1    = ', pt1
      write(6,'(A9,F10.3)') 'ptvbl  = ', ptvbl
      write(6,'(A9,F10.3)') 'tt1    = ', tt1
      write(6,'(A9,F10.3)') 'pout   = ', pout
      write(6,'(A9,F10.5)') 'ft     = ', ft
      write(6,'(A9,F10.5)') 'rfd    = ', rfd
      write(6,'(A9,F10.5)') 'rfpin  = ', rfpin
      write(6,'(A9,F10.5)') 'rfpout = ', rfpout
      write(6,'(A9,F10.5)') 'ermax  = ', ermax

      
      return
      end
      subroutine read_cgns
c     ====================

      include 'cgnslib_f.h'
      include 'param3D.h'
      include 'common3D.h'

      character*32 basename,zonename
      dimension isize(3,3),irmin(3),irmax(3)

      write(6,*) 'read_cgns:'

      call cg_open_f('Demo3D_result.cgns',CG_MODE_READ,index_file,ier)
      write(6,*) 'read open: index_file =',index_file
      if(ier.ne.CG_OK) call cg_error_exit_f

      index_base = 1
      call cg_base_read_f(index_file,index_base,basename,
     &                   icelldim,iphysdim,ier)
      write(6,*) 'read base: index_base =',index_base
      if(ier.ne.CG_OK) call cg_error_exit_f

      write(6,*) 'read base: basename = ',basename
      write(6,*) 'read base: Celldim =',icelldim
      write(6,*) 'read base: Physdim =',iphysdim

      index_zone = 1
      call cg_zone_read_f(index_file,index_base,index_zone,zonename,
     &                    isize,ier)
      if(ier.ne.CG_OK) call cg_error_exit_f

      write(6,*) 'read zone: index_zone =',index_zone
      write(6,*) 'read zone: zonename =',zonename

c     index ranges
      irmin(1) = 1
      irmin(2) = 1
      irmin(3) = 1
      irmax(1) = isize(1,1)
      irmax(2) = isize(2,1)
      irmax(3) = isize(3,1)

      write(6,*) 'irmax(1) = ',irmax(1)
      write(6,*) 'irmax(2) = ',irmax(2)
      write(6,*) 'irmax(3) = ',irmax(3)

c     read cgns file
      call cg_coord_read_f(index_file,index_base,index_zone,
     &                    'CoordinateX',RealSingle,irmin,irmax,x,ier)
      call cg_coord_read_f(index_file,index_base,index_zone,
     &                    'CoordinateY',RealSingle,irmin,irmax,y,ier)
      call cg_coord_read_f(index_file,index_base,index_zone,
     &                    'CoordinateZ',RealSingle,irmin,irmax,z,ier)

      if(istart.eq.1) then

c      we know there is only one FlowSolution_t
       index_flow = 1

       call cg_field_read_f(index_file,index_base,index_zone,index_flow,
     &                     'Density',RealSingle,irmin,irmax,u1,ier)

       call cg_field_read_f(index_file,index_base,index_zone,index_flow,
     &                     'Pressure',RealSingle,irmin,irmax,p,ier)

       call cg_field_read_f(index_file,index_base,index_zone,index_flow,
     &                     'VelocityX',RealSingle,irmin,irmax,vx,ier)

       call cg_field_read_f(index_file,index_base,index_zone,index_flow,
     &                     'VelocityY',RealSingle,irmin,irmax,vy,ier)

       call cg_field_read_f(index_file,index_base,index_zone,index_flow,
     &                     'VelocityZ',RealSingle,irmin,irmax,vz,ier)

      endif


c close CGNS file
      call cg_close_f(index_file,ier)
      write(6,'('' Successfully read Demo3D_start.cgns'')')

      im = irmax(1)
      jm = irmax(2)
      km = irmax(3)


      return
      end
      subroutine sfl(m,q,qm)
c     ======================

      include 'param3D.h'
      include 'common3D.h'

c     sum up single balances to control volume

      dimension q(id,jd,kd),qm(id,jd,kd)

      if(m.eq.0) then
       f8 = 1.0
       f4 = 1.0
       f2 = 1.0
      else
       f8 = 0.125
       f4 = 0.25
       f2 = 0.50
      endif


c     regular field
      do   k = 2,km-1
       do  j = 2,jm-1
        DO i = 2,im-1

         qm(i,j,k)  =   (q(i-1,j-1,k-1)
     &                  +q(i-1,j-1,k  )
     &                  +q(i  ,j-1,k  )
     &                  +q(i  ,j-1,k-1)
     &                  +q(i-1,j  ,k-1)
     &                  +q(i-1,j  ,k  )
     &                  +q(i  ,j  ,k  )
     &                  +q(i  ,j  ,k-1))
     &                  *f8
        enddo
       enddo
      enddo

c     values at boundaries

      do   j = 2,jm-1
       do  i = 2,im-1
c       left wall
        qm(i,j,1)  = (q(i-1,j-1,1)
     &               +q(i-1,j  ,1)
     &               +q(i  ,j  ,1)
     &               +q(i  ,j-1,1))
     &               *f4
c       right wall
        qm(i,j,km) = (q(i-1,j-1,km-1)
     &               +q(i-1,j  ,km-1)
     &               +q(i  ,j  ,km-1)
     &               +q(i  ,j-1,km-1))
     &               *f4
       enddo
      enddo

      do   k = 2,km-1
       do  i = 2,im-1
c       lower wall
        qm(i,1 ,k) = (q(i-1,1 ,k-1)
     &               +q(i-1,1 ,k  )
     &               +q(i  ,1 ,k  )
     &               +q(i  ,1 ,k-1))
     &               *f4
c       upper wall
        qm(i,jm,k) = (q(i-1,jm-1,k-1)
     &               +q(i-1,jm-1,k  )
     &               +q(i  ,jm-1,k  )
     &               +q(i  ,jm-1,k-1))
     &               *f4
       enddo
      enddo

      do   k = 2,km-1
       do  j = 2,jm-1
c       inlet
        qm(1,j,k)  =  (q(1,j-1,k-1)
     &                +q(1,j-1,k  )
     &                +q(1,j  ,k  )
     &                +q(1,j  ,k-1))
     &                *f4
c       outlet
        qm(im,j,k) =  (q(im-1,j-1,k-1)
     &                +q(im-1,j-1,k  )
     &                +q(im-1,j  ,k  )
     &                +q(im-1,j  ,k-1))
     &                *f4
       enddo
      enddo


      do i = 2,im-1
c      lower left corner
       qm(i,1 ,1 )  =  (q(i-1,1,1)
     &                 +q(i  ,1,1))*f2

c      upper left corner
       qm(i,1 ,km)  =  (q(i-1,1,km-1)
     &                 +q(i  ,1,km-1))*f2

c      upper right corner
       qm(i,jm,km)  =  (q(i-1,jm-1,km-1)
     &                 +q(i  ,jm-1,km-1))*f2

c      lower right corner
       qm(i,jm,1 )  =  (q(i-1,jm-1,1)
     &                 +q(i  ,jm-1,1))*f2
      enddo

      do k = 2,km-1
c      inlet lower
       qm(1 ,1,k)  =  (q(1,1,k-1)
     &                +q(1,1,k  ))*f2
c      outlet lower
       qm(im,1,k)  =  (q(im-1,1,k-1)
     &                +q(im-1,1,k  ))*f2
     &
c      inlet upper
       qm(1 ,jm,k) =  (q(1,jm-1,k-1)
     &                +q(1,jm-1,k  ))*f2
c      outlet upper
       qm(im,jm,k) =  (q(im-1,jm-1,k-1)
     &                +q(im-1,jm-1,k  ))*f2
      enddo

      do j = 2,jm-1
c      inlet left corner
       qm(1 ,j,1)  =  (q(1,j-1,1)
     &                +q(1,j  ,1))*f2
c      outlet left corner
       qm(im,j,1)  =  (q(im-1,j-1,1)
     &                +q(im-1,j  ,1))*f2
c      inlet right corner
       qm(1 ,j,km) =  (q(1,j-1,km-1)
     &                +q(1,j  ,km-1))*f2
c      outlet right corner
       qm(im,j,km) =  (q(im-1,j-1,km-1)
     &                +q(im-1,j  ,km-1))*f2
      enddo

c     corners inlet
      qm(1 ,1 ,1 ) = q(1 ,1   ,1   )
      qm(1 ,jm,1 ) = q(1 ,jm-1,1   )
      qm(1 ,jm,km) = q(1 ,jm-1,km-1)
      qm(1 ,1 ,km) = q(1 ,1   ,km-1)

c     corners outlet
      qm(im,1 ,1 ) = q(im-1,1   ,1   )
      qm(im,jm,1 ) = q(im-1,jm-1,1   )
      qm(im,jm,km) = q(im-1,jm-1,km-1)
      qm(im,1 ,km) = q(im-1,1   ,km-1)


      return
      end
      subroutine start
c     ================

c     Initialization of variables

      include 'param3D.h'
      include 'common3D.h'

      write(11,*) 'ptin profile'

      pi = 4.*atan(1.)

c     inlet total pressure profile
      do j = 1,jm
       ptin(j) = ptvbl*pt1
     &   +((1.0-ptvbl)*pt1)*sin(real(j-1)/real(jm-1)*pi)
       write(11,*)  j,ptin(j)
      enddo

c     some gas relations

      cpzurg = cpzucv/(cpzucv-1.0)
      cvzurg = 1.0/(cpzucv-1.0)
      rgzucp = 1.0/cpzurg
      rgzucv = 1.0/cvzurg

c     initialization of density, pressure, velocity components
      if(istart.eq.0) then
       do   k = 1,km
        do  j = 1,jm
         do i = 1,im

          u1(i,j,k) = u1star
          p (i,j,k) = pstart
          vx(i,j,k) = vxstar
          vy(i,j,k) = vystar
          vz(i,j,k) = vzstar

         enddo
        enddo
       enddo
      endif

      do   k = 1,km
       do  j = 1,jm
        do i = 1,im

         u2(i,j,k) = u1(i,j,k)*vx(i,j,k)
         u3(i,j,k) = u1(i,j,k)*vy(i,j,k)
         u4(i,j,k) = u1(i,j,k)*vz(i,j,k)
         u5(i,j,k) = p(i,j,k)/rgzucv + 0.5*u1(i,j,k)*
     &               (vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2)
        enddo

         pinold(j,k) = p (1 ,j,k)

       enddo
      enddo

c      write(11,*) 'start:'
c      call write_prime

      return
      end
      subroutine update
c     =================

c     Iteration of conservative variables to convergence

      include 'param3D.h'
      include 'common3D.h'

      dimension flowi(id)

      rfdmin = rfd
      nconv  = 0

      do n = 1,nmax

       nstep = n

       if(istart.eq.0.and.n.lt.nrfd) then
        rfd = 1.0 - real(n)/real(nrfd)*(1.0-rfdmin)
       endif

c      calculate new timesteps
       call cfln


c      update density
c      **************

c      calculate flux variables

       do   k = 1,km
        do  j = 1,jm
         do i = 1,im

          ar1(i,j,k) = u1(i,j,k)*vx(i,j,k)
          ar2(i,j,k) = u1(i,j,k)*vy(i,j,k)
          ar3(i,j,k) = u1(i,j,k)*vz(i,j,k)

         enddo
        enddo
       enddo

c      calculate flux balances
       call flbal(1)

c      update u1
       call lax(u1,fblp,u1e,u1m)


c      claculate massflow at inlet and outlet
       do   i = 1,im,im-1
        flowi(i) = 0.0
        do  k = 1,km-1
         do j = 1,jm-1
          flowi(i) = flowi(i) + flai(i,j,k)
         enddo
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

       do   k = 1,km
        do  j = 1,jm
         do i = 1,im

          ar1(i,j,k) = u2(i,j,k)*vx(i,j,k)+p(i,j,k)
          ar2(i,j,k) = u2(i,j,k)*vy(i,j,k)
          ar3(i,j,k) = u2(i,j,k)*vz(i,j,k)

         enddo
        enddo
       enddo

c      calculate flux balances
       call flbal(2)

c      update u2
       call lax(u2,fblp,u2e,u2m)



c      update y-momentum
c      *****************

c      calculate flux variables

       do   k = 1,km
        do  j = 1,jm
         do i = 1,im

          ar1(i,j,k) = u3(i,j,k)*vx(i,j,k)
          ar2(i,j,k) = u3(i,j,k)*vy(i,j,k)+p(i,j,k)
          ar3(i,j,k) = u3(i,j,k)*vz(i,j,k)

         enddo
        enddo
       enddo

c      calculate flux balances
       call flbal(3)

c      update u3
       call lax(u3,fblp,u3e,u3m)



c      update z-momentum
c      *****************

c      calculate flux variables

       do   k = 1,km
        do  j = 1,jm
         do i = 1,im

          ar1(i,j,k) = u4(i,j,k)*vx(i,j,k)
          ar2(i,j,k) = u4(i,j,k)*vy(i,j,k)
          ar3(i,j,k) = u4(i,j,k)*vz(i,j,k)+p(i,j,k)

         enddo
        enddo
       enddo

c      calculate flux balances
       call flbal(4)

c      update u4
       call lax(u4,fblp,u4e,u4m)


c      update total energy
c      *******************

c      calculate flux variables

       do   k = 1,km
        do  j = 1,jm
         do i = 1,im

          ar1(i,j,k) = (u5(i,j,k)+p(i,j,k))*vx(i,j,k)
          ar2(i,j,k) = (u5(i,j,k)+p(i,j,k))*vy(i,j,k)
          ar3(i,j,k) = (u5(i,j,k)+p(i,j,k))*vz(i,j,k)

         enddo
        enddo
       enddo


c      calculate flux balances
       call flbal(5)

c      update u5
       call lax(u5,fblp,u5e,u5m)


c      new primary variables

       do   k = 1,km
        do  j = 1,jm
         do i = 1,im

          qc(i,j,k) = amax1(0.00001,vx(i,j,k))

          vx(i,j,k) = u2(i,j,k)/u1(i,j,k)
          vy(i,j,k) = u3(i,j,k)/u1(i,j,k)
          vz(i,j,k) = u4(i,j,k)/u1(i,j,k)

          p (i,j,k) = (u5(i,j,k)
     &                -0.5*(vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2)
     &                *u1(i,j,k))*rgzucv
         enddo

         pinold(j,k) = p (1 ,j,k)

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
       kemax = 0
       evmax = 0.0

       do   k = 1,km
        do  j = 1,jm
         do i = 1,im
          erijk = 1.0-vx(i,j,k)/qc(i,j,k)
          if(abs(erijk).gt.1000.0) then
           write(6 ,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(11,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(6 ,*) 'erijk:',erijk,' Demo3D diverges!'
           write(11,*) 'erijk:',erijk,' Demo3D diverges!'
           goto 999
          endif

          if(erijk.gt.evmax) then
           evmax = erijk
           iemax = i
           jemax = j
           kemax = k
          endif
         enddo
        enddo
       enddo

       write(6 , *) 'n = ',n,' at',iemax,jemax,kemax,'evmax = ',evmax
       write(11, *) 'n = ',n,' at',iemax,jemax,kemax,'evmax = ',evmax
       write(13, *) evmax

c      check convergence

       if(abs(flrate-1.0).lt.ermax.and.n.gt.1000) then
        nconv = nconv +1
        write(6 ,*) 'nconv = ',nconv
        write(11,*) 'nconv = ',nconv
        if(nconv.eq.10) goto 999
       endif

       read(14,*) istop
       if(istop.eq.1) then
        write(11,*) ' '
        write(11,*) 'Run has been stopped by user'
        write(11,*) ' '
        goto 999
       endif
       rewind 14

      enddo


  999 continue


      return
      end
      subroutine write_cgns
c     =====================

      include 'cgnslib_f.h'
      include 'param3D.h'
      include 'common3D.h'

      dimension isize(3,3)

      character basename*32,zonename*32,solname*32

      write(6,*) 'write_cgns:'

c     Write results to cgns file
c     open CGNS file for write

      call cg_open_f('Demo3D_result.cgns',CG_MODE_MODIFY,index_file,ier)
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


c     write grid coordinates (user must use SIDS-standard names here)
      call cg_coord_write_f(index_file,index_base,index_zone,
     &                     RealSingle,'CoordinateX',x,index_coord,ier)
      write(6,*) 'index_coord X =',index_coord

      call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateY',y,index_coord,ier)
      write(6,*) 'index_coord Y =',index_coord

      call cg_coord_write_f(index_file,index_base,index_zone,
     &                     RealSingle,'CoordinateZ',z,index_coord,ier)
      write(6,*) 'index_coord Z =',index_coord


c    define flow solution node name (user can give any name)
      solname = 'FlowSolution'

c    create flow solution node
      call cg_sol_write_f(index_file,index_base,index_zone,solname,
     &                    Vertex,index_flow,ier)

c    write flow solution (user must use SIDS-standard names here)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Density',u1,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Pressure',p,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityX',vx,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityY',vy,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityZ',vz,index_field,ier)



c close CGNS file
      call cg_close_f(index_file,ier)
      write(6,'(''Demo3D_result.cgns written'')')



      return
      end
      subroutine write_geo
c     ====================

      include 'param3D.h'
      include 'common3D.h'


      do   k = 1,km-1
       do  j = 1,jm-1
        do i = 1,im-1

         ar1(i,j,k)   =  aix(i+1,j  ,k  )-aix(i,j,k)
     &                 + ajx(i  ,j+1,k  )-ajx(i,j,k)
     &                 + akx(i  ,j  ,k+1)-akx(i,j,k)
     &                 + aiy(i+1,j  ,k  )-aiy(i,j,k)
     &                 + ajy(i  ,j+1,k  )-ajy(i,j,k)
     &                 + aky(i  ,j  ,k+1)-aky(i,j,k)
     &                 + aiz(i+1,j  ,k  )-aiz(i,j,k)
     &                 + ajz(i  ,j+1,k  )-ajz(i,j,k)
     &                 + akz(i  ,j  ,k+1)-akz(i,j,k)
        enddo
       enddo
      enddo


      write(11,*) 'area check'
      do   k = 1,km-1
       do  j = 1,jm-1
        do i = 1,im-1
         write(11,*) i,j,k,ar1(i,j,k)
        enddo
       enddo
      enddo

      return
      end
      subroutine write_prime
c     ======================

      include 'param3D.h'
      include 'common3D.h'


      write(11,*) 'write prime'
      do   k = 1,km,4
       do  j = 1,jm,4
        do i = 1,im,4
         write(11,*) i,j,k,vx(i,j,k),vy(i,j,k),vz(i,j,k),
     &                      p(i,j,k),u1(i,j,k)
        enddo
       enddo
      enddo

      write(11,*) 'write conservatives'
      do   k = 1,km,4
       do  j = 1,jm,4
        do i = 1,im,4
         write(11,*) i,j,k,u2(i,j,k),u3(i,j,k),u4(i,j,k),
     &                     u5(i,j,k)
        enddo
       enddo
      enddo

      return
      end
      subroutine voltet(x1,x2,x3,x4,x5,x6,x7,x8,
     &                  y1,y2,y3,y4,y5,y6,y7,y8,
     &                  z1,z2,z3,z4,z5,z6,z7,z8,vol)
c     ==============================================

      t1 =  (x4-x5)*(y8*z7-y7*z8)+(y8-y4)*(x8*z7-x7*z8)
     &     +(z4-z5)*(x8*y7-x7*y8)
     &     +x4*(-y5*z7+y7*z5+y5*z8-y8*z5)
     &     +y4*( x5*z7-x7*z5-x5*z8+x8*z5)
     &     +z4*(-x5*y7+x7*y5+x5*y8-x8*y5)

      t2 =  (x5-x7)*(y2*z6-y6*z2)+(y7-y5)*(x2*z6-x6*z2)
     &     +(z5-z7)*(x2*y6-x6*y2)
     &     +x5*(-y7*z6+y6*z7+y7*z2-y2*z7)
     &     +y5*( x7*z6-x6*z7-x7*z2+x2*z7)
     &     +z5*(-x7*y6+x6*y7+x7*y2-x2*y7)

      t3 =  (x4-x2)*(y7*z3-y3*z7)+(y2-y4)*(x7*z3-x3*z7)
     &     +(z4-z2)*(x7*y3-x3*y7)
     &     +x4*(-y2*z3+y3*z2+y2*z7-y7*z2)
     &     +y4*( x2*z3-x3*z2-x2*z7+x7*z2)
     &     +z4*(-x2*y3+x3*y2+x2*y7-x7*y2)

      t4 =  (x4-x1)*(y5*z2-y2*z5)+(y1-y4)*(x5*z2-x2*z5)
     &     +(z4-z1)*(x5*y2-x2*y5)
     &     +x4*(-y1*z2+y2*z1+y1*z5-y5*z1)
     &     +y4*( x1*z2-x2*z1-x1*z5+x5*z1)
     &     +z4*(-x1*y2+x2*y1+x1*y5-x5*y1)

      t5 =  (x5-x4)*(y2*z7-y7*z2)+(y4-y5)*(x2*z7-x7*z2)
     &     +(z5-z4)*(x2*y7-x7*y2)
     &     +x5*(-y4*z7+y7*z4+y4*z2-y2*z4)
     &     +y5*( x4*z7-x7*z4-x4*z2+x2*z4)
     &     +z5*(-x4*y7+x7*y4+x4*y2-x2*y4)

      vol  =  (t1+t2+t3+t4+t5)/6.0

      return
      end
