
c     fortran code of a Demo-2D-Solver (euler method)

c     Author: Bertram Stubert
c     Time of Project-Start: 21.02.2024
    
c     Goal of this code is a demonstration for interested people to get
c     an idea of a 1D-aero-code functionality as simple as possible

c     Application deals with the calculation of a
c     two-block-configuration building a thin curved plate.


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
      open (unit = 15, file = 'Demo2D_psprof.dat')

c     konwn indices
      im = id
      jm = jd

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
      call write_rel_cgns
      call write_boundary
      call write_abs_cgns
      

      write(11,*) ' Demo3D ended'
      write(11,*) ' ************'


      call cpu_time(time)

      write(11,*) ' '
      write(11,'(a13,f10.2)') 'cpu-seconds: ', time
      write(11,'(a13,f10.2)') 'cpu-minutes: ', time/60.0
      write(11,*) ' '

      close (10)
      close (11)
      close (12)
      close (13)
      close (14)
      close (15)

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
        aq       = amax1(100.0,cpzucv*p1(i,j)/u11(i,j))
        ar11(i,j) = sqrt(aq)
        aq       = amax1(100.0,cpzucv*p2(i,j)/u12(i,j))
        ar12(i,j) = sqrt(aq)
       enddo
      enddo


       do  j = 1,jm-1
        do i = 1,im-1

c        block 1

c        volume mean values for speed of sound and velocity

         rav = 0.25*(ar11(i  ,j  )
     &              +ar11(i  ,j+1)
     &              +ar11(i+1,j  )
     &              +ar11(i+1,j+1))

         vxv = 0.25*(vx1 (i  ,j  )
     &              +vx1 (i  ,j+1)
     &              +vx1 (i+1,j  )
     &              +vx1 (i+1,j+1))

         vyv = 0.25*(vy1 (i  ,j  )
     &              +vy1 (i  ,j+1)
     &              +vy1 (i+1,j  )
     &              +vy1 (i+1,j+1))


c     area components at center of single cell volume

         rix = 0.5*(aix1(i+1,j  )+aix1(i,j))
         riy = 0.5*(aiy1(i+1,j  )+aiy1(i,j))

         rjx = 0.5*(ajx1(i  ,j+1)+ajx1(i,j))
         rjy = 0.5*(ajy1(i  ,j+1)+ajy1(i,j))


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

         rla1(i,j) = amax1(rdti,rdtj)


c        block 2

c        volume mean values for speed of sound and velocity

         rav = 0.25*(ar12(i  ,j  )
     &              +ar12(i  ,j+1)
     &              +ar12(i+1,j  )
     &              +ar12(i+1,j+1))

         vxv = 0.25*(vx2 (i  ,j  )
     &              +vx2 (i  ,j+1)
     &              +vx2 (i+1,j  )
     &              +vx2 (i+1,j+1))

         vyv = 0.25*(vy2 (i  ,j  )
     &              +vy2 (i  ,j+1)
     &              +vy2 (i+1,j  )
     &              +vy2 (i+1,j+1))


c     area components at center of single cell volume

         rix = 0.5*(aix2(i+1,j  )+aix2(i,j))
         riy = 0.5*(aiy2(i+1,j  )+aiy2(i,j))

         rjx = 0.5*(ajx2(i  ,j+1)+ajx2(i,j))
         rjy = 0.5*(ajy2(i  ,j+1)+ajy2(i,j))


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

         rla2(i,j) = amax1(rdti,rdtj)

       enddo
      enddo

c     delta(volume)/delta(time) at nodes for multiplying flux-balances

c     sum up for control volumes
      call sfl(0,rla1,step1,rla2,step2)

c     values at boundary elements

c     walls
      do  i = 2,im-1
       step1(i,jm) = 2.0*step1(i,jm)
       step2(i,1)  = 2.0*step2(i,1)
      enddo

      do  i = ile,ite
c      suction side
       step1(i,1 ) = 2.0*step1(i,1)
c      pressure side
       step2(i,jm) = 2.0*step2(i,jm)
      enddo

c     inlet and outlet
      do j = 1,jm
       step1(1 ,j) = 2.0*step1(1 ,j)
       step1(im,j) = 2.0*step1(im,j)
      enddo

      do j = 1,jm-1
       step2(1 ,j) = 2.0*step2(1 ,j)
       step2(im,j) = 2.0*step2(im,j)
      enddo

c     invert to correct step values delta(time)/delta(volume)


      do  j = 1,jm
       do i = 1,im
        step1(i,j) = 1.0/step1(i,j)
       enddo
      enddo

      do  j = 1,jm-1
       do i = 1,im
        step2(i,j) = 1.0/step2(i,j)
       enddo
      enddo

      do i = ile+1,ite-1
        step2(i,jm) = 1.0/step2(i,jm)
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

         flai1(i,j) = (ar11(i,j  )
     &                +ar11(i,j+1))
     &                *aix1(i,j)*0.5
     &               +(ar21(i,j  )
     &                +ar21(i,j+1))
     &                *aiy1(i,j)*0.5


         flai2(i,j) = (ar12(i,j  )
     &                +ar12(i,j+1))
     &                *aix2(i,j)*0.5
     &               +(ar22(i,j  )
     &                +ar22(i,j+1))
     &                *aiy2(i,j)*0.5
       enddo
      enddo

c     calculate flux through area type aj

      do  j = 1,jm
       do i = 1,im-1

         flaj1(i,j) = (ar11(i  ,j)
     &                +ar11(i+1,j))
     &                *ajx1(i,j)*0.5
     &               +(ar21(i  ,j)
     &                +ar21(i+1,j))
     &                *ajy1(i,j)*0.5

         flaj2(i,j) = (ar12(i  ,j)
     &                +ar12(i+1,j))
     &                *ajx2(i,j)*0.5
     &               +(ar22(i  ,j)
     &                +ar22(i+1,j))
     &                *ajy2(i,j)*0.5
       enddo
      enddo


c     boundary condition at walls j=1, j=jm

      if(m.eq.1.or.m.eq.4) then
       do i = ile,ite-2
         flaj1(i,1 ) = 0.0
         flaj2(i,jm) = 0.0
       enddo
       do i = 1,im-1
         flaj1(i,jm) = 0.0
         flaj2(i,1)  = 0.0
       enddo
      endif

      if(m.eq.2) then

c     boundary condition at walls j=1, j=jm

       do i = ile,ite-1
         flaj1(i,1)  = (p1(i  ,1)
     &                 +p1(i+1,1))
     &                 *ajx1(i,1)*0.5

         flaj2(i,jm) = (p2(i  ,jm)
     &                 +p2(i+1,jm))
     &                 *ajx2(i,jm)*0.5
       enddo

       do i = 1,im-1
         flaj1(i,jm) = (p1(i  ,jm)
     &                 +p1(i+1,jm))
     &                 *ajx1(i,jm)*0.5

         flaj2(i,1)  = (p2(i  ,1)
     &                 +p2(i+1,1))
     &                 *ajx2(i,1)*0.5
       enddo
      endif


      if(m.eq.3) then

c     boundary condition at walls j=1, j=jm

       do i = ile,ite-1
         flaj1(i,1)  = (p1(i  ,1)
     &                 +p1(i+1,1))
     &                 *ajy1(i,1)*0.5

         flaj2(i,jm) = (p2(i  ,jm)
     &                 +p2(i+1,jm))
     &                 *ajy2(i,jm)*0.5
       enddo

       do i = 1,im-1
         flaj1(i,jm) = (p1(i  ,jm)
     &                 +p1(i+1,jm))
     &                 *ajy1(i,jm)*0.5

         flaj2(i,1)  = (p2(i  ,1)
     &                 +p2(i+1,1))
     &                 *ajy2(i,1)*0.5
       enddo
      endif


c     fluxbalance for single volume

      do  j = 1,jm-1
       do i = 1,im-1

         fblc1(i,j) = flai1(i,j)-flai1(i+1,j  )
     &               +flaj1(i,j)-flaj1(i  ,j+1)

         fblc2(i,j) = flai2(i,j)-flai2(i+1,j  )
     &               +flaj2(i,j)-flaj2(i  ,j+1)
       enddo
      enddo

c     sum up for flux balances of node control volumes
      call sfl(0,fblc1,fblp1,fblc2,fblp2)


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
        aix1(i,j) = y1(i,j+1) - y1(i,j)
        aix2(i,j) = y2(i,j+1) - y2(i,j)
       enddo
      enddo

      do  j = 1,jm-1
       do i = 1,im
        aiy1(i,j) = -(x1(i,j+1) - x1(i,j))
        aiy2(i,j) = -(x2(i,j+1) - x2(i,j))
       enddo
      enddo


c     area components aj

      do  j = 1,jm
       do i = 1,im-1
        ajx1(i,j) = -(y1(i+1,j) - y1(i,j))
        ajx2(i,j) = -(y2(i+1,j) - y2(i,j))
       enddo
      enddo

      do  j = 1,jm
       do i = 1,im-1
        ajy1(i,j) = x1(i+1,j) - x1(i,j)
        ajy2(i,j) = x2(i+1,j) - x2(i,j)
       enddo
      enddo


c     calculate single volumes

      do  j = 1,jm-1
       do i = 1,im-1

          call gfl(x1 (i  ,j  ),x1 (i+1,j  ),
     &             x1 (i+1,j+1),x1 (i  ,j+1),
     &             y1 (i  ,j  ),y1 (i+1,j  ),
     &             y1 (i+1,j+1),y1 (i  ,j+1),
     &             volc1(i,j))

          call gfl(x2 (i  ,j  ),x2 (i+1,j  ),
     &             x2 (i+1,j+1),x2 (i  ,j+1),
     &             y2 (i  ,j  ),y2 (i+1,j  ),
     &             y2 (i+1,j+1),y2 (i  ,j+1),
     &             volc2(i,j))
       enddo
      enddo

c     sum up for node control volumes
      call sfl(0,volc1,volp1,volc2,volp2)


c      call write_geo



      return
      end
      subroutine gfl(x1,x2,x3,x4,y1,y2,y3,y4,a)
c     =========================================

c     Calculation of cross prduct of two vectors

      a = 0.5*((x3-x1)*(y2-y4)-(x2-x4)*(y3-y1))

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
         x1   (i,j) = 0.0
         y1   (i,j) = 0.0
         u11  (i,j) = 0.0
         u21  (i,j) = 0.0
         u31  (i,j) = 0.0
         u41  (i,j) = 0.0
         vx1  (i,j) = 0.0
         vy1  (i,j) = 0.0
         p1   (i,j) = 0.0
         qc1  (i,j) = 0.0
         u1e1 (i,j) = 0.0
         u2e1 (i,j) = 0.0
         u3e1 (i,j) = 0.0
         u4e1 (i,j) = 0.0
         u1m1 (i,j) = 0.0
         u2m1 (i,j) = 0.0
         u3m1 (i,j) = 0.0
         u4m1 (i,j) = 0.0

         step1(i,j) = 0.0
         flai1(i,j) = 0.0
         flaj1(i,j) = 0.0
         rla1 (i,j) = 0.0
         fblc1(i,j) = 0.0
         fblp1(i,j) = 0.0
         aix1 (i,j) = 0.0
         aiy1 (i,j) = 0.0
         ajx1 (i,j) = 0.0
         ajy1 (i,j) = 0.0
         ar11 (i,j) = 0.0
         ar21 (i,j) = 0.0

         x2   (i,j) = 0.0
         y2   (i,j) = 0.0
         u12  (i,j) = 0.0
         u22  (i,j) = 0.0
         u32  (i,j) = 0.0
         u42  (i,j) = 0.0
         vx2  (i,j) = 0.0
         vy2  (i,j) = 0.0
         p2   (i,j) = 0.0
         qc2  (i,j) = 0.0
         u1e2 (i,j) = 0.0
         u2e2 (i,j) = 0.0
         u3e2 (i,j) = 0.0
         u4e2 (i,j) = 0.0
         u1m2 (i,j) = 0.0
         u2m2 (i,j) = 0.0
         u3m2 (i,j) = 0.0
         u4m2 (i,j) = 0.0

         step2(i,j) = 0.0
         flai2(i,j) = 0.0
         flaj2(i,j) = 0.0
         rla2 (i,j) = 0.0
         fblc2(i,j) = 0.0
         fblp2(i,j) = 0.0
         aix2 (i,j) = 0.0
         aiy2 (i,j) = 0.0
         ajx2 (i,j) = 0.0
         ajy2 (i,j) = 0.0
         ar12 (i,j) = 0.0
         ar22 (i,j) = 0.0

        enddo

        pinol1(j) = 0.0
        pinol2(j) = 0.0

      enddo

      return
      end
      subroutine inlet
c     ================

      include 'param2D.h'
      include 'common2D.h'

      dimension ain(jd)

      arin   = 0.0
      pinsum = 0.0

      do j = 1,jm-1

       ain(j)  = sqrt(aix1(1,j)**2+aiy1(1,j)**2)
       arin    = arin+ain(j)
       pinsum  = pinsum + 0.5*(p1(2,j)+p1(2,j+1))*ain(j)
       ain(j)  = sqrt(aix2(1,j)**2+aiy2(1,j)**2)
       arin    = arin+ain(j)
       pinsum  = pinsum + 0.5*(p2(2,j)+p2(2,j+1))*ain(j)
      enddo

      pinsum = pinsum/arin

c     control output for inlet static pressure
      write( 6,*) 'mean static pressure inlet: ',pinsum
      write(11,*) 'mean static pressure inlet: ',pinsum

c     boundary condition at inflow

      do j = 1,jm

c       from isentropic relations
        p1(1,j)   = (1.0-rfpin)*p1(2,j)+rfpin*pinol1(j)
        pinol1(j) = p1(1,j)
        ts1       = tt1*(p1(1,j)/pt1)**rgzucp
        rma       = sqrt(2.0*cvzurg*(tt1/ts1-1.0))
        u11(1,j)  = p1(1,j)/(rg*ts1)
c       axial inflow
        vx1(1,j)  = rma*sqrt(cpzucv*rg*ts1)
        vy1(1,j)  = 0.0
        u21(1,j)  = u11(1,j)*vx1(1,j)
        u31(1,j)  = u11(1,j)*vy1(1,j)
        u41(1,j)  = p1(1,j)/rgzucv
     &             +0.5*u11(1,j)*(vx1(1,j)**2+vy1(1,j)**2)

        p2(1,j)   = (1.0-rfpin)*p2(2,j)+rfpin*pinol2(j)
        pinol2(j) = p2(1,j)
        ts1       = tt1*(p2(1,j)/pt1)**rgzucp
        rma       = sqrt(2.0*cvzurg*(tt1/ts1-1.0))
        u12(1,j)  = p2(1,j)/(rg*ts1)
c       axial inflow
        vx2(1,j)  = rma*sqrt(cpzucv*rg*ts1)
        vy2(1,j)  = 0.0
        u22(1,j)  = u12(1,j)*vx2(1,j)
        u32(1,j)  = u12(1,j)*vy2(1,j)
        u42(1,j)  = p2(1,j)/rgzucv
     &             +0.5*u12(1,j)*(vx2(1,j)**2+vy2(1,j)**2)

      enddo

      return
      end
      subroutine lax(q1,fbal1,qe1,qm1,q2,fbal2,qe2,qm2)
c     =================================================


      include 'param2D.h'
      include 'common2D.h'

      dimension q1(id,jd),qe1(id,jd),qm1(id,jd)
      dimension fbal1(id,jd)
      dimension q2(id,jd),qe2(id,jd),qm2(id,jd)
      dimension fbal2(id,jd)

c     volume averages single volumes

      do  j = 1,jm-1
       do i = 1,im-1

         qe1(i,j) = 0.25*(q1(i  ,j  )+q1(i  ,j+1)
     &                   +q1(i+1,j+1)+q1(i+1,j  ))

         qe2(i,j) = 0.25*(q2(i  ,j  )+q2(i  ,j+1)
     &                   +q2(i+1,j+1)+q2(i+1,j  ))
       enddo
      enddo


c     sum up volume averages to nodes
      call sfl(1,qe1,qm1,qe2,qm2)


c     update conservative variable

c     regular nodes
      do  j = 1,jm
       do i = 1,im
         q1(i,j)  = (1.0-rfd)*q1(i,j)+rfd*qm1(i,j)
     &                +ft*fbal1(i,j)*step1(i,j)
       enddo
      enddo

      do  j = 1,jm-1
       do i = 1,im
         q2(i,j)  = (1.0-rfd)*q2(i,j)+rfd*qm2(i,j)
     &                +ft*fbal2(i,j)*step2(i,j)
       enddo
      enddo

      do i = ile+1,ite-1
         q2(i,jm)  = (1.0-rfd)*q2(i,jm)+rfd*qm2(i,jm)
     &                +ft*fbal2(i,jm)*step2(i,jm)
      enddo


c     block connectivity
      do i = 1,ile
       q2(i,jm)  = q1(i,1)
      enddo
      do i = ite,im
       q2(i,jm)  = q1(i,1)
      enddo

      return
      end
      subroutine outlet
c     =================

      include 'param2D.h'
      include 'common2D.h'

      dimension aout(jd)

c     boundary condition at outflow

      do j = 2,jm-1

       p1(im,j) = (1.0-rfpout)*p1(im,j)+rfpout*pout
       p2(im,j) = (1.0-rfpout)*p2(im,j)+rfpout*pout

      enddo

      p1(im,1) = (1.0-rfpout)*p1(im,1)+rfpout*pout

c     periodic boundary
      p2(im,1)  = p1(im,jm)
      p2(im,jm) = p1(im,1)


      asout  = 0.0
      pimsum = 0.0

      do j = 1,jm-1

        aout(j)  = sqrt(aix1(im,j)**2+aiy1(im,j)**2)
        asout    = asout+aout(j)

        pimsum = pimsum+0.5*(p1(im,j  )
     &                      +p1(im,j+1))*aout(j)

        aout(j)  = sqrt(aix2(im,j)**2+aiy2(im,j)**2)
        asout    = asout+aout(j)

        pimsum = pimsum+0.5*(p2(im,j  )
     &                      +p2(im,j+1))*aout(j)
      enddo

      pimsum   = pimsum/asout

c     control output outlet static pressure
      write(6,*)  'mean static pressure outlet:',pimsum,' pout:',pout
      write(11,*) 'mean static pressure outlet:',pimsum,' pout:',pout


      return
      end
      subroutine post
c     ===============

      include 'param2D.h'
      include 'common2D.h'

      dimension flowi(id)

c     possible output
c      write(11,*) ' '
c      write(11,*) 'update vx,vy,p,u1'
c      call write_prime

c     calculate massflow for all i-planes
      do  j = 1,jm
       do i = 1,im

          ar11(i,j) = u11(i,j)*vx1(i,j)
          ar21(i,j) = u11(i,j)*vy1(i,j)

          ar12(i,j) = u12(i,j)*vx2(i,j)
          ar22(i,j) = u12(i,j)*vy2(i,j)

       enddo
      enddo

      call flbal(1)

      write(11,*) '   '
      write(11,*) 'massflow ratio at i-planes: '
      write(11,*) '   '

      do   i = 1,im
       flowi(i) = 0.0
       do j = 1,jm-1
        flowi(i) = flowi(i) + flai1(i,j)+flai2(i,j)
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
        ts        = p1(1,j)/(u11(1,j)*rg)
        aq        = cpzucv*rg*ts
        rmaq      = (vx1(1,j)**2+vy1(1,j)**2)/aq
        tt        = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
        pt        = p1(1,j)/((ts/tt)**cpzurg)
        ptinlt    = ptinlt+pt*flai1(1,j)
        ttinlt    = ttinlt+tt*flai1(1,j)

        ts        = p2(1,j)/(u12(1,j)*rg)
        aq        = cpzucv*rg*ts
        rmaq      = (vx2(1,j)**2+vy2(1,j)**2)/aq
        tt        = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
        pt        = p2(1,j)/((ts/tt)**cpzurg)
        ptinlt    = ptinlt+pt*flai2(1,j)
        ttinlt    = ttinlt+tt*flai2(1,j)

c       outlet
        ts        = p1(im,j)/(u11(im,j)*rg)
        aq        = cpzucv*rg*ts
        rmaq      = (vx1(im,j)**2+vy1(im,j)**2)/aq
        tt        = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
        pt        = p1(im,j)/((ts/tt)**cpzurg)
        ptout     = ptout+pt*flai1(im,j)
        ttout     = ttout+tt*flai1(im,j)

        ts        = p2(im,j)/(u12(im,j)*rg)
        aq        = cpzucv*rg*ts
        rmaq      = (vx2(im,j)**2+vy2(im,j)**2)/aq
        tt        = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
        pt        = p2(im,j)/((ts/tt)**cpzurg)
        ptout     = ptout+pt*flai2(im,j)
        ttout     = ttout+tt*flai2(im,j)
      enddo

      ptinlt = ptinlt/flowin
      ttinlt = ttinlt/flowin

      ptout  = ptout/flowot
      ttout  = ttout/flowot

      ptv    = ptout/ptinlt
      ttv    = ttout/ttinlt

c     lift and drag

      fdrag = 0.0
      flift = 0.0

      do i = ile,ite-1
       fdrag = fdrag + 0.5*(p1(i,1) +p1(i+1,1)) *ajx1(i,1)*0.001
       fdrag = fdrag + 0.5*(p2(i,jm)+p2(i+1,jm))*ajx2(i,jm)*0.001
       flift = flift + 0.5*(p1(i,1) +p1(i+1,1)) *ajy1(i,1)*0.001
       flift = flift + 0.5*(p2(i,jm)+p2(i+1,jm))*ajy2(i,jm)*0.001
      enddo

      xmom1 = 0.0
      xmom2 = 0.0
      do j = 1,jm-1
       xmom1 = xmom1+0.5*(u11(1,j+1)+u11(1,j))
     &             *(0.5*(vx1(1,j+1)+vx1(1,j)))**2*aix1(1,j)
       xmom1 = xmom1+0.5*(u12(1,j+1)+u12(1,j))
     &             *(0.5*(vx2(1,j+1)+vx2(1,j)))**2*aix2(1,j)

       xmom2 = xmom2+0.5*(u11(im,j+1)+u11(im,j))
     &             *(0.5*(vx1(im,j+1)+vx1(im,j)))**2*aix1(im,j)
       xmom2 = xmom2+0.5*(u12(im,j+1)+u12(im,j))
     &             *(0.5*(vx2(im,j+1)+vx2(im,j)))**2*aix2(im,j)
      enddo

      write(11,*) '   '
      write(11,*) 'mass averaged total pressure and temperature:'
      write(11,*) '   '

      write(11,'(A17,F12.5)') 'pt inlet   = ',ptinlt
      write(11,'(A17,F12.5)') 'tt inlet   = ',ttinlt
      write(11,'(A17,F12.5)') 'pt outlet  = ',ptout
      write(11,'(A17,F12.5)') 'tt outlet  = ',ttout
      write(11,'(A17,F12.5)') 'ptout/ptin = ',ptv
      write(11,'(A17,F12.5)') 'ttout/ttin = ',ttv
      write(11,'(A17,F12.5)') 'lift       = ',flift
      write(11,'(A17,F12.5)') 'drag       = ',fdrag
      write(11,'(A17,F12.5)') 'xmom_in    = ',xmom1
      write(11,'(A17,F12.5)') 'xmom_out   = ',xmom2

      write( 6,'(A17,F12.5)') 'massflow inlet  = ',flowin
      write( 6,'(A17,F12.5)') 'massflow outlet = ',flowot
      write( 6,'(A17,F12.5)') 'pt inlet   = ',ptinlt
      write( 6,'(A17,F12.5)') 'tt inlet   = ',ttinlt
      write( 6,'(A17,F12.5)') 'pt outlet  = ',ptout
      write( 6,'(A17,F12.5)') 'tt outlet  = ',ttout
      write( 6,'(A17,F12.5)') 'ptout/ptin = ',ptv
      write( 6,'(A17,F12.5)') 'ttout/ttin = ',ttv
      write( 6,'(A17,F12.5)') 'xmom_in    = ',xmom1
      write( 6,'(A17,F12.5)') 'xmom_out   = ',xmom2

      write(11,*) '   '
      write(11,*) '   '
      write(11,*) '   '

c     blade surface pressure distribution
      do i = ile,ite
       write(15,*) x1(i,1),p1(i,1)
      enddo

      do i = ite,ile,-1
       write(15,*) x2(i,jm),p2(i,jm)
      enddo



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
      write(11,'(A9,F10.3)') 'ptvbl  = ', ptvbl
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
      subroutine read_coord
c     =====================

      include 'cgnslib_f.h'
      include 'param2D.h'
      include 'common2D.h'

      write(11,*) 'read_coord','  im=',im,'jm=',jm
      do  j = 1,jm
       do i = 1,im
        read(14,*) x1(i,j),y1(i,j)
       enddo
      enddo
      do  j = 1,jm
       do i = 1,im
        read(14,*) x2(i,j),y2(i,j)
       enddo
      enddo

c     possible simple test
c     shoebox
c      do  j = 1,jm
c       do i = 2,im
c        y1(i,j) = y1(1,j)
c        y2(i,j) = y2(1,j)
c       enddo
c      enddo

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

          u11(i,j) = u1star
          p1 (i,j) = pstart
          vx1(i,j) = vxstar
          vy1(i,j) = vystar

          u12(i,j) = u1star
          p2 (i,j) = pstart
          vx2(i,j) = vxstar
          vy2(i,j) = vystar

        enddo
       enddo

      do  j = 1,jm
       do i = 1,im

         u21(i,j) = u11(i,j)*vx1(i,j)
         u31(i,j) = u11(i,j)*vy1(i,j)
         u41(i,j) = p1(i,j)/rgzucv
     &              +0.5*u11(i,j)*(vx1(i,j)**2+vy1(i,j)**2)

         u22(i,j) = u12(i,j)*vx2(i,j)
         u32(i,j) = u12(i,j)*vy2(i,j)
         u42(i,j) = p2(i,j)/rgzucv
     &             +0.5*u12(i,j)*(vx2(i,j)**2+vy2(i,j)**2)
       enddo

        pinol1(j) = p1(1,j)
        pinol2(j) = p2(1,j)

      enddo


      return
      end
      subroutine sfl(m,dudt1,dqdt1,dudt2,dqdt2)
c     =========================================

      include 'param2D.h'
      include 'common2D.h'

c     sum up single balances to control volume

      dimension dudt1(id,jd),dqdt1(id,jd)
      dimension dudt2(id,jd),dqdt2(id,jd)

      if(m.eq.0) then
       f4 = 1.0
       f2 = 1.0
      else
       f4 = 0.25
       f2 = 0.50
      endif


c     regular field
       do  j = 2,jm-1
        dO i = 2,im-1

         dqdt1(i,j) = (dudt1(i-1,j-1)
     &                +dudt1(i-1,j  )
     &                +dudt1(i  ,j  )
     &                +dudt1(i  ,j-1))*f4
         dqdt2(i,j) = (dudt2(i-1,j-1)
     &                +dudt2(i-1,j  )
     &                +dudt2(i  ,j  )
     &                +dudt2(i  ,j-1))*f4
       enddo
      enddo

c     values at boundaries

c     walls
      do  i = 2,im-1

        dqdt1(i,jm) = (dudt1(i-1,jm-1)
     &                +dudt1(i  ,jm-1))*f2
        dqdt2(i,1)  = (dudt2(i-1,1)
     &                +dudt2(i  ,1))*f2
      enddo

c     block connectivity
      do  i = 2,ile

        dqdt1(i,1) = (dudt1(i-1,1)
     &               +dudt1(i  ,1)
     &               +dudt2(i-1,jm-1)
     &               +dudt2(i  ,jm-1))*f4
      enddo

      do  i = ite,im-1

        dqdt1(i,1) = (dudt1(i-1,1)
     &               +dudt1(i  ,1)
     &               +dudt2(i-1,jm-1)
     &               +dudt2(i  ,jm-1))*f4
      enddo

      do  i = ile+1,ite-1
c       suction side
        dqdt1(i,1)  = (dudt1(i-1,1)
     &                +dudt1(i  ,1))*f2

c       pressure side
        dqdt2(i,jm) = (dudt2(i-1,jm-1)
     &                +dudt2(i  ,jm-1))*f2
      enddo

      do  j = 2,jm-1
c       inlet
        dqdt1(1,j)  = (dudt1(1,j-1)
     &                +dudt1(1,j  ))*f2
        dqdt2(1,j)  = (dudt2(1,j-1)
     &                +dudt2(1,j  ))*f2

c       outlet
        dqdt1(im,j) = (dudt1(im-1,j-1)
     &                +dudt1(im-1,j  ))*f2
        dqdt2(im,j) = (dudt2(im-1,j-1)
     &                +dudt2(im-1,j  ))*f2
      enddo

c     inlet block connectivity and walls
      dqdt1(1,1)  = (dudt1(1,1)
     &              +dudt2(1,jm-1))*f2
      dqdt1(1,jm) =  dudt1(1,jm-1)
      dqdt2(1,1)  =  dudt2(1,1)


c     outlet block connectivity and walls
      dqdt1(im,1)  = (dudt1(im-1,1)
     &               +dudt2(im-1,jm-1))*f2
      dqdt1(im,jm) =  dudt1(im-1,jm-1)
      dqdt2(im,1)  =  dudt2(im-1,1)


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

          ar11(i,j) = u11(i,j)*vx1(i,j)
          ar21(i,j) = u11(i,j)*vy1(i,j)

          ar12(i,j) = u12(i,j)*vx2(i,j)
          ar22(i,j) = u12(i,j)*vy2(i,j)

        enddo
       enddo

c      calculate flux balances
       call flbal(1)

c      update u1
       call lax(u11,fblp1,u1e1,u1m1,u12,fblp2,u1e2,u1m2)


c      calculate massflow at inlet and outlet
       do   i = 1,im,im-1
        flowi(i) = 0.0
        do j = 1,jm-1
          flowi(i) = flowi(i) + flai1(i,j)+flai2(i,j)
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

          ar11(i,j) = u21(i,j)*vx1(i,j)+p1(i,j)
          ar21(i,j) = u21(i,j)*vy1(i,j)

          ar12(i,j) = u22(i,j)*vx2(i,j)+p2(i,j)
          ar22(i,j) = u22(i,j)*vy2(i,j)

        enddo
       enddo

c      calculate flux balances
       call flbal(2)

c      update u2
       call lax(u21,fblp1,u2e1,u2m1,u22,fblp2,u2e2,u2m2)


c      update y-momentum
c      *****************

c      calculate flux variables

       do  j = 1,jm
        do i = 1,im

          ar11(i,j) = u31(i,j)*vx1(i,j)
          ar21(i,j) = u31(i,j)*vy1(i,j)+p1(i,j)

          ar12(i,j) = u32(i,j)*vx2(i,j)
          ar22(i,j) = u32(i,j)*vy2(i,j)+p2(i,j)

        enddo
       enddo

c      calculate flux balances
       call flbal(3)

c      update u3
       call lax(u31,fblp1,u3e1,u3m1,u32,fblp2,u3e2,u3m2)


c      update total energy
c      *******************

c      calculate flux variables

       do  j = 1,jm
        do i = 1,im

          ar11(i,j) = (u41(i,j)+p1(i,j))*vx1(i,j)
          ar21(i,j) = (u41(i,j)+p1(i,j))*vy1(i,j)

          ar12(i,j) = (u42(i,j)+p2(i,j))*vx2(i,j)
          ar22(i,j) = (u42(i,j)+p2(i,j))*vy2(i,j)

        enddo
       enddo


c      calculate flux balances
       call flbal(4)

c      update u4
       call lax(u41,fblp1,u4e1,u4m1,u42,fblp2,u4e2,u4m2)


c      new primary variables

       do  j = 1,jm
        do i = 1,im

          qc1(i,j) = amax1(0.00001,abs(vx1(i,j)))
          qc2(i,j) = amax1(0.00001,abs(vx2(i,j)))

          vx1(i,j) = u21(i,j)/u11(i,j)
          vy1(i,j) = u31(i,j)/u11(i,j)

          vx2(i,j) = u22(i,j)/u12(i,j)
          vy2(i,j) = u32(i,j)/u12(i,j)

          p1(i,j) = (u41(i,j) - 0.5*(vx1(i,j)**2+vy1(i,j)**2)
     &              *u11(i,j))*rgzucv
          p2(i,j) = (u42(i,j) - 0.5*(vx2(i,j)**2+vy2(i,j)**2)
     &              *u12(i,j))*rgzucv
        enddo

         pinol1(j) = p1 (1,j)
         pinol2(j) = p2 (1,j)

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
          erij = 1.0-sqrt(vx1(i,j)**2+vy1(i,j)**2)/qc1(i,j)
          if(abs(erij).gt.1000.0) then
           write(6 ,*) 'vx:',vx1(i,j),'vy:',vy1(i,j),'qc:',qc1(i,j)
           write(6 ,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(11,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(6 ,*) 'erij:',erij,' Demo3D diverges! in block 1'
           write(11,*) 'erij:',erij,' Demo3D diverges! in block 1'
           goto 999
          endif

          if(erij.gt.evmax) then
           iwbl  = 1
           evmax = erij
           iemax = i
           jemax = j
          endif
        enddo
       enddo

       do  j = 1,jm
        do i = 1,im
          erij = 1.0-sqrt(vx2(i,j)**2+vy2(i,j)**2)/qc2(i,j)
          if(abs(erij).gt.1000.0) then
           write(6 ,*) 'vx:',vx2(i,j),'vy:',vy2(i,j),'qc:',qc2(i,j)
           write(6 ,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(11,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(6 ,*) 'erij:',erij,' Demo3D diverges! in block 2'
           write(11,*) 'erij:',erij,' Demo3D diverges! in block 2'
           goto 999
          endif

          if(erij.gt.evmax) then
           iwbl  = 2
           evmax = erij
           iemax = i
           jemax = j
          endif
        enddo
       enddo

       write(6 , *) 'n = ',n,' at',iemax,jemax,iwbl,'evmax = ',evmax
       write(11, *) 'n = ',n,' at',iemax,jemax,iwbl,'evmax = ',evmax
       write(13, *) evmax

c      check convergence
       if(abs(flrate-1.0).lt.ermax.and.n.gt.1000) then
        nconv = nconv +1
        write(6 ,*) 'nconv = ',nconv
        write(11,*) 'nconv = ',nconv
        if(nconv.eq.10) goto 999
       endif

      enddo

  999 continue


      call write_prime

      return
      end
      subroutine write_boundary
c     =========================

      include 'cgnslib_f.h'
      include 'param2D.h'
      include 'common2D.h'


      character*32 basename,zonename,famname
      character*32 donorname,boconame,connectname

      dimension isize(3,3),ipnts(3,3),itransfrm(3),ipntsdonor(3,3)

      call cg_open_f
     &     ('Demo2D_result_rel.cgns',CG_MODE_MODIFY,index_file,ier)
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

      index_zone = 2
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

      write(6,*) 'boundaries Block 1'
c     *******************************


      write(6,*) 'boundary Inlet'
c     **************************
      famname = 'Inlet'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family Inlet'
      write(6,*) 'Inlet: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

c     index space for Inlet
      ipnts(1,1) = 1
      ipnts(2,1) = 1
      ipnts(3,1) = 1
      ipnts(1,2) = 1
      ipnts(2,2) = jd
      ipnts(3,2) = kd

      index_zone = 1
      boconame   = 'Block1_ilo'

      call cg_boco_write_f(index_file,index_base,index_zone,
     &   boconame,BCTunnelInflow,PointRange,2,ipnts,index_bc,ier)

      index_zone = 2
      boconame   = 'Block2_ilo'

      call cg_boco_write_f(index_file,index_base,index_zone,
     &   boconame,BCTunnelInflow,PointRange,2,ipnts,index_bc,ier)


      write(6,*) 'boundary Outlet'
c     ****************************

      famname = 'Outlet'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family Outlet'
      write(6,*) 'Outlet: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

      ipnts(1,1) = id
      ipnts(2,1) = 1
      ipnts(3,1) = 1
      ipnts(1,2) = id
      ipnts(2,2) = jd
      ipnts(3,2) = kd

      index_zone = 1
      boconame   = 'Block1_ihi'

      call cg_boco_write_f(index_file,index_base,index_zone,
     &   boconame,BCTunnelOutflow,PointRange,2,ipnts,index_bc,ier)

      index_zone = 2
      boconame   = 'Block2_ihi'

      call cg_boco_write_f(index_file,index_base,index_zone,
     &   boconame,BCTunnelOutflow,PointRange,2,ipnts,index_bc,ier)


      famname = 'WallLeft'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family WallLeft'
      write(6,*) 'WallLeft: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

      write(6,*) 'boundary j = 1'
c     ***************************

      ipnts(1,1) = 1
      ipnts(2,1) = 1
      ipnts(3,1) = 1
      ipnts(1,2) = id
      ipnts(2,2) = 1
      ipnts(3,2) = kd

      index_zone = 1
      boconame   = 'Block1_jlo'

      call cg_boco_write_f(index_file,index_base,index_zone,
     &   boconame,BCSymmetryPlane,PointRange,2,ipnts,index_bc,ier)

      index_zone = 2
      boconame   = 'Block2_jlo'

      call cg_boco_write_f(index_file,index_base,index_zone,
     &   boconame,BCSymmetryPlane,PointRange,2,ipnts,index_bc,ier)


      famname = 'WallRight'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family WallRight'
      write(6,*) 'WallLeft: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

      write(6,*) 'boundary j = jd'
c     *****************************

      ipnts(1,1) = 1
      ipnts(2,1) = jd
      ipnts(3,1) = 1
      ipnts(1,2) = id
      ipnts(2,2) = jd
      ipnts(3,2) = kd

      index_zone = 1
      boconame   = 'Block1_jhi'

      call cg_boco_write_f(index_file,index_base,index_zone,
     &   boconame,BCSymmetryPlane,PointRange,2,ipnts,index_bc,ier)

      index_zone = 2
      boconame   = 'Block2_jhi'

      call cg_boco_write_f(index_file,index_base,index_zone,
     &   boconame,BCSymmetryPlane,PointRange,2,ipnts,index_bc,ier)


      famname = 'SuctionSide'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family WallRight'
      write(6,*) 'WallLeft: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

      write(6,*) 'boundary k = 1 suction side'
c     ****************************************

      ipnts(1,1) = ile
      ipnts(2,1) = 1
      ipnts(3,1) = 1
      ipnts(1,2) = ite
      ipnts(2,2) = jd
      ipnts(3,2) = 1

      index_zone = 1
      boconame   = 'Block1 suction side'

      call cg_boco_write_f(index_file,index_base,index_zone,
     &   boconame,BCWallViscous,PointRange,2,ipnts,index_bc,ier)


      famname = 'PressureSide'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family WallRight'
      write(6,*) 'WallLeft: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

      ipnts(1,1) = ile
      ipnts(2,1) = 1
      ipnts(3,1) = kd
      ipnts(1,2) = ite
      ipnts(2,2) = jd
      ipnts(3,2) = kd

      index_zone = 2
      boconame   = 'Block2 pressure side'

      call cg_boco_write_f(index_file,index_base,index_zone,
     &   boconame,BCWallViscous,PointRange,2,ipnts,index_bc,ier)



c close CGNS file
      call cg_close_f(index_file,ier)
      write(6,'('' boundary panels written to Demo2D_result.cgns'')')

      return
      end
      subroutine write_rel_cgns
c     =========================

      include 'cgnslib_f.h'
      include 'param2D.h'
      include 'common2D.h'

      dimension x3d1 (id,jd,kd),y3d1  (id,jd,kd),z3d1 (id,jd,kd)
      dimension vx3d1(id,jd,kd),vy3d1 (id,jd,kd),vz3d1(id,jd,kd)
      dimension p3d1 (id,jd,kd),u13d1 (id,jd,kd)
      dimension x3d2 (id,jd,kd),y3d2  (id,jd,kd),z3d2 (id,jd,kd)
      dimension vx3d2(id,jd,kd),vy3d2 (id,jd,kd),vz3d2(id,jd,kd)
      dimension p3d2 (id,jd,kd),u13d2 (id,jd,kd)

      dimension isize(3,3)

      character basename*32,zonename*32,solname*32

      km = 3

      do   k = 1,km
       do  j = 1,jm
        do i = 1,im

         x3d1 (i,j,k) = x1(i,j)
         y3d1 (i,j,k) = y1(i,j)
         z3d1 (i,j,k) = real(k)*0.001
c        relatvsystem
         vx3d1(i,j,k) = vx1(i,j)

         vy3d1(i,j,k) = vy1(i,j)
         vz3d1(i,j,k) = 0.0
         p3d1 (i,j,k) = p1(i,j)
         u13d1(i,j,k) = u11(i,j)

         x3d2 (i,j,k) = x2(i,j)
         y3d2 (i,j,k) = y2(i,j)
         z3d2 (i,j,k) = real(k)*0.001
c        relatvsystem
         vx3d2(i,j,k) = vx2(i,j)

         vy3d2(i,j,k) = vy2(i,j)
         vz3d2(i,j,k) = 0.0
         p3d2 (i,j,k) = p2(i,j)
         u13d2(i,j,k) = u12(i,j)

        enddo
       enddo
      enddo



c     WRITE X, Y, Z GRID POINTS TO CGNS FILE
c     open CGNS file for write

      call cg_open_f
     &     ('Demo2D_result_rel.cgns',CG_MODE_WRITE,index_file,ier)

      write(6,*) 'index_file =',index_file

      if (ier .ne. CG_OK) call cg_error_exit_f

c create base (user can give any name)

      basename = 'Base 1'
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


c      write Block 1
c      *************

c create zone
        call cg_zone_write_f(index_file,index_base,zonename,isize,
     &                       Structured,index_zone,ier)

        write(6,*) 'index_zone =',index_zone

c write grid coordinates (user must use SIDS-standard names here)
        call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateX',x3d1,index_coord,ier)
        write(6,*) 'index_coord X =',index_coord

        call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateY',y3d1,index_coord,ier)
        write(6,*) 'index_coord Y =',index_coord

        call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateZ',z3d1,index_coord,ier)
        write(6,*) 'index_coord Z =',index_coord

c    define flow solution node name (user can give any name)
      solname = 'FlowSolution'

c    create flow solution node
      call cg_sol_write_f(index_file,index_base,index_zone,solname,
     &                    Vertex,index_flow,ier)

c    write flow solution (user must use SIDS-standard names here)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Density',u13d1,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Pressure',p3d1,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityX',vx3d1,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityY',vy3d1,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityZ',vz3d1,index_field,ier)


c      write Block 2
c      *************

       zonename = 'Block 2'


c create zone
        call cg_zone_write_f(index_file,index_base,zonename,isize,
     &                       Structured,index_zone,ier)

        write(6,*) 'index_zone =',index_zone

c write grid coordinates (user must use SIDS-standard names here)
        call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateX',x3d2,index_coord,ier)
        write(6,*) 'index_coord X =',index_coord

        call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateY',y3d2,index_coord,ier)
        write(6,*) 'index_coord Y =',index_coord

        call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateZ',z3d2,index_coord,ier)
        write(6,*) 'index_coord Z =',index_coord


c    define flow solution node name (user can give any name)
      solname = 'FlowSolution'

c    create flow solution node
      call cg_sol_write_f(index_file,index_base,index_zone,solname,
     &                    Vertex,index_flow,ier)

c    write flow solution (user must use SIDS-standard names here)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Density',u13d2,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Pressure',p3d2,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityX',vx3d2,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityY',vy3d2,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityZ',vz3d2,index_field,ier)



c close CGNS file
      call cg_close_f(index_file,ier)
      write(6,'('' Results wrote to file Demo2D_result.cgns'')')



      return
      end
      subroutine write_abs_cgns
c     =========================

      include 'cgnslib_f.h'
      include 'param2D.h'
      include 'common2D.h'

      dimension x3d1 (id,jd,kd),y3d1  (id,jd,kd),z3d1 (id,jd,kd)
      dimension vx3d1(id,jd,kd),vy3d1 (id,jd,kd),vz3d1(id,jd,kd)
      dimension p3d1 (id,jd,kd),u13d1 (id,jd,kd)
      dimension x3d2 (id,jd,kd),y3d2  (id,jd,kd),z3d2 (id,jd,kd)
      dimension vx3d2(id,jd,kd),vy3d2 (id,jd,kd),vz3d2(id,jd,kd)
      dimension p3d2 (id,jd,kd),u13d2 (id,jd,kd)

      dimension isize(3,3)

      character basename*32,zonename*32,solname*32

      km = 3

      do   k = 1,km
       do  j = 1,jm
        do i = 1,im

         x3d1 (i,j,k) = x1(i,j)
         y3d1 (i,j,k) = y1(i,j)
         z3d1 (i,j,k) = real(k)*0.001

c        absolutsystem
         vx3d1(i,j,k) = vx1(i,j)-vx1(1,1)

         vy3d1(i,j,k) = vy1(i,j)
         vz3d1(i,j,k) = 0.0

         ts           = p1(i,j)/(rg*u11(i,j))
         pt           = p1(i,j)*(tt1/ts)**cpzurg
         p3d1 (i,j,k) = pt
         u13d1(i,j,k) = u11(i,j)

         x3d2 (i,j,k) = x2(i,j)
         y3d2 (i,j,k) = y2(i,j)
         z3d2 (i,j,k) = real(k)*0.001

c        absolutsystem
         vx3d2(i,j,k) = vx2(i,j)-vx1(1,1)

         vy3d2(i,j,k) = vy2(i,j)
         vz3d2(i,j,k) = 0.0
         ts           = p2(i,j)/(rg*u12(i,j))
         pt           = p2(i,j)*(tt1/ts)**cpzurg
         p3d2 (i,j,k) = pt
         u13d2(i,j,k) = u12(i,j)

        enddo
       enddo
      enddo



c     WRITE X, Y, Z GRID POINTS TO CGNS FILE
c     open CGNS file for write

      call cg_open_f
     &     ('Demo2D_result_abs.cgns',CG_MODE_WRITE,index_file,ier)

      write(6,*) 'index_file =',index_file

      if (ier .ne. CG_OK) call cg_error_exit_f

c create base (user can give any name)

      basename = 'Base 1'
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


c      write Block 1
c      *************

c create zone
        call cg_zone_write_f(index_file,index_base,zonename,isize,
     &                       Structured,index_zone,ier)

        write(6,*) 'index_zone =',index_zone

c write grid coordinates (user must use SIDS-standard names here)
        call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateX',x3d1,index_coord,ier)
        write(6,*) 'index_coord X =',index_coord

        call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateY',y3d1,index_coord,ier)
        write(6,*) 'index_coord Y =',index_coord

        call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateZ',z3d1,index_coord,ier)
        write(6,*) 'index_coord Z =',index_coord

c    define flow solution node name (user can give any name)
      solname = 'FlowSolution'

c    create flow solution node
      call cg_sol_write_f(index_file,index_base,index_zone,solname,
     &                    Vertex,index_flow,ier)

c    write flow solution (user must use SIDS-standard names here)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Density',u13d1,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Pressure',p3d1,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityX',vx3d1,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityY',vy3d1,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityZ',vz3d1,index_field,ier)


c      write Block 2
c      *************

       zonename = 'Block 2'


c create zone
        call cg_zone_write_f(index_file,index_base,zonename,isize,
     &                       Structured,index_zone,ier)

        write(6,*) 'index_zone =',index_zone

c write grid coordinates (user must use SIDS-standard names here)
        call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateX',x3d2,index_coord,ier)
        write(6,*) 'index_coord X =',index_coord

        call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateY',y3d2,index_coord,ier)
        write(6,*) 'index_coord Y =',index_coord

        call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateZ',z3d2,index_coord,ier)
        write(6,*) 'index_coord Z =',index_coord


c    define flow solution node name (user can give any name)
      solname = 'FlowSolution'

c    create flow solution node
      call cg_sol_write_f(index_file,index_base,index_zone,solname,
     &                    Vertex,index_flow,ier)

c    write flow solution (user must use SIDS-standard names here)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Density',u13d2,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Pressure',p3d2,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityX',vx3d2,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityY',vy3d2,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityZ',vz3d2,index_field,ier)



c close CGNS file
      call cg_close_f(index_file,ier)
      write(6,'('' Results wrote to file Demo2D_result_abs.cgns'')')



      return
      end
      subroutine write_geo
c     ====================

      include 'param2D.h'
      include 'common2D.h'


       do  j = 1,jm-1
        do i = 1,im-1

         ar11(i,j)  =  aix1(i+1,j  )-aix1(i,j)
     &               + ajx1(i  ,j+1)-ajx1(i,j)
     &               + aiy1(i+1,j  )-aiy1(i,j)
     &               + ajy1(i  ,j+1)-ajy1(i,j)
         ar12(i,j)  =  aix2(i+1,j  )-aix2(i,j)
     &               + ajx2(i  ,j+1)-ajx2(i,j)
     &               + aiy2(i+1,j  )-aiy2(i,j)
     &               + ajy2(i  ,j+1)-ajy2(i,j)
       enddo
      enddo

      write(11,*) 'area aix,aiy'
      do  i = 1,im
       do j = 1,jm-1
         write(11,*) i,j,aix1(i,j),aiy1(i,j),aix2(i,j),aiy2(i,j)
       enddo
      enddo

      write(11,*) 'area ajx,ajy'
      do  i = 1,im-1
       do j = 1,jm
         write(11,*) i,j,ajx1(i,j),ajy1(i,j),ajx2(i,j),ajy2(i,j)
       enddo
      enddo

      write(11,*) 'volc'
      do  i = 1,im-1
       do j = 1,jm-1
         write(11,*) i,j,volc1(i,j),volc2(i,j)
       enddo
      enddo

      write(11,*) 'area check'

      do  j = 1,jm-1
       do i = 1,im-1
         write(11,*) i,j,ar11(i,j),ar12(i,j)
       enddo
      enddo

      return
      end
      subroutine write_prime
c     ======================

      include 'param2D.h'
      include 'common2D.h'


      write(11,*) 'write prime block 1',' im=',im,'jm=',jm
          j = 1
       do i = 1,im,2
         write(11,*) i,j,vx1(i,j),vy1(i,j),p1(i,j),u11(i,j)
       enddo
          j = (jm-1)/2+1
       do i = 1,im,2
         write(11,*) i,j,vx1(i,j),vy1(i,j),p1(i,j),u11(i,j)
       enddo
          j = jm
       do i = 1,im,2
         write(11,*) i,j,vx1(i,j),vy1(i,j),p1(i,j),u11(i,j)
       enddo


      write(11,*) 'write prime block 2',' im=',im,'jm=',jm
          j = 1
       do i = 1,im,2
         write(11,*) i,j,vx2(i,j),vy2(i,j),p2(i,j),u12(i,j)
       enddo
          j = (jm-1)/2+1
       do i = 1,im,2
         write(11,*) i,j,vx2(i,j),vy2(i,j),p2(i,j),u12(i,j)
       enddo
          j = jm
       do i = 1,im,2
         write(11,*) i,j,vx2(i,j),vy2(i,j),p2(i,j),u12(i,j)
       enddo




      return
      end
