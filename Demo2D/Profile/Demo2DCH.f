
c     fortran code of a Demo-2D-Solver (euler method) C-type-grid

c     Author: Bertram Stubert
c     Time of Project-Start: 26.11.2023
    
c     Goal of this code is a demonstration for interested people to get
c     an idea of a 2D-aero-code functionality as possible simple

c     Application deals with the calculation of a
c     two-block-configuration building a wing profile
c     within a wind tunnel with slip walls.
      
      
      include 'param2D.h'
      include 'common2D.h'

      open (unit = 10, file = 'Demo2D_input.dat')
      open (unit = 11, file = 'Demo2D_output.dat')
      open (unit = 12, file = 'Demo2D_konti.dat')
      open (unit = 13, file = 'Demo2D_evmax.dat')
      open (unit = 14, file = 'prof2_grid.dat')
      open (unit = 15, file = 'Demo2D_psprof.dat')



c     known indices
      im  = id
      jm  = jd
      imh = idh
      jmh = jdh

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

c      call write_test

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
        aq       = amax1(100.0,cpzucv*p(i,j)/u1(i,j))
        ar1(i,j) = sqrt(aq)
       enddo
      enddo

      do  j = 1,jmh
       do i = 1,imh
        aq        = amax1(100.0,cpzucv*ph(i,j)/u1h(i,j))
        ar1h(i,j) = sqrt(aq)
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

       do  j = 1,jmh-1
        do i = 1,imh-1

c        volume mean values for speed of sound and velocity

         rav = 0.25*(ar1h(i  ,j  )
     &              +ar1h(i  ,j+1)
     &              +ar1h(i+1,j  )
     &              +ar1h(i+1,j+1))

         vxv = 0.25*(vxh (i  ,j  )
     &              +vxh (i  ,j+1)
     &              +vxh (i+1,j  )
     &              +vxh (i+1,j+1))

         vyv = 0.25*(vyh (i  ,j  )
     &              +vyh (i  ,j+1)
     &              +vyh (i+1,j  )
     &              +vyh (i+1,j+1))


c     area components at center of single cell volume

         rix = 0.5*(aixh(i+1,j  )+aixh(i,j))
         riy = 0.5*(aiyh(i+1,j  )+aiyh(i,j))

         rjx = 0.5*(ajxh(i  ,j+1)+ajxh(i,j))
         rjy = 0.5*(ajyh(i  ,j+1)+ajyh(i,j))


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

         rlah(i,j) = amax1(rdti,rdtj)

       enddo
      enddo

c     delta(volume)/delta(time) at nodes for multiplying flux-balances

c     sum up for control volumes
      call sfl(rla,step,rlah,steph,0)

c     values at boundary elements

c     outer walls and inlet
      do  i = 1,ile1
       step(i,jm)  = 2.0*step(i,jm)
      enddo
      do  i = ile2,im
       step(i,jm)  = 2.0*step(i,jm)
      enddo

c     profil surface
      do  i = ite1,ite2
       step(i,1)   = 2.0*step(i,1)
      enddo

c     outlet
      do j = 2,jm-1
       step(1 ,j) = 2.0*step(1 ,j)
       step(im,j) = 2.0*step(im,j)
      enddo

      step(1,1) = 2.0*step(1,1)

      do i = 1,imh-1
        steph(i,1)   = 2.0*steph(i,1)
        steph(i,jmh) = 2.0*steph(i,jmh)
      enddo

c     invert to correct step values delta(time)/delta(volume)

c     regular field
      do  j = 2,jm
       do i = 1,im
        step(i,j) = 1.0/step(i,j)
       enddo
      enddo

c     block connectivity and blade wall
      do i = 1,ite2-1
        step(i,1) = 1.0/step(i,1)
      enddo

c     h-block
      do  j = 1,jmh
       do i = 1,imh-1
        steph(i,j) = 1.0/steph(i,j)
       enddo
      enddo


c      write(11,*) 'step'
c      call write_q(step)
c      write(11,*) 'rla'
c      call write_q(rla)
c      write(11,*) 'steph'
c      call write_qh(steph)
c      write(11,*) 'rlah'
c      call write_qh(rlah)

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
       do i = 1,ile1-1
         flaj(i,jm) = 0.0
       enddo
       do i = ile2,im-1
         flaj(i,jm) = 0.0
       enddo
       do i = ite1,ite2-1
         flaj(i,1) = 0.0
       enddo
      endif

      if(m.eq.2) then

c     boundary condition at walls j=1, j=jm

       do i = 1,ile1-1
         flaj(i,jm) = (p(i,jm)+p(i+1,jm))*ajx(i,jm)*0.5
       enddo
       do i = ile2,im-1
         flaj(i,jm) = (p(i,jm)+p(i+1,jm))*ajx(i,jm)*0.5
       enddo
       do i = ite1,ite2-1
         flaj(i,1)  = (p(i,1)+p(i+1,1))*ajx(i,1)*0.5
       enddo
      endif

      if(m.eq.3) then

c     boundary condition at walls j=1, j=jm

       do i = 1,ile1-1
         flaj(i,jm) = (p(i,jm)+p(i+1,jm))*ajy(i,jm)*0.5
       enddo
       do i = ile2,im-1
         flaj(i,jm) = (p(i,jm)+p(i+1,jm))*ajy(i,jm)*0.5
       enddo
       do i = ite1,ite2-1
         flaj(i,1)  = (p(i,1)+p(i+1,1))*ajy(i,1)*0.5
       enddo
      endif


c     fluxbalance for single volume

      do  j = 1,jm-1
       do i = 1,im-1

         fblc(i,j) = flai(i,j)-flai(i+1,j  )
     &              +flaj(i,j)-flaj(i  ,j+1)
       enddo
      enddo

c     h-block

c     calculate flux through area type ai

      do  j = 1,jmh-1
       do i = 1,imh

         flaih(i,j) = (ar1h(i,j  )
     &                +ar1h(i,j+1))
     &                *aixh(i,j)*0.5
     &               +(ar2h(i,j  )
     &                +ar2h(i,j+1))
     &                *aiyh(i,j)*0.5
       enddo
      enddo

c     calculate flux through area type aj

      do  j = 1,jmh
       do i = 1,imh-1

         flajh(i,j) = (ar1h(i  ,j)
     &                +ar1h(i+1,j))
     &                *ajxh(i,j)*0.5
     &               +(ar2h(i  ,j)
     &                +ar2h(i+1,j))
     &                *ajyh(i,j)*0.5
       enddo
      enddo

c     boundary condition at walls j=1, j=jm

      if(m.eq.1.or.m.eq.4) then
       do i = 1,imh-1
         flajh(i,1)   = 0.0
         flajh(i,jmh) = 0.0
       enddo
      endif

      if(m.eq.2) then

c     boundary condition at walls j=1, j=jm

       do i = 1,imh-1
         flajh(i,1)   = (ph(i,1)  +ph(i+1,1))  *ajxh(i,1) *0.5
         flajh(i,jmh) = (ph(i,jmh)+ph(i+1,jmh))*ajxh(i,jmh)*0.5
       enddo
      endif

      if(m.eq.3) then

c     boundary condition at walls j=1, j=jm

       do i = 1,imh-1
         flajh(i,1)   = (ph(i,1)  +ph(i+1,1))  *ajyh(i,1) *0.5
         flajh(i,jmh) = (ph(i,jmh)+ph(i+1,jmh))*ajyh(i,jmh)*0.5
       enddo
      endif


c     fluxbalance for single volume

      do  j = 1,jm-1
       do i = 1,im-1

         fblc(i,j) = flai(i,j)-flai(i+1,j  )
     &              +flaj(i,j)-flaj(i  ,j+1)
       enddo
      enddo

      do  j = 1,jmh-1
       do i = 1,imh-1

         fblch(i,j) = flaih(i,j)-flaih(i+1,j  )
     &               +flajh(i,j)-flajh(i  ,j+1)
       enddo
      enddo

c     sum up for flux balances of node control volumes
      call sfl(fblc,fblp,fblch,fblph,0)


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

      do  j = 1,jmh-1
       do i = 1,imh
        aixh(i,j) = yh(i,j+1) - yh(i,j)
       enddo
      enddo

      do  j = 1,jmh-1
       do i = 1,imh
        aiyh(i,j) = -(xh(i,j+1) - xh(i,j))
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

      do  j = 1,jmh
       do i = 1,imh-1
        ajxh(i,j) = -(yh(i+1,j) - yh(i,j))
       enddo
      enddo

      do  j = 1,jmh
       do i = 1,imh-1
        ajyh(i,j) = xh(i+1,j) - xh(i,j)
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

      do  j = 1,jmh-1
       do i = 1,imh-1

          call gfl(xh  (i  ,j  ),xh  (i+1,j  ),
     &             xh  (i+1,j+1),xh  (i  ,j+1),
     &             yh  (i  ,j  ),yh  (i+1,j  ),
     &             yh  (i+1,j+1),yh  (i  ,j+1),
     &             volch(i,j))
       enddo
      enddo

c     sum up for node control volumes
      call sfl(volc,volp,volch,volph,0)

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
         pinold(i) = 0.0
        enddo
      enddo

      do  j = 1,jmh
       do i = 1,imh
         xh   (i,j) = 0.0
         yh   (i,j) = 0.0
         u1h  (i,j) = 0.0
         u2h  (i,j) = 0.0
         u3h  (i,j) = 0.0
         u4h  (i,j) = 0.0
         vxh  (i,j) = 0.0
         vyh  (i,j) = 0.0
         ph   (i,j) = 0.0
         qch  (i,j) = 0.0
         u1eh (i,j) = 0.0
         u2eh (i,j) = 0.0
         u3eh (i,j) = 0.0
         u4eh (i,j) = 0.0
         u1mh (i,j) = 0.0
         u2mh (i,j) = 0.0
         u3mh (i,j) = 0.0
         u4mh (i,j) = 0.0

         steph(i,j) = 0.0
         flaih(i,j) = 0.0
         flajh(i,j) = 0.0
         rlah (i,j) = 0.0
         fblch(i,j) = 0.0
         fblph(i,j) = 0.0
         aixh (i,j) = 0.0
         aiyh (i,j) = 0.0
         ajxh (i,j) = 0.0
         ajyh (i,j) = 0.0
         ar1h (i,j) = 0.0
         ar2h (i,j) = 0.0
        enddo
         pinold(j)  = 0.0
      enddo


      return
      end
      subroutine inlet
c     ================

      include 'param2D.h'
      include 'common2D.h'

      dimension ain(jdh)

c     boundary condition at inflow

c     from isentropic relations
      do j = 1,jmh
       ph(1,j)   = (1.0-rfpin)*ph(2,j)+rfpin*pinold(j)
       pinold(j) = ph(1,j)
       ts1       = tt1*(ph(1,j)/pt1)**rgzucp
       rma       = sqrt(2.0*cvzurg*(tt1/ts1-1.0))
       vin       = rma*sqrt(cpzucv*rg*ts1)
       u1h(1,j)  = ph(1,j)/(rg*ts1)

c      inflow direction
       vxh(1,j)   = vin
       vyh(1,j)   = 0.0
       u2h(1,j)   = u1h(1,j)*vxh(1,j)
       u3h(1,j)   = u1h(1,j)*vyh(1,j)
       u4h(1,j)   = ph(1,j)/rgzucv
     &             +0.5*u1h(1,j)*(vxh(1,j)**2+vyh(1,j)**2)
      enddo

      sain    = 0.0
      pinsum  = 0.0

      do j = 1,jmh-1
       ain(j)  = sqrt(ajxh(1,j)**2+ajyh(1,j)**2)
       sain    = sain+ain(j)
       pinsum  = pinsum + 0.5*(ph(1,j+1)+ph(1,j))*ain(j)
      enddo

      pinsum = pinsum/sain

      write(6,*) 'mean static pressure inlet: ',pinsum


      return
      end
      subroutine lax(q,fbal,qe,qm,qh,fbalh,qeh,qmh)
c     =============================================


      include 'param2D.h'
      include 'common2D.h'

      dimension q(id,jd),qe(id,jd),qm(id,jd)
      dimension fbal(id,jd)
      dimension qh(idh,jdh),qeh(idh,jdh),qmh(idh,jdh)
      dimension fbalh(idh,jdh)

c     volume averages single volumes

      do  j = 1,jm-1
       do i = 1,im-1

         qe(i,j) = 0.25*(q(i  ,j  )+q(i  ,j+1)
     &                  +q(i+1,j+1)+q(i+1,j  ))

       enddo
      enddo

      do  j = 1,jmh-1
       do i = 1,imh-1

         qeh(i,j) = 0.25*(qh(i  ,j  )+qh(i  ,j+1)
     &                   +qh(i+1,j+1)+qh(i+1,j  ))

       enddo
      enddo

c     sum up volume averages to nodes
      call sfl(qe,qm,qeh,qmh,1)

c     update conservative variable

c     regular field
      do  j = 2,jm-1
       do i = 1,im

         q(i,j)  = (1.0-rfd)*q(i,j)+rfd*qm(i,j)
     &                +ft*fbal(i,j)*step(i,j)

       enddo
      enddo

c     outer boundary and block connectivity
      do i = 1,im
         q(i,jm)  = (1.0-rfd)*q(i,jm)+rfd*qm(i,jm)
     &                +ft*fbal(i,jm)*step(i,jm)
      enddo

c     inner boundary
      do i = 1,ite2-1
         q(i,1)  = (1.0-rfd)*q(i,1)+rfd*qm(i,1)
     &                +ft*fbal(i,1)*step(i,1)
      enddo

c     block connectivity
        i2 = ite1+1
      do i = ite2,im
        i2 = i2-1
        q(i,1) = q(i2,1)
      enddo

c     h-block
      do  j = 1,jmh
       do i = 1,imh-1

         qh(i,j)  = (1.0-rfd)*qh(i,j)+rfd*qmh(i,j)
     &                +ft*fbalh(i,j)*steph(i,j)

       enddo
      enddo

      do j = 1,jmh

       qh(imh,j) = q(ile2-j+1,jm)

      enddo


      return
      end
      subroutine outlet
c     =================

      include 'param2D.h'
      include 'common2D.h'

      dimension aout(jd)

c     boundary condition at outflow

       do j = 1,jm

        p(1 ,j) = (1.0-rfpout)*p(1 ,j)+rfpout*pout
        p(im,j) = (1.0-rfpout)*p(im,j)+rfpout*pout

       enddo


      asout  = 0.0
      pimsum = 0.0

      do j = 1,jm-1

        aout(j)  = sqrt(aix(1,j)**2+aiy(1,j)**2)
        asout    = asout+aout(j)

        pimsum = pimsum+0.5*(p(1,j  )
     &                      +p(1,j+1))*aout(j)

      enddo

      do j = 1,jm-1

        aout(j)  = sqrt(aix(im,j)**2+aiy(im,j)**2)
        asout    = asout+aout(j)

        pimsum = pimsum+0.5*(p(im,j  )
     &                      +p(im,j+1))*aout(j)

      enddo

      pimsum   = pimsum/asout

      write(6,*)  'mean static pressure outlet:',pimsum,' pout:',pout
c      write(11,*) 'mean static pressure outlet:',pimsum,' pout:',pout


      return
      end
      subroutine post
c     ===============

      include 'param2D.h'
      include 'common2D.h'

      dimension flowi(id),ain(jdh)

c     mean static pressure at inlet
      pinsum = 0.0
      sain   = 0.0
      do j = 1,jmh-1
       ain(j)  = sqrt(aixh(1,j)**2+aiyh(1,j)**2)
       sain    = sain+ain(j)
       pinsum  = pinsum + 0.5*(ph(1,j+1)+ph(1,j))*ain(j)
      enddo

      pinsum = pinsum/sain

      write(11,*) 'mean static pressure inlet: ',pinsum
      write(11,*) '   '

c     calculate massflow

      do  j = 1,jm
       do i = 1,im

          ar1(i,j) = u1(i,j)*vx(i,j)
          ar2(i,j) = u1(i,j)*vy(i,j)

       enddo
      enddo

      do  j = 1,jmh
       do i = 1,imh

          ar1h(i,j) = u1h(i,j)*vxh(i,j)
          ar2h(i,j) = u1h(i,j)*vyh(i,j)

       enddo
      enddo

      call flbal(1)

      write(11,*) '   '
      write(11,*) 'massflow ratio at i-planes: '
      write(11,*) '   '

c     massflow at inlet
      flowin = 0.0
      do   j = 1,jmh-1
       flowin = flowin + flaih(1,j)
      enddo


      do  i = 1,imh
       flowih = 0.0
       do j = 1,jmh-1
        flowih = flowih+flaih(i,j)
       enddo
       write(11,*) i,flowih,flowih/flowin
      enddo

          i2 = im+1
      do   i = 1,ile1
          i2 = i2-1
       flowi(i) = 0.0
       do j = 1,jm-1
        flowi(i) = flowi(i)+flai(i,j)-flai(i2,j)
       enddo
       write(11,*) i,flowi(i),flowi(i)/flowin
      enddo

c     massflow at outlet
      flowot = flowi(1)

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

      do j = 1,jmh-1

c       inlet
        ts        = (ph(1,j+1)+ph(1,j))/((u1h(1,j+1)+u1h(1,j))*rg)
        aq        = cpzucv*rg*ts
        rmaq      = ((0.5*(vxh(1,j+1)+vxh(1,j)))**2
     &              +(0.5*(vyh(1,j+1)+vyh(1,j)))**2)/aq
        tt        = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
        pt        = 0.5*(ph(1,j+1)+ph(1,j))/((ts/tt)**cpzurg)
        ptinlt    = ptinlt+pt*flaih(1,j)
        ttinlt    = ttinlt+tt*flaih(1,j)

      enddo


      do j = 1,jm-1

c       outlet 1
        ts        = (p(1,j+1)+p(1,j))/((u1(1,j+1)+u1(1,j))*rg)
        aq        = cpzucv*rg*ts
        rmaq      = ((0.5*(vx(1,j+1)+vx(1,j)))**2
     &              +(0.5*(vy(1,j+1)+vy(1,j)))**2)/aq
        tt        = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
        pt        = p(1,j)/((ts/tt)**cpzurg)
        ptout     = ptout+pt*flai(1,j)
        ttout     = ttout+tt*flai(1,j)

c       outlet 2
        ts        = (p(im,j+1)+p(im,j))/((u1(im,j+1)+u1(im,j))*rg)
        aq        = cpzucv*rg*ts
        rmaq      = ((0.5*(vx(im,j+1)+vx(im,j)))**2
     &              +(0.5*(vy(im,j+1)+vy(im,j)))**2)/aq
        tt        = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
        pt        = p(im,j)/((ts/tt)**cpzurg)
        ptout     = ptout-pt*flai(im,j)
        ttout     = ttout-tt*flai(im,j)

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
      write(11,*) 'Forces:   '

      fdrag = 0.0
      flift = 0.0
      do i = ite1,ite2-1
       fdrag = fdrag-0.5*(p(i,1)+p(i+1,1))*ajx(i,1)*0.003
       flift = flift+0.5*(p(i,1)+p(i+1,1))*ajy(i,1)*0.003
      enddo

      xmom1 = 0.0
      do i = ile1,ile2-1
       xmom1 = xmom1+0.5*(u1(i+1,jm)+u1(i,jm))
     &             *(0.5*(vx(i+1,jm)+vx(i,jm)))**2*ajx(i,jm)
      enddo

      xmom2 = 0.0
      do j = 1,jm-1
       xmom2 = xmom2+0.5*(u1(1,j+1)+u1(1,j))
     &             *(0.5*(vx(1,j+1)+vx(1,j)))**2*aix(1,j)

       xmom2 = xmom2-0.5*(u1(im,j+1)+u1(im,j))
     &             *(0.5*(vx(im,j+1)+vx(im,j)))**2*aix(im,j)
      enddo

      write(11,*) '   '
      write(11,'(A17,F12.5)') 'lift       = ',flift
      write(11,'(A17,F12.5)') 'drag       = ',fdrag
      write(11,'(A17,F12.5)') 'xmom_in    = ',xmom1
      write(11,'(A17,F12.5)') 'xmom_out   = ',xmom2
      write(11,*) '   '

c     blade surface pressure distribution
      do i = ite1,ite2
       write(15,*) x(i,1),p(i,1)
      enddo

c      write(11,*) 'write vx'
c      call write_q(vx)

c      write(11,*) 'write p'
c      call write_q(p)

c      write(11,*) 'write vxh'
c      call write_qh(vxh)

c      write(11,*) 'write ph'
c      call write_qh(ph)

      return
      end
      subroutine readst
c     =================

      include 'param2D.h'
      include 'common2D.h'

      character*80 titel


      namelist /l1/ nmax,nrfd,istart,
     &              rfd,ft,rg,cpzucv,ermax,
     &              rfpin,rfpout,pt1,tt1,alfa,pout,
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
      write(11,'(A9,F10.5)') 'u1star = ', u1star
      write(11,'(A9,F10.3)') 'pt1    = ', pt1
      write(11,'(A9,F10.3)') 'tt1    = ', tt1
      write(11,'(A9,F10.3)') 'alfa   = ', alfa
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
      write(6,'(A9,F10.5)') 'u1star = ', u1star
      write(6,'(A9,F10.3)') 'pt1    = ', pt1
      write(6,'(A9,F10.3)') 'tt1    = ', tt1
      write(6,'(A9,F10.3)') 'alfa   = ', alfa
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

c      write(11,*) 'read_coord','im=',im,'jm=',jm
      do  j = 1,jm
       do i = 1,im
        read(14,*) x(i,j),y(i,j)
       enddo
      enddo

c     generate h-inlet block
      do  j = 1,jmh
       do i = 1,imh
        xh(i,j) = x(ile1,  jm)-real(imh-i)
     &          *(x(ile1-1,jm)-x(ile1,jm))
        yh(i,j) = y(ile2-j+1,jm)
       enddo
      enddo

      return
      end
      subroutine sfl(dudt,dqdt,dudth,dqdth,m)
c     =======================================

      include 'param2D.h'
      include 'common2D.h'

c     sum up single balances to control volume

      dimension dudt(id,jd),dqdt(id,jd)
      dimension dudth(idh,jdh),dqdth(idh,jdh)

c     simple summing up or average
      if(m.eq.1) then
       f4 = 0.25
       f3 = 0.333333
       f2 = 0.5
      else
       f4 = 1.0
       f3 = 1.0
       f2 = 1.0
      endif

c     regular field
       do  j = 2,jm-1
        dO i = 2,im-1

         dqdt(i,j)  = (dudt(i-1,j-1)
     &                +dudt(i-1,j  )
     &                +dudt(i  ,j  )
     &                +dudt(i  ,j-1))*f4
       enddo
      enddo

       do  j = 2,jmh-1
        do i = 2,imh-1

         dqdth(i,j)  = (dudth(i-1,j-1)
     &                 +dudth(i-1,j  )
     &                 +dudth(i  ,j  )
     &                 +dudth(i  ,j-1))*f4
       enddo
      enddo

c     values at boundaries

c     walls
      do  i = 2,ile1-1
         dqdt(i,jm)  = (dudt(i-1,jm-1)
     &                 +dudt(i  ,jm-1))*f2
      enddo

      do  i = ile2+1,im-1
         dqdt(i,jm)  = (dudt(i-1,jm-1)
     &                 +dudt(i  ,jm-1))*f2
      enddo

      do  i = 2,imh-1
         dqdth(i,1)   = (dudth(i-1,1)
     &                  +dudth(i  ,1))*f2
         dqdth(i,jmh) = (dudth(i-1,jmh-1)
     &                  +dudth(i  ,jmh-1))*f2
      enddo

c     inlet
      do  j = 2,jmh-1
       dqdth(1,j)  = (dudth(1,j-1)
     &               +dudth(1,j))*f2
      enddo

c     outlet
      do  j = 2,jm-1
        dqdt(1,j)  = (dudt(1,j-1)
     &               +dudt(1,j))*f2

        dqdt(im,j) = (dudt(im-1,j-1)
     &               +dudt(im-1,j))*f2
      enddo

c     profile wall
      do  i = ite1+1,ite2-1
        dqdt(i,1)  = (dudt(i-1,1)
     &               +dudt(i  ,1))*f2
      enddo

c     block connectivity c-h
         jh = jmh
      do  i = ile1+1,ile2-1
         jh = jh-1
        dqdt(i,jm)  = (dudt(i-1,jm-1)
     &                +dudt(i  ,jm-1)
     &                +dudth(imh-1,jh)
     &                +dudth(imh-1,jh-1))*f4
      enddo

      dqdt(ile1,jm) = (dudt(ile1-1,jm-1)
     &                +dudt(ile1  ,jm-1)
     &                +dudth(imh-1,jmh-1))*f3
      dqdt(ile2,jm) = (dudt(ile2-1,jm-1)
     &                +dudt(ile2  ,jm-1)
     &                +dudth(imh-1,1))*f3

c     block connectivity c-c
c      do  i = 2,ite1-1
      do  i = 2,ite1
         i2 = im-i+1

        dqdt(i,1) = (dudt(i2-1,1)
     &              +dudt(i -1,1)
     &              +dudt(i   ,1)
     &              +dudt(i2-1,1))*f4
      enddo

c     edges inlet
      dqdth(1,1)   = dudth(1,1)
      dqdth(1,jmh) = dudth(1,jmh-1)

c     edges outlet
      dqdt( 1,jm)  =  dudt(1   ,jm-1)
      dqdt(im,jm)  =  dudt(im-1,jm-1)

      dqdt(1,1)    = (dudt(1,1)
     &               +dudt(im-1,1))*f2

c     trailing edge
c      dqdt(ite1,1) = (dudt(ite1-1,1)
c     &               +dudt(ite2  ,1))*f2


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

        do  j = 1,jmh
         do i = 1,imh

          u1h(i,j) = u1star
          ph (i,j) = pstart
          vxh(i,j) = vxstar
          vyh(i,j) = vystar

        enddo
       enddo


      do  j = 1,jm
       do i = 1,im

         u2(i,j) = u1(i,j)*vx(i,j)
         u3(i,j) = u1(i,j)*vy(i,j)
         u4(i,j) = p(i,j)/rgzucv+0.5*u1(i,j)*(vx(i,j)**2+vy(i,j)**2)

       enddo
      enddo

      do  j = 1,jmh
       do i = 1,imh

         u2h(i,j) = u1h(i,j)*vxh(i,j)
         u3h(i,j) = u1h(i,j)*vyh(i,j)
         u4h(i,j) = ph(i,j)/rgzucv
     &              +0.5*u1h(i,j)*(vxh(i,j)**2+vyh(i,j)**2)
       enddo
      enddo

      do j = 1,jmh
       pinold(j) = ph(1,j)
      enddo


c      write(11,*) 'start:'
c      call write_prime

c     write bilances
c      call write_fluxes

      return
      end
      subroutine update
c     =================

c     Iteration of conservative variables to convergence

      include 'param2D.h'
      include 'common2D.h'

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

       do  j = 1,jmh
        do i = 1,imh

          ar1h(i,j) = u1h(i,j)*vxh(i,j)
          ar2h(i,j) = u1h(i,j)*vyh(i,j)

        enddo
       enddo

c      calculate flux balances
       call flbal(1)

c      update u1
       call lax(u1,fblp,u1e,u1m,u1h,fblph,u1eh,u1mh)


c      massflow at inlet
       flowin = 0.0
       do   j = 1,jmh-1
        flowin = flowin + flaih(1,j)
       enddo
c      write(11,*) 'massflow at inlet: ',flowin


c      massflow at outlet
       flowot = 0.0
       do j = 1,jm-1
         flowot = flowot+flai(1,j)-flai(im,j)
       enddo

       flrate = flowot/flowin
c      write(11,*) flowin,flowot,flrate


       write(6 ,*)  'n = ',nstep,'massflow out/in',flrate,' rfd = ',rfd
c      write(11,*)  'n = ',nstep,'massflow out/in',flrate,' rfd = ',rfd



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


       do  j = 1,jmh
        do i = 1,imh

          ar1h(i,j) = u2h(i,j)*vxh(i,j)+ph(i,j)
          ar2h(i,j) = u2h(i,j)*vyh(i,j)

        enddo
       enddo

c      calculate flux balances
       call flbal(2)

c      update u2
       call lax(u2,fblp,u2e,u2m,u2h,fblph,u2eh,u2mh)


c      update y-momentum
c      *****************


c      calculate flux variables

       do  j = 1,jm
        do i = 1,im

          ar1(i,j) = u3(i,j)*vx(i,j)
          ar2(i,j) = u3(i,j)*vy(i,j)+p(i,j)

        enddo
       enddo

       do  j = 1,jmh
        do i = 1,imh

          ar1h(i,j) = u3h(i,j)*vxh(i,j)
          ar2h(i,j) = u3h(i,j)*vyh(i,j)+ph(i,j)

        enddo
       enddo

c      calculate flux balances
       call flbal(3)

c      update u3
       call lax(u3,fblp,u3e,u3m,u3h,fblph,u3eh,u3mh)


c      update total energy
c      *******************

c      calculate flux variables

       do  j = 1,jm
        do i = 1,im

          ar1(i,j) = (u4(i,j)+p(i,j))*vx(i,j)
          ar2(i,j) = (u4(i,j)+p(i,j))*vy(i,j)

        enddo
       enddo

       do  j = 1,jmh
        do i = 1,imh

          ar1h(i,j) = (u4h(i,j)+ph(i,j))*vxh(i,j)
          ar2h(i,j) = (u4h(i,j)+ph(i,j))*vyh(i,j)

        enddo
       enddo


c      calculate flux balances
       call flbal(4)

c      update u4
       call lax(u4,fblp,u4e,u4m,u4h,fblph,u4eh,u4mh)


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

       do  j = 1,jmh
        do i = 1,imh-1

          qch(i,j) = amax1(0.00001,abs(vxh(i,j)))

          vxh(i,j) = u2h(i,j)/u1h(i,j)
          vyh(i,j) = u3h(i,j)/u1h(i,j)

          ph (i,j) = (u4h(i,j) - 0.5*(vxh(i,j)**2+vyh(i,j)**2)
     &              *u1h(i,j))*rgzucv

        enddo
       enddo

c      connectivity c-h
       do  j = 1,jmh

        vxh(imh,j) = vx(ile2-j+1,jm)
        vyh(imh,j) = vy(ile2-j+1,jm)
        ph (imh,j) = p (ile2-j+1,jm)

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
           write(6 ,*) 'vxc:',vx(i,j),'vyc:',vy(i,j),'qc:',qc(i,j)
           write(6 ,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(11,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(6 ,*) 'erijc:',erij,' Demo3D diverges! at: ',i,j
           write(11,*) 'erijc:',erij,' Demo3D diverges!'
           goto 999
          endif

          if(erij.gt.evmax) then
           evmax = erij
           iemax = i
           jemax = j
          endif
        enddo
       enddo



       do  j = 1,jmh
        do i = 1,imh-1
          erij = 1.0-sqrt(vxh(i,j)**2+vyh(i,j)**2)/qch(i,j)
          if(abs(erij).gt.1000.0) then
           write(6 ,*) 'ih,jh ',i,j
           write(6 ,*) 'vxh:',vxh(i,j),'vyh:',vyh(i,j),'qch:',qch(i,j)
           write(6 ,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(11,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(6 ,*) 'erijh:',erij,' Demo3D diverges! at: ',i,j
           write(11,*) 'erijh:',erij,' Demo3D diverges!'
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
c       write(11, *) 'n = ',n,' at',iemax,jemax,'evmax = ',evmax
       write(13, *) evmax

c      check convergence
       if(abs(flrate-1.0).lt.ermax.and.n.gt.1000) then
        nconv = nconv +1
        write(6 ,*) 'nconv = ',nconv
        write(11,*) 'nconv = ',nconv
        if(nconv.eq.10) then
         write(6 ,*) 'job converged '
         goto 999
        endif
       endif


      enddo


  999 continue
      write(11,*) ' '
      write(11,*) 'nsteps calculated = ',nstep
      write(11,*) ' '
      return
      end
      subroutine write_boundary
c     =========================

      include 'cgnslib_f.h'
      include 'param2D.h'

      character*32 basename,zonename,famname

      dimension isize(3,3),ipnts(3,3)

      call cg_open_f
     &    ('Demo2D_result_rel.cgns',CG_MODE_MODIFY,index_file,ier)
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

c     write Outlet panel

      famname = 'Outlet'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family Outlet'
      write(6,*) 'Outlet: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f


c     index space for Outlet
      ipnts(1,1) = ilo
      ipnts(2,1) = jlo
      ipnts(3,1) = klo
      ipnts(1,2) = ilo
      ipnts(2,2) = jhi
      ipnts(3,2) = khi

      call cg_boco_write_f(index_file,index_base,index_zone,'Outlet1',
     &             BCTunnelOutflow,PointRange,2,ipnts,index_bc,ier)
      ipnts(1,1) = ihi
      ipnts(2,1) = jlo
      ipnts(3,1) = klo
      ipnts(1,2) = ihi
      ipnts(2,2) = jhi
      ipnts(3,2) = khi

      call cg_boco_write_f(index_file,index_base,index_zone,'Outlet2',
     &             BCTunnelOutflow,PointRange,2,ipnts,index_bc,ier)

c     write Blade panel

      famname = 'Blade'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family Blade'
      write(6,*) 'Blade: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

c     index space for Blade
      ipnts(1,1) = ite1
      ipnts(2,1) = jlo
      ipnts(3,1) = klo
      ipnts(1,2) = ite2
      ipnts(2,2) = jlo
      ipnts(3,2) = khi

      call cg_boco_write_f(index_file,index_base,index_zone,'Blade',
     &             BCWallViscous,PointRange,2,ipnts,index_bc,ier)

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

      call cg_boco_write_f(index_file,index_base,index_zone,'Inlet',
     &             BCTunnelInflow,PointRange,2,ipnts,index_bc,ier)


c     write WallHub panel

      famname = 'WallHubH'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family WallHubH'
      write(6,*) 'WallHubH: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

c     index space for WallHubH
      ipnts(1,1) = ilo
      ipnts(2,1) = jlo
      ipnts(3,1) = klo
      ipnts(1,2) = ihi
      ipnts(2,2) = jlo
      ipnts(3,2) = khi

      call cg_boco_write_f(index_file,index_base,index_zone,'WallHubH',
     &             BCWallViscous,PointRange,2,ipnts,index_bc,ier)

c     write WallTip panel

      famname = 'WallTipH'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family WallTipH'
      write(6,*) 'WallTipH: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

c     index space for WallTipH
      ipnts(1,1) = ilo
      ipnts(2,1) = jhi
      ipnts(3,1) = klo
      ipnts(1,2) = ihi
      ipnts(2,2) = jhi
      ipnts(3,2) = khi

      call cg_boco_write_f(index_file,index_base,index_zone,'WallTipH',
     &             BCWallViscous,PointRange,2,ipnts,index_bc,ier)

c     write WallLeft panel

      famname = 'WallLeftH'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family WallLeftH'
      write(6,*) 'WallLeftH: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

c     index space for WallLeft
      ipnts(1,1) = ilo
      ipnts(2,1) = jlo
      ipnts(3,1) = klo
      ipnts(1,2) = ihi
      ipnts(2,2) = jhi
      ipnts(3,2) = klo

      call cg_boco_write_f(index_file,index_base,index_zone,'WallLeftH',
     &             BCWallViscous,PointRange,2,ipnts,index_bc,ier)


c     write WallRight panel

      famname = 'WallRightH'
      call cg_family_write_f(index_file,index_base,famname,
     &                       index_fam,ier)
      write(6,*) 'create family WallRightH'
      write(6,*) 'WallRightH: index_fam =',index_fam
      if(ier.ne.CG_OK) call cg_error_exit_f

c     index space for WallRight
      ipnts(1,1) = ilo
      ipnts(2,1) = jlo
      ipnts(3,1) = khi
      ipnts(1,2) = ihi
      ipnts(2,2) = jhi
      ipnts(3,2) = khi

      call cg_boco_write_f(index_file,index_base,index_zone,
     & 'WallRightH',BCWallViscous,PointRange,2,ipnts,index_bc,ier)

c close CGNS file
      call cg_close_f(index_file,ier)
      write(6,'('' Successfully wrote boundary panels to cgns'')')

      return
      end
      subroutine write_rel_cgns
c     =========================

      include 'cgnslib_f.h'
      include 'param2D.h'
      include 'common2D.h'

      dimension x3d (id,jd,kd),y3d  (id,jd,kd),z3d (id,jd,kd)
      dimension vx3d(id,jd,kd),vy3d (id,jd,kd),vz3d(id,jd,kd)
      dimension p3d (id,jd,kd),u13d (id,jd,kd)
      dimension x3dh (idh,jdh,kdh),y3dh (idh,jdh,kdh),z3dh (idh,jdh,kdh)
      dimension vx3dh(idh,jdh,kdh),vy3dh(idh,jdh,kdh),vz3dh(idh,jdh,kdh)
      dimension p3dh (idh,jdh,kdh),u13dh(idh,jdh,kdh)

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

      call cg_open_f
     &    ('Demo2D_result_rel.cgns',CG_MODE_WRITE,index_file,ier)

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


      do   k = 1,km
       do  j = 1,jmh
        do i = 1,imh

         x3dh (i,j,k) = xh(i,j)
         y3dh (i,j,k) = yh(i,j)
         z3dh (i,j,k) = real(k)*0.001
         vx3dh(i,j,k) = vxh(i,j)
         vy3dh(i,j,k) = vyh(i,j)
         vz3dh(i,j,k) = 0.0
         p3dh (i,j,k) = ph(i,j)
         u13dh(i,j,k) = u1h(i,j)

        enddo
       enddo
      enddo

      zonename = 'Block 2'

c vertex size
      isize(1,1) = imh
      isize(2,1) = jmh
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
     &                    RealSingle,'CoordinateX',x3dh,index_coord,ier)
      write(6,*) 'index_coord X =',index_coord

      call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateY',y3dh,index_coord,ier)
      write(6,*) 'index_coord Y =',index_coord

      call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateZ',z3dh,index_coord,ier)
      write(6,*) 'index_coord Z =',index_coord


c    define flow solution node name (user can give any name)
      solname = 'FlowSolution'

c    create flow solution node
      call cg_sol_write_f(index_file,index_base,index_zone,solname,
     &                    Vertex,index_flow,ier)

c    write flow solution (user must use SIDS-standard names here)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Density',u13dh,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Pressure',p3dh,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityX',vx3dh,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityY',vy3dh,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityZ',vz3dh,index_field,ier)



c close CGNS file
      call cg_close_f(index_file,ier)
      write(6,'('' Successfully wrote grid to file
     &           Demo2D_result_rel.cgns'')')



      return
      end
      subroutine write_abs_cgns
c     =========================

      include 'cgnslib_f.h'
      include 'param2D.h'
      include 'common2D.h'

      dimension x3d (id,jd,kd),y3d  (id,jd,kd),z3d (id,jd,kd)
      dimension vx3d(id,jd,kd),vy3d (id,jd,kd),vz3d(id,jd,kd)
      dimension p3d (id,jd,kd),u13d (id,jd,kd)
      dimension x3dh (idh,jdh,kdh),y3dh (idh,jdh,kdh),z3dh (idh,jdh,kdh)
      dimension vx3dh(idh,jdh,kdh),vy3dh(idh,jdh,kdh),vz3dh(idh,jdh,kdh)
      dimension p3dh (idh,jdh,kdh),u13dh(idh,jdh,kdh)

      dimension isize(3,3)

      character basename*32,zonename*32,solname*32

      km = 3

      do   k = 1,km
       do  j = 1,jm
        do i = 1,im

         x3d (i,j,k) = x(i,j)
         y3d (i,j,k) = y(i,j)
         z3d (i,j,k) = real(k)*0.001
c        absolut vx
         vx3d(i,j,k) = vx(i,j)-vxh(1,(jmh-1)/2+1)
         vy3d(i,j,k) = vy(i,j)
         vz3d(i,j,k) = 0.0
         ts          = p(i,j)/(rg*u1(i,j))
         pt          = p(i,j)*(tt1/ts)**cpzurg
         p3d (i,j,k) = pt
         u13d(i,j,k) = u1(i,j)

        enddo
       enddo
      enddo


c     WRITE X, Y, Z GRID POINTS TO CGNS FILE
c     open CGNS file for write

      call cg_open_f
     &    ('Demo2D_result_abs.cgns',CG_MODE_WRITE,index_file,ier)

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



      do   k = 1,km
       do  j = 1,jmh
        do i = 1,imh

         x3dh (i,j,k) = xh(i,j)
         y3dh (i,j,k) = yh(i,j)
         z3dh (i,j,k) = real(k)*0.001
         vx3dh(i,j,k) = vxh(i,j)-vxh(1,(jmh-1)/2+1)
         vy3dh(i,j,k) = vyh(i,j)
         vz3dh(i,j,k) = 0.0
         ts           = ph(i,j)/(rg*u1h(i,j))
         pt           = ph(i,j)*(tt1/ts)**cpzurg
         p3dh (i,j,k) = pt
         u13dh(i,j,k) = u1h(i,j)

        enddo
       enddo
      enddo

      zonename = 'Block 2'

c vertex size
      isize(1,1) = imh
      isize(2,1) = jmh
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
     &                    RealSingle,'CoordinateX',x3dh,index_coord,ier)
      write(6,*) 'index_coord X =',index_coord

      call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateY',y3dh,index_coord,ier)
      write(6,*) 'index_coord Y =',index_coord

      call cg_coord_write_f(index_file,index_base,index_zone,
     &                    RealSingle,'CoordinateZ',z3dh,index_coord,ier)
      write(6,*) 'index_coord Z =',index_coord


c    define flow solution node name (user can give any name)
      solname = 'FlowSolution'

c    create flow solution node
      call cg_sol_write_f(index_file,index_base,index_zone,solname,
     &                    Vertex,index_flow,ier)

c    write flow solution (user must use SIDS-standard names here)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Density',u13dh,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'Pressure',p3dh,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityX',vx3dh,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityY',vy3dh,index_field,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_flow,
     & RealSingle,'VelocityZ',vz3dh,index_field,ier)



c close CGNS file
      call cg_close_f(index_file,ier)
      write(6,'('' Successfully wrote grid to file
     &             Demo2D_result_abs.cgns'')')



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


       do  j = 1,jmh-1
        do i = 1,imh-1

         ar1h(i,j)  =  aixh(i+1,j  )-aixh(i,j)
     &               + ajxh(i  ,j+1)-ajxh(i,j)
     &               + aiyh(i+1,j  )-aiyh(i,j)
     &               + ajyh(i  ,j+1)-ajyh(i,j)
       enddo
      enddo

      write(11,*) 'area aixh,aiyh'
      do  i = 1,imh
       do j = 1,jmh-1
         write(11,*) i,j,aixh(i,j),aiyh(i,j)
       enddo
      enddo

      write(11,*) 'area ajxh,ajyh'
      do  i = 1,imh-1
       do j = 1,jmh
         write(11,*) i,j,ajxh(i,j),ajyh(i,j)
       enddo
      enddo

      write(11,*) 'volc'
      do  i = 1,imh-1
       do j = 1,jmh-1
         write(11,*) i,j,volch(i,j)
       enddo
      enddo

      write(11,*) 'area check'

      do  j = 1,jmh-1
       do i = 1,imh-1
         write(11,*) i,j,ar1h(i,j)
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
          j = 1
       do i = 1,im,4
         write(11,*) i,j,vx(i,j),vy(i,j),p(i,j),u1(i,j)
       enddo
c      enddo

      write(11,*) 'write conservatives'

c      do  j = 1,jm
          j = 1
       do i = 1,im,4
         write(11,*) i,j,u1(i,j),u2(i,j),u3(i,j),u4(i,j)
       enddo
c      enddo

      return
      end
      subroutine write_test
c     =====================

      include 'param2D.h'
      include 'common2D.h'


      write(11,*) 'write test'
      write(11,*) 'im=',im,'jm=',jm

      i = ile1
      write(11,*) 'ile1,jm ',ile1,jm
      write(11,*) vx(i,jm),vy(i,jm),p(i,jm),u1(i,jm)

      i = (im-1)/2+1
      write(11,*) '(im-1)/2+1,jm '
      write(11,*) vx(i,jm),vy(i,jm),p(i,jm),u1(i,jm)

      i = ile2
      write(11,*) 'ile2,jm ',ile2
      write(11,*) vx(i,jm),vy(i,jm),p(i,jm),u1(i,jm)

      j = 1
      i = 1
      write(11,*) 'i=1,j=1 '
      write(11,*) vx(i,j),vy(i,j),p(i,j),u1(i,j)
      i = im
      write(11,*) 'i=im,j=1 '
      write(11,*) vx(i,j),vy(i,j),p(i,j),u1(i,j)

      j = jm
      i = 1
      write(11,*) 'i=1,j=jm '
      write(11,*) vx(i,j),vy(i,j),p(i,j),u1(i,j)
      i = im
      write(11,*) 'i=im,j=jm '
      write(11,*) vx(i,j),vy(i,j),p(i,j),u1(i,j)

      j = 1
      i = ite1
      write(11,*) 'i=ite1,j=1'
      write(11,*) vx(i,j),vy(i,j),p(i,j),u1(i,j)
      i = ite2
      write(11,*) 'i=ite2,j=1'
      write(11,*) vx(i,j),vy(i,j),p(i,j),u1(i,j)

      return
      end
      subroutine write_q(q)
c     =====================

      include 'param2D.h'
      include 'common2D.h'

      dimension q(id,jd)


      write(11,*) 'i,j=1,j=jmid,j=jm'
      jmid = (jm-1)/2+1
      do i = 1,im,2
      write(11,*) i,q(i,1),q(i,jmid),q(i,jm)
      enddo



      return
      end
      subroutine write_qh(qh)
c     =======================

      include 'param2D.h'
      include 'common2D.h'

      dimension qh(idh,jdh)


      write(11,*) 'ih,jh'

      jmid = (jmh-1)/2+1

      do  j = 1,jmh,4
       do i = 1,imh,4
       write(11,*) i,j,qh(i,j)
       enddo
      enddo



      return
      end

