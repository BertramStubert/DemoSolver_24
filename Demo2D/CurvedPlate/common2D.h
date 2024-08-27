      common /intg01/ im,jm,km,nmax,nstep,nrfd

      common /bcnd01/ pt1,tt1,pout
      common /bcnd02/ ptinlt,ttinlt,ptv,ttv

      common /read01/ rfd,ft,rg,ermax,rfpin,rfpout,ptvbl
      common /read03/ pstart,vxstar,vystar,u1star

      common /ther01/ cpzucv,cp,cv,rcp,rcv
      common /ther02/ rgzucv,rgzucp,cvzurg,cpzurg


      common /grid11/ x1  (id,jd),y1  (id,jd)

      common /aero11/ u11 (id,jd),u21 (id,jd),u31 (id,jd)
      common /aero12/ u41 (id,jd),u51 (id,jd)
      common /aero13/ u1e1(id,jd),u2e1(id,jd),u3e1(id,jd)
      common /aero14/ u4e1(id,jd),u5e1(id,jd)
      common /aero15/ u1m1(id,jd),u2m1(id,jd),u3m1(id,jd)
      common /aero16/ u4m1(id,jd),u5m1(id,jd)

      common /bcnd13/ pinol1(jd),vxinl1(jd),vximl1(jd)
      common /bcnd14/ poutl1(jd),u1iml1(jd),ptin1(jd)

      common /solv11/ rla1 (id,jd),step1(id,jd)
      common /solv12/ flai1(id,jd),flaj1(id,jd),flak1(id,jd)
      common /solv14/ fblc1(id,jd),fblp1(id,jd)

      common /velo11/ vx1 (id,jd),vy1 (id,jd),vz1 (id,jd)
      common /velo12/ p1  (id,jd),qc1 (id,jd)

      common /geom11/ aix1(id,jd),aiy1(id,jd),aiz1(id,jd)
      common /geom12/ ajx1(id,jd),ajy1(id,jd),ajz1(id,jd)
      common /geom14/ volc1(id,jd),volp1(id,jd)

      common /arra11/ ar11(id,jd),ar21(id,jd)


      common /grid21/ x2  (id,jd),y2  (id,jd)

      common /aero21/ u12 (id,jd),u22 (id,jd),u32 (id,jd)
      common /aero22/ u42 (id,jd),u52 (id,jd)
      common /aero23/ u1e2(id,jd),u2e2(id,jd),u3e2(id,jd)
      common /aero24/ u4e2(id,jd),u5e2(id,jd)
      common /aero25/ u1m2(id,jd),u2m2(id,jd),u3m2(id,jd)
      common /aero26/ u4m2(id,jd),u5m2(id,jd)

      common /bcnd23/ pinol2(jd),vxinl2(jd),vximl2(jd)
      common /bcnd24/ poutl2(jd),u1iml2(jd),ptin2(jd)

      common /solv21/ rla2 (id,jd),step2(id,jd)
      common /solv22/ flai2(id,jd),flaj2(id,jd),flak2(id,jd)
      common /solv24/ fblc2(id,jd),fblp2(id,jd)

      common /velo21/ vx2 (id,jd),vy2 (id,jd),vz2 (id,jd)
      common /velo22/ p2  (id,jd),qc2 (id,jd)

      common /geom21/ aix2(id,jd),aiy2(id,jd),aiz2(id,jd)
      common /geom22/ ajx2(id,jd),ajy2(id,jd),ajz2(id,jd)
      common /geom24/ volc2(id,jd),volp2(id,jd)

      common /arra21/ ar12(id,jd),ar22(id,jd)


