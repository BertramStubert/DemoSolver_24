      common /intg01/ im,jm,km,nmax,nstep,nrfd,istart
      common /grid01/ x  (id,jd),y  (id,jd)

      common /aero01/ u1 (id,jd),u2 (id,jd),u3 (id,jd)
      common /aero02/ u4 (id,jd),u5 (id,jd)
      common /aero03/ u1e(id,jd),u2e(id,jd),u3e(id,jd)
      common /aero04/ u4e(id,jd),u5e(id,jd)
      common /aero05/ u1m(id,jd),u2m(id,jd),u3m(id,jd)
      common /aero06/ u4m(id,jd),u5m(id,jd)

      common /bcnd01/ pt1,tt1,pout,pinold(jd),rfpin,rfpout
      common /bcnd02/ ptinlt,ttinlt,ptv,ttv

      common /solv01/ rla (id,jd),step(id,jd)
      common /solv02/ flai(id,jd),flaj(id,jd),flak(id,jd)
      common /solve4/ fblc(id,jd),fblp(id,jd)

      common /velo01/ vx (id,jd),vy (id,jd),vz (id,jd)
      common /velo02/ p  (id,jd),qc (id,jd)

      common /geom01/ aix(id,jd),aiy(id,jd),aiz(id,jd)
      common /geom02/ ajx(id,jd),ajy(id,jd),ajz(id,jd)
      common /geom04/ volc(id,jd),volp(id,jd)

      common /array1/ ar1(id,jd),ar2(id,jd)

      common /read01/ rfd,ft,rg,ermax
      common /read03/ pstart,vxstar,vystar,u1star

      common /ther01/ cpzucv,cp,cv,rcp,rcv
      common /ther02/ rgzucv,rgzucp,cvzurg,cpzurg





