      common /intg01/ im,jm,km,nmax,nstep,nrfd,istart
      common /grid01/ x  (id,jd,kd),y  (id,jd,kd),z  (id,jd,kd)

      common /aero01/ u1 (id,jd,kd),u2 (id,jd,kd),u3 (id,jd,kd)
      common /aero02/ u4 (id,jd,kd),u5 (id,jd,kd),p  (id,jd,kd)
      common /aero03/ u1e(id,jd,kd),u2e(id,jd,kd),u3e(id,jd,kd)
      common /aero04/ u4e(id,jd,kd),u5e(id,jd,kd)
      common /aero05/ u1m(id,jd,kd),u2m(id,jd,kd),u3m(id,jd,kd)
      common /aero06/ u4m(id,jd,kd),u5m(id,jd,kd)

      common /bcnd01/ pt1,tt1,pout
      common /bcnd02/ ptinlt,ttinlt,ptv,ttv
      common /bcnd03/ pinold(jd,jd),vxinld(jd,kd),vximld(jd,kd)
      common /bcnd04/ poutld(jd,kd),u1imld(jd,kd)
      common /bcnd05/ ptin(jd,kd)

      common /solv01/ rla (id,jd,kd),step(id,jd,kd)
      common /solv02/ flai(id,jd,kd),flaj(id,jd,kd),flak(id,jd,kd)
      common /solve4/ fblc(id,jd,kd),fblp(id,jd,kd)

      common /velo01/ vx (id,jd,kd),vy (id,jd,kd),vz (id,jd,kd)
      common /velo02/ qc (id,jd,kd)

      common /geom01/ aix(id,jd,kd),aiy(id,jd,kd),aiz(id,jd,kd)
      common /geom02/ ajx(id,jd,kd),ajy(id,jd,kd),ajz(id,jd,kd)
      common /geom03/ akx(id,jd,kd),aky(id,jd,kd),akz(id,jd,kd)
      common /geom04/ volc(id,jd,kd),volp(id,jd,kd)

      common /array1/ ar1(id,jd,kd),ar2(id,jd,kd), ar3(id,jd,kd)

      common /read01/ rfd,ft,rg,ermax,rfpin,rfpout,ptvbl
      common /read03/ pstart,vxstar,vystar,vzstar,u1star

      common /ther01/ cpzucv,cp,cv,rcp,rcv
      common /ther02/ rgzucv,rgzucp,cvzurg,cpzurg





