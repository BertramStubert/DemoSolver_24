      common /intg01/ im,jm,km,nmax,nstep,nrfd
      common /intg21/ imh,jmh,kmh

      common /grid01/ x  (id,jd),y  (id,jd)
      common /grid21/ xh (idh,jdh),yh (idh,jdh)

      common /aero01/ u1 (id,jd),u2 (id,jd),u3 (id,jd)
      common /aero02/ u4 (id,jd),u5 (id,jd)
      common /aero03/ u1e(id,jd),u2e(id,jd),u3e(id,jd)
      common /aero04/ u4e(id,jd),u5e(id,jd)
      common /aero05/ u1m(id,jd),u2m(id,jd),u3m(id,jd)
      common /aero06/ u4m(id,jd),u5m(id,jd)

      common /aero21/ u1h (idh,jdh),u2h (idh,jdh),u3h (idh,jdh)
      common /aero22/ u4h (idh,jdh),u5h (idh,jdh)
      common /aero23/ u1eh(idh,jdh),u2eh(idh,jdh),u3eh(idh,jdh)
      common /aero24/ u4eh(idh,jdh),u5eh(idh,jdh)
      common /aero25/ u1mh(idh,jdh),u2mh(idh,jdh),u3mh(idh,jdh)
      common /aero26/ u4mh(idh,jdh),u5mh(idh,jdh)

      common /bcnd01/ pt1,tt1,pout
      common /bcnd02/ ptinlt,ttinlt,ptv,ttv
      common /bcnd03/ pinold(jdh)

      common /solv01/ rla (id,jd),step(id,jd)
      common /solv02/ flai(id,jd),flaj(id,jd)
      common /solve4/ fblc(id,jd),fblp(id,jd)

      common /solv21/ rlah (idh,jdh),steph(idh,jdh)
      common /solv22/ flaih(idh,jdh),flajh(idh,jdh)
      common /solv24/ fblch(idh,jdh),fblph(idh,jdh)

      common /velo01/ vx (id,jd),vy (id,jd),vz (id,jd)
      common /velo02/ p  (id,jd),qc (id,jd)

      common /velo21/ vxh(idh,jdh),vyh(idh,jdh),vzh(idh,jdh)
      common /velo22/ ph (idh,jdh),qch(idh,jdh)

      common /geom01/ aix(id,jd),aiy(id,jd)
      common /geom02/ ajx(id,jd),ajy(id,jd)
      common /geom04/ volc(id,jd),volp(id,jd)

      common /geom21/ aixh (idh,jdh),aiyh (idh,jdh)
      common /geom22/ ajxh (idh,jdh),ajyh (idh,jdh)
      common /geom24/ volch(idh,jdh),volph(idh,jdh)

      common /array1/ ar1 (id,jd),  ar2(id,jd)
      common /array2/ ar1h(idh,jdh),ar2h(idh,jdh)

      common /read01/ rfd,ft,rg,ermax,rfpin,rfpout,ptvbl
      common /read03/ pstart,vxstar,vystar,u1star

      common /ther01/ cpzucv,cp,cv,rcp,rcv
      common /ther02/ rgzucv,rgzucp,cvzurg,cpzurg





