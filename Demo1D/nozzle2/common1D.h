      common /intg01/ im,nmax,nstep,nrfd
      common /grid01/ x  (id)

      common /aero01/ u1 (id),u2 (id),u3 (id)
      common /aero03/ u1e(id),u2e(id),u3e(id)
      common /aero05/ u1m(id),u2m(id),u3m(id)

      common /bcnd01/ pt1,tt1,pout
      common /bcnd02/ ptinlt,ttinlt,ptv,ttv
      common /bcnd03/ pinold,vxinld,vximld
      common /bcnd04/ poutld,u1imld

      common /solv01/ rla (id),step(id)
      common /solv02/ flai(id)
      common /solve4/ fblc(id),fblp(id)

      common /velo01/ vx (id)
      common /velo02/ p  (id),qc (id)

      common /geom01/ aix(id)

      common /geom04/ volc(id),volp(id)

      common /array1/ ar1(id)

      common /read01/ rfd,ft,rg,ermax,rfpin,rfpout,ptvbl
      common /read03/ pstart,vxstar,u1star

      common /ther01/ cpzucv,cp,cv,rcp,rcv
      common /ther02/ rgzucv,rgzucp,cvzurg,cpzurg





