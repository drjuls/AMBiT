 
      subroutine bsplvb(t,jhigh,index,x,left,biatx)
      implicit doubleprecision(a-h,o-z)
      parameter(JMAX=20)
      dimension biatx(jhigh),t(1),deltal(jmax),deltar(jmax)
      data j/1/
      go to (10,20),index
 10   j = 1
      biatx(1) = 1.0
      if(j .ge. jhigh)
     *   go to 99
 20   continue
         jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
         saved = 0.0
         do 26 i = 1,j
             term = biatx(i)/(deltar(i) + deltal(jp1-i))
             biatx(i) = saved + deltar(i)*term
             saved = deltal(jp1-i)*term
 26      continue
         biatx(jp1) = saved
         j = jp1
         if(j .lt. jhigh)
     *go to 20 
 99   RETURN
      end


      subroutine bsplvd(t,k,x,left,dbiatx,nderiv)
      implicit doubleprecision(a-h,o-z)
*******  changes from de Boor
      parameter(KX=15)
      dimension a(KX,KX),dbiatx(KX,KX)
*  also a missing in arg list
*****************
      dimension t(1)
 
      mhigh = max0(min0(nderiv,k),1)
      kp1 = k + 1
      call bsplvb(t,kp1-mhigh,1,x,left,dbiatx)
      if(mhigh .eq. 1) go to 99
 
      ideriv = mhigh
      do 15 m = 2,mhigh
          jp1mid = 1
          do 11 j = ideriv,k
             dbiatx(j,ideriv) = dbiatx(jp1mid,1)
             jp1mid = jp1mid + 1
 11       continue
          ideriv = ideriv - 1
          call bsplvb(t,kp1-ideriv,2,x,left,dbiatx)
 15   continue
 
      jlow = 1
      do 20 i = 1,k
         do 19 j = jlow,k
            a(j,i) = 0.0
 19      continue
         jlow = i
         a(i,i) = 1.0
 20   continue
 
      do 40 m = 2,mhigh
         kp1mm = kp1 - m
         fkp1mm = kp1mm
         il = left
         i = k
 
         do 25 ldummy = 1,kp1mm
            factor = fkp1mm/(t(il+kp1mm) - t(il))
            do 24 j = 1,i
               a(i,j) = (a(i,j) -a(i-1,j))*factor
 24         continue
            il = il - 1
            i = i - 1
 25      continue
 
         do 36 i = 1,k
              sum = 0.0
              jlow = max0(i,m)
              do 35 j = jlow,k
                 sum = sum + a(j,i)*dbiatx(j,m)
 35           continue
              dbiatx(i,m) = sum
 36      continue
 40   continue 
 99   RETURN
      end


      function bvalue(t,bcoef,n,k,x,jderiv)
      implicit doubleprecision(a-h,o-z)
      parameter(KMAX=20)
      dimension bcoef(n),t(1),aj(KMAX),dl(KMAX),dr(KMAX)
 
      bvalue = 0.0
      if(jderiv .ge. k)                go to 99
      call interv (t,n+k,x,i,mflag)
      if(mflag .ne. 0)                 go to 99
 
      km1 = k - 1
      if(km1 . gt. 0)                  go to 1
      bvalue = bcoef(i)
                                       go to 99
 
 1    jcmin = 1
      imk = i - k
      if(imk .ge. 0)                    go to 8
      jcmin = 1 - imk
 
      do 5 j = 1,i
         dl(j) = x - t(i+1-j)
 5    continue
      do 6 j = i,km1
         aj(k-j)  = 0.0
         dl(j) = dl(i)
 6    continue
                                      go to 10
 
 8    do 9 j = 1,km1
         dl(j) = x - t(i+1-j)
 9    continue
 
 10   jcmax = k
      nmi = n - i
      if(nmi .ge. 0)                  go to 18
      jcmax = k + nmi
      do 15 j = 1,jcmax
          dr(j) = t(i+j) - x
 15   continue
      do 16 j = jcmax,km1
          aj(j+1) = 0.0
          dr(j) = dr(jcmax)
 16   continue
                                     go to 20
 
 18   do 19 j = 1,km1
        dr(j) = t(i+j) - x
 19   continue
 
 20   do 21 jc = jcmin,jcmax
          aj(jc) = bcoef(imk + jc)
 21   continue
 
      if(jderiv .eq. 0)             go to 30
      do 23 j = 1,jderiv
         kmj = k - j
         fkmj = kmj
         ilo = kmj
         do 22 jj = 1,kmj
            aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
            ilo = ilo - 1
 22      continue
 23   continue
 
 30   if(jderiv .eq. km1)           go to 39
      do 33 j = jderiv+1,km1
         kmj = k - j
         ilo = kmj
         do 32 jj = 1,kmj
           aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
           ilo = ilo - 1
 32      continue
 33   continue
 39   bvalue = aj(1)
 99   RETURN
      end


      subroutine gauss(n,x,w)
      implicit doubleprecision(a-h,o-z)
******************************************************************
*
*  Gaussian coordinates and weights for the interval [0..1]
*        adapted from "setgau" in Numerical Recipes
*******************************************************************
      dimension x(n),w(n)
      data eps /1d-15/
      pi = dacos(-1.0d0)
      m = (n+1)/2
      do 120 i = 1,m
         z = cos(pi*(i-0.25d0)/(n+0.5d0))
 100     continue
            p1 = 1.d0
            p2 = 0.d0
            do 110 j = 1,n
               p3 = p2
               p2 = p1
               p1 = (( 2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
 110        continue
            pp = n*(z*p1-p2)/(z*z-1.d0)
            z1 = z
            z = z1-p1/pp
         if(dabs(z-z1).gt.eps) go to 100
         x(i) = 0.5d0*(1.0d0 - z)
         x(n+1-i) = 0.5d0*(1.0d0 + z)
         w(i) = 1.d0/((1.d0-z*z)*pp*pp)
         w(n+1-i) = w(i)
 120  continue
      return
      end


      subroutine interv (xt,lxt,x,left,mflag)
      implicit doubleprecision(a-h,o-z)
      dimension xt(lxt)
      data ilo /1/
 
      ihi = ilo + 1
      if(ihi .lt. lxt)                go to 20
      if(x .ge. xt(lxt))              go to 110
      if(lxt . le. 1)                 go to 90
      ilo = lxt - 1
      ihi = lxt
 
 20   if(x .ge. xt(ihi))              go to 40
      if(x .ge. xt(ilo))              go to 100
 
      istep = 1
 31   continue
         ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)              go to 35
         if (x .ge. xt(ilo))          go to 50
         istep = istep*2
      go to 31
 35   ilo = 1
      if(x .lt. xt(1))                go to 90
                                      go to 50
 
 40   istep = 1
 41   continue
         ilo = ihi
         ihi = ilo + istep
         if(ihi .ge. lxt)             go to 45
         if(x .lt. xt(ihi))           go to 50
         istep = istep*2
      go to 41
 
 45   if (x .ge. xt(lxt))             go to 110
      ihi = lxt
 50   continue
            middle = (ilo + ihi)/2
            if(middle .eq. ilo)       go to 100
            if(x .lt. xt(middle))     go to 53
            ilo = middle
          go to 50
 53       ihi = middle
      go to 50
 
 90   mflag = -1
      left = 1
      RETURN
 
 100  mflag = 0
      left = ilo
      RETURN
 
 110  mflag = 1
      left = lxt
      RETURN
      end
