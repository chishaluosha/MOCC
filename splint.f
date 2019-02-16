c P.C. Stancil 10-8-92
c Evaluates a function using spline array from spline.f
c From Numerical Recipes p.89 (1986)
c
c $Id: splint.f,v 1.1 2014/02/15 18:57:13 lyonsdav Exp $
c
c SUMMARY
c
c REVISION HISTORY
c $Log: splint.f,v $
c Revision 1.1  2014/02/15 18:57:13  lyonsdav
c Initial revision
c
c

      subroutine splint(xa,ya,y2a,n,x,y)
c
      integer n,klo,khi,k
      double precision x,y,a,b,h
      double precision xa(n),ya(n),y2a(n)
c
      klo=1
      khi=n
c
 10   if ((khi-klo) .gt. 1) then
         k=(khi+klo)/2
         if (xa(k) .gt. x) then
            khi=k
         else
            klo=k
         endif
         goto 10
      endif
c
      h=xa(khi)-xa(klo)
      if (h .eq. 0.0d0) pause 'Bad XA input'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+ 
     .  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*h**2/6.0d0
c
      return
      end 

