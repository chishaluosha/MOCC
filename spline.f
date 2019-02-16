c P.C. Stancil 10-8-92
c Determines parameter array for spline
c From Numerical Recipes p.88 (1986)
c
c $Id: spline.f,v 1.1 2014/02/15 18:57:13 lyonsdav Exp $
c
c SUMMARY
c
c REVISION HISTORY
c $Log: spline.f,v $
c Revision 1.1  2014/02/15 18:57:13  lyonsdav
c Initial revision
c
c

c
      subroutine spline(x,y,n,yp1,ypn,y2)
c
      parameter (nmax=300)
      integer i,n,k
      double precision yp1,ypn,sig,p,qn,un
      double precision x(n),y(n),y2(n),u(nmax)
c
      if (yp1 .gt. 0.99d30) then
         y2(1)=0.0d0
         u(1)=0.0d0
      else
         y2(1)=-0.5d0
         u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
c
      do 10 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2
         y2(i)=(sig-1.0d0)/p
         u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     .         /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
 10   continue
c
      if (ypn .gt. 0.99d30) then
         qn=0.0d0
         un=0.0d0
      else
         qn=0.5d0
         un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
c
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
c
      do 20 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
 20   continue
c
      return
      end      

