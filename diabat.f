c--------------------------------------------------
c  n   number of internuclear distance points
c  nch number of channels   
c  r   internuclear distance vector
c  v   diabatic potential array
c  -----------------------------------------------
c
c $Id: diabat.f,v 1.1 2014/02/15 18:57:13 lyonsdav Exp $
c
c SUMMARY
c
c REVISION HISTORY
c $Log: diabat.f,v $
c Revision 1.1  2014/02/15 18:57:13  lyonsdav
c Initial revision
c
c

      implicit none
      integer n,i,nch,ng,ns,j
      parameter (n=212,nch=2)
      double precision r(n),a(n),om,a2(n),h,rs,e1s,e2s,
     .   u11(n),u22(n),u12(n),pom,e12(n),e22(n),e1(n),e2(n)
      parameter (h=0.05d0)
      external gleg,csspln,splint
c
      om=0.0d0
      open(unit=8, file='dia.out')
      open(unit=9, file='intdatam')
c
c     read d/dr data and fit to spline
      call csspln (n,r,e1,e12,e2,e22,a,a2)
c
      ng=1000
      ns=dint((r(n)-r(1))/h)
c      print*, ns
      u12(n)=0.0d0
      u11(n)=e1(n)
      u22(n)=e2(n)
      rs=r(n)
      j=1
      do 30 i=ns,1,-1
          rs=rs-h
          call gleg(n,rs,h,r,a,a2,pom,ng)
          om=om+pom
c          write(*,*) i,rs,r(n-j),pom
          if(dabs(rs-r(n-j)).le.0.001d0) then
c	       print*,'in the loop'
              call splint(r,e1,e12,n,rs,e1s)
              call splint(r,e2,e22,n,rs,e2s)
              u11(n-j)=e1s*dcos(om)**2+e2s*dsin(om)**2
              u22(n-j)=e1s*dsin(om)**2+e2s*dcos(om)**2
              u12(n-j)=(e2s-e1s)*dsin(om)*dcos(om)
              write(9,*) r(n-j), om
              j=j+1
          endif
30    continue

      do 40 i=1,n-1
          write(8,*) r(i), u11(i), u22(i), u12(i)
40    continue
      write(8,*) r(n),e1(n),e2(n),u12(n)     
      end
c
       subroutine gleg(ns,rs,h,es,css,cs2,int,n)
       INTEGER n, i, ns
       DOUBLE PRECISION int,x,y,a,b,rs,h
       double precision csx,csy,isum
       double precision es(ns),css(ns),cs2(ns)
       DOUBLE PRECISION xi(1000),wi(1000)
       external gauleg,splint
       a = rs
       b = rs+h 
       call gauleg(-1.0d0,1.0d0,xi,wi,n)
       isum = 0.0d0
       do 10 i=1,n/2
       x=(b+a)/2 + (b-a)/2*xi(i)
       y=(b+a)/2 - (b-a)/2*xi(i)
          call splint(es,css,cs2,ns,x,csx)
          call splint(es,css,cs2,ns,y,csy)
       isum = (csx+csy)*wi(i) +isum
10     continue
       int = (b-a)/2.0d0*isum
c       write (6,*) 'int= ', int
       return
       END
C
      subroutine gauleg(x1,x2,x,w,n)
c      implicit real *16 (a-h, o-z)
          integer n,m,i,j
          double precision x1,x2,x(n),w(n),xm,xl,p1,p2
          double precision pi,z,pp,p3,z1
          parameter (eps=3.0d-14)
          parameter (pi=3.141592653589793238462643d0)

          m=(n+1)/2
          xm=0.5d0*(x2+x1)
          xl=0.5d0*(x2-x1)

          do 10 i=1,m
              z=dcos(pi*(i-0.25d0)/(n+0.5d0))
1             continue
              p1=1.0d0
              p2=0.0d0
              do 20 j=1,n
                  p3=p2
                  p2=p1
                  p1=((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
20            continue

              pp=n*(z*p1-p2)/(z*z-1.0d0)
              z1=z
              z=z1-p1/pp

              if (abs(z-z1) .gt. eps) goto 1
              x(i)=xm-xl*z
              x(n+1-i)=xm+xl*z
              w(i)=2.0d0*xl/((1.0d0-z*z)*pp*pp)
              w(n+1-i)=w(i)
10        continue
          return
      end

