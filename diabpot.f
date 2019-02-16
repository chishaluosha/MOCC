        program diabpot
c
c  to prepare the input data for diabatic potentials for MOCC.f   
c
c  read 1. adiabatic potentials;        file name: <pots>
c       2. radial coupling elements;    file name: <rads>
c       3. rotational coupling element; file name: <rots>
c  produce the diabatic potential       file name: <pot.diabatic>
c
c  modified by LBZ on Feb. 10, 2004, in Univ. of Georgia,
c  on the basis of the PCS's version 
c
c    npot: point number of potentials
c     nrr: point number of radial or rotational coupling elements
c          it is assumed that both have the same points. 
c
        implicit real*8(a-h,o-z)
        namelist/inpdiab/nch,npot,nrr,mrad,mrot,h
c
        open(3,file='input.diabatic',status='old')  
        open(4,file='pot.diabatic',status='unknown')  
        open(14,file='cci.mat',status='unknown')  
        read(3,nml=inpdiab)
        mrr=mrad+mrot
        call diabpot1(nch,npot,nrr,mrad,mrr,h)
c
        stop '(DIABPOT normal exit!)'
        end
c********************************************************************zlb
        subroutine diabpot1(n,npot,nrr,mrad,mrr,h)
        implicit double precision (a-h,o-z)
czlb        parameter(n=7,ndat=40,ndata=40)
        integer line(mrr),coul(mrr)
        dimension e(n,n),a(n,n),delta(n,n),vs(npot,n,n)
        common/dimen/n1,n2,n3
        double precision w(n,n),v(n,n),u(n,n),ui(n,n),
     1      uiu(n,n),lw(n,n),lv(n,n),lm(n,n),vl(nrr,n,n)
c
        dimension rp(npot),pot(npot,n),pot2d(npot,n)
        dimension rcp(nrr),coup(nrr,mrr),coup2d(nrr,mrr)
        dimension vsl(nrr,n,n)
c
        n1=n
        n2=n
        n3=n
        nch=n
        nc=nrr
czlb ------------------------------------------------       
c         call input1(ndat)
c         call input2(ndata)
        call inputpot(npot,nch,rp,pot,pot2d)
        call inputcp(nrr,mrad,mrr,rcp,coup,coup2d,line,coul)
c         r=30.000d0
c         h=0.00001d0
c         max=2800002
        r=rcp(nrr)
        max=dabs(rcp(1)-r)/h+2+h
czlb ------------------------------------------------        
        r=r+h
c
c intialize u 
c
        do 1 i=1,n
        do 1 j=1,n
        delta(i,j)=0.0d0
        if(i.eq.j) delta(i,j)=1.0d0 
1       u(i,j)=delta(i,j)
c
c start loop integrate backwards
c
        do 100 m=1,max
c
czlb        call adiab(r,e,a,lm,n,npot,nrr)
        call adiab(r,e,a,lm,nch,npot,nrr,mrad,mrr,rp,rcp,pot,pot2d,
     1             coup,coup2d,line,coul)
construct diabatic potential at point r-h
        r=r-h        
        do 2 i=1,n
        do 2 j=1,n
        sum=0.0d0
        do 21 k=1,n
c notice global phase of a = -1
21      sum=sum+(delta(i,k)-h*a(i,k))*u(k,j)
c21      sum=sum+u(i,k)*(delta(k,j)+h*a(k,j))
2       w(i,j)=sum
c reset u matrix
        do 22 i=1,n
        do 22 j=1,n
22      u(i,j)=w(i,j)       

        do 25 i=1,n
        do 25 j=1,n
25      ui(i,j)=delta(i,j)
c find u inverse
        call gaussl(w,ui,n)      

c construct diabatic term

        do 31 i=1,n
        do 31 j=1,n
        sum=0.0d0
        suml=0.0d0
        do 301 k=1,n
c301     sum=sum+e(i,k)*ui(k,j)
        sum=sum+e(i,k)*u(k,j)
301     suml=suml+lm(i,k)*u(k,j)
        w(i,j)=sum         
31      lw(i,j)=suml
c
        do 32 i=1,n
        do 32 j=1,n
        sum=0.0d0
        sumu=0.0d0
        suml=0.0d0
        do 302 k=1,n
        sumu=sumu+u(k,i)*u(k,j)
c302     sum=sum+u(i,k)*w(k,j)       
        sum=sum+ui(i,k)*w(k,j)
302     suml=suml+ui(i,k)*lw(k,j)
c new v matrix
        lv(i,j)=suml
        uiu(i,j)=sumu
32      v(i,j)=sum
                

c write out discrete grid
c          if (r .le. rcp(nc)) then
c          if (r .le. (rcp(nc)+h)) then
           ratio=rcp(nc)/r
           if(dabs(ratio-1.0d0).le.(h/10.0)) then
             write(14,120) nc, rcp(nc),r
             do 111 ii=1,n 
             vs(nc,ii,ii)=v(ii,ii)
                vsl(nc,ii,ii)= vs(nc,ii,ii)
                write(14,140) (uiu(ii,j), j=1,n)
                do 112 jj=ii+1,n
                vs(nc,jj,ii)=v(jj,ii)
                vl(nc,jj,ii)=lv(jj,ii)
czlb
                vsl(nc,jj,ii)= vs(nc,jj,ii)+vl(nc,jj,ii)
czlb
 112            continue
 111         continue
c
c             write(19,999) rcp(nc),v(1,1),v(2,2),v(3,3)
c     .                     ,lv(1,2),v(1,3),lv(2,3)
            nc=nc-1
           endif
303     continue     
100     continue
120     format(i5,2f20.10)
140     format(9f13.9)
c
        do i =1,nrr
          do j =1,nch
            do k=1,j
              vsl(i,k,j)= vsl(i,j,k)
            enddo
          enddo
        enddo 
c
        do i =1,nrr
          write(4,*) '-----------------------------',
     1               '--------------------------'
c                      1234567890123456789012345678901234567890
          write(4,400)'   Matrix of Diabatic Potentials at R = ',rcp(i)
          do j = 1,nch
            do k=1,j 
              if(i.eq.(nrr-1).and.abs(vsl(i,j,k)).lt.1.0d-20) 
     .        vsl(i,j,k) =1.0D-10
              if(i.eq.nrr.and.abs(vsl(i,j,k)).lt.1.d-20) 
     .          vsl(i,j,k)=1.0D-11
            enddo
            write(4,500) (vsl(i,j,k),k=1,j)
          enddo
c
          do kz=1,nch
          do jz=kz,nch
           lchn=jz*10+kz
           write(lchn,*) rcp(i),vsl(i,jz,kz)
          enddo
          enddo
c
        enddo 
c
400     format(a40,1pE16.8)
500     format(10(1pD20.12))
        return
        end               
c****************************************************************************zlb
        subroutine adiab(x,e,a,lv,nch,npot,nrr,mrad,mrr,rp,rcp,
     1                   pot,pot2d,coup,coup2d,line,coul)
        implicit real*8(a-h,o-z)
        integer line(mrr),coul(mrr)
        double precision e(nch,nch),a(nch,nch),lv(nch,nch)
        dimension rp(npot),pot(npot,nch),pot2d(npot,nch)
        dimension rcp(nrr),coup(nrr,mrr),coup2d(nrr,mrr)
        dimension pott(npot),pott2d(npot)
c
c initialize potential matrix
c
        do 100 i=1,nch
        do 100 j=1,nch
100     e(i,j)=0.0d0
c        
        do j=1,nch       
          do i=1,npot
            pott(i)=pot(i,j)
            pott2d(i)=pot2d(i,j)
          enddo
          call splint(rp,pott,pott2d,npot,x,e(j,j))
        enddo
c
c initialize a matrix
c
        do 200 i=1,nch
        do 200 j=1,nch
        lv(i,j)=0.0d0
200     a(i,j)=0.0d0
c
        do j=1,mrad      
          do i=1,npot
            pott(i)=coup(i,j)
            pott2d(i)=coup2d(i,j)
          enddo
          call splint(rcp,pott,pott2d,nrr,x,zrcp)
          a(line(j),coul(j))=zrcp
        enddo
c
        do j=mrad+1,mrr      
          do i=1,npot
            pott(i)=coup(i,j)
            pott2d(i)=coup2d(i,j)
          enddo
          call splint(rcp,pott,pott2d,nrr,x,zrcp)
          lv(line(j),coul(j))=zrcp
        enddo
c
c make matrix symmetric
c
        do 300 j=1,nch
        do 300 i=1,j
        lv(i,j)=lv(j,i)
        a(i,j)=-a(j,i)
c        if (x.eq.20.0d0) print*, i,j, a(i,j),a(j,i),lv(i,j),lv(j,i)
300     continue
        return
        end
c****************************************************************************zlb
        subroutine inputcp(nrr,mrad,mrr,rcp,pot,pot2d,line,coul)
        implicit real*8(a-h,o-z) 
        integer line(mrr),coul(mrr)
        dimension pott(nrr),ypott(nrr)
        dimension rcp(nrr),pot(nrr,mrr),pot2d(nrr,mrr)
c
c input reads in x, y values spline nodes
c
c     nrr: point number of radial or rotational coupling elements,
c          it is assumed that both have the same points. 
c    mrad: number of nonzero radial coupling 
c    mrot: number of nonzero rotational coupling 
c     mrr: mrad+mrot
c    
        open(unit=8, file='rads')
        open(unit=9, file='rots')
c
        a=0.0d0
        b=0.0d0
        do 100 j=1,mrad
            read(8,*) line(j),coul(j)
          do 20 i=1,nrr
            read(8,*) rcp(i),pot(i,j)
            pott(i)=pot(i,j)
20        continue
          call spline(rcp,pott,nrr,a,b,ypott)
          do 30 i=1,nrr
            pot2d(i,j)=ypott(i)
30        continue
100     continue
c
        do 200 j=mrad+1,mrr
            read(9,*) line(j),coul(j)
          do 40 i=1,nrr
            read(9,*) rcp(i),pot(i,j)
            pott(i)=pot(i,j)
40        continue
          call spline(rcp,pott,nrr,a,b,ypott)
          do 50 i=1,nrr
            pot2d(i,j)=ypott(i)
50        continue
200     continue
        return
        end
c****************************************************************************zlb
        subroutine inputpot(npot,nch,rp,pot,pot2d)
        implicit real*8 (a-h,o-z)
c
c input reads in x, y values spline nodes
c
c       purpose -- read adiabatic potentials and calculate the 2nd derivative.  
c
c     nch: number of channels 
c    npot: point number of potentials
c    pott: 1st derivative of adiabatic potentials
c   pot2d: 2nd derivative of adiabatic potentials 
c
        dimension pott(npot),ypott(npot)
        dimension rp(npot),pot(npot,nch),pot2d(npot,nch)
c
        open(unit=7,file='pots')
c
        yp1=0.0d0
        yp2=0.0d0
        do 100 j=1,nch
          do 20 i=1,npot
            read(7,*) rp(i),pot(i,j)
            pott(i)=pot(i,j)
20        continue
          call spline(rp,pott,npot,yp1,yp2,ypott)
          do 30 i=1,npot
            pot2d(i,j)=ypott(i)
30        continue
c
100     continue
        return
        end
c****************************************************************************zlb
      subroutine spline(x,y,n,yp1,ypn,y2)
      parameter (nmax=100002)
      implicit real*8(a-h,o-z)
      double precision x(n),y(n),y2(n),u(nmax)
      if(yp1.gt. .99d30) then
      y2(1)=0.0d0
      u(1)=0.0d0
      else
      y2(1)=-0.5d0
      u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.d0
      y2(i)=(sig-1.d0)/p
      u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn .gt. .99d30) then
      qn=0.d0
      un=0.d0
      else
      qn=0.5d0
      un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end
c****************************************************************************zlb
      subroutine splint(xa,ya,y2a,n,x,y)
      implicit real*8(a-h,o-z)
      dimension xa(n),ya(n),y2a(n)
      klo=1
      khi=n
1     if (khi-klo .gt. 1) then
      k=(khi+klo)/2
      if (xa(k).gt.x) then
        khi=k
      else
        klo=k
      endif
      go to 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) stop
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     x  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
c (disregard b.z.)  slpoe=(ya(khi)-ya(klo))/h
c        y=ya(klo)+slpoe*(x-xa(klo))
      return
      end
c****************************************************************************zlb
      SUBROUTINE GAUSSL (A,B,N)                                        
C*****SOLVES N BY N MATRIX EQUATION AX=B BY GAUSSIAN ELIMINATION.
C*****ANSWER IS RETURNED IN B, A IS DESTROYED.
      implicit real*8(a-h,o-z)
      common/dimen/n1,n2,n3
      DIMENSION A(n,n),B(n,n)                                          
      IF(N.EQ.1) GO TO 280                                             
      NS = N - 1                                                       
      DO 200 I=1,NS                                                    
      IG = I + 1                                                       
      AT = A(I,I)                                                      
      DO 120 J=IG,N                                                    
120   A(I,J) = A(I,J)/AT                                               
      DO 130 J=1,N                                                     
130   B(I,J) = B(I,J)/AT                                                
      DO 150 K=1,N                                                     
      IF(K-I) 151,150,151                                              
  151 AT = A(K,I)                                                      
      DO 140 L=IG,N                                                     
140   A(K,L) = A(K,L) - A(I,L)*AT                                      
      DO 145 L=1,N                                                     
145   B(K,L) = B(K,L) - B(I,L) *AT                                     
150   CONTINUE                                                         
200   CONTINUE                                                         
      AT = A(N,N)                                                      
      DO 230 J=1,N                                                     
230   B(N,J) = B(N,J)/AT                                              
      DO 250 K=1,NS                                                    
      DO 250 L=1,N                                                     
      AT = A(K,N)                                                      
250   B(K,L) = B(K,L) - B(N,L)*AT                                      
      RETURN                                                           
  280 B(1,1)= B(1,1)/A(1,1)                                            
      A(1,1)= 1.0D0                                                    
      RETURN                                                           
      END
c****************************************************************************zlb