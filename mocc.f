      program mocc   
c
c modified by LBZ on Feb 11, 2004, in Univ. of Georgia,
c based on PCS's version
c     
c    nch: number of channels 
c   ndat: number of potential data points per channel
c   npar: number of partial waves per pass
c
      open(unit=3, file='input.capture')
      open(unit=16, file='pot.diabatic',status='old')
      read(3,*) nch,ndat,npar
      write(*,*) 'Point 1'
      npot=nch*(nch+1)/2
      call mocc1(nch,ndat,npar,npot)
c
      stop '(MOCC normal exit!)'
      end
c********************************************************************ZLB
c---------------------------------------------------------- 
c 12-17-03 modified for general case, pcs
c 2-9-95   modified for SiH2+, two channels, pcs
c 11-5-94  modified for arbitrary number of channels, pcs
c 10-24-94 modified for singlet sigma+ manifold of
c          SiHe4+, three channels,  pcs
c
c Driver for n-channel electron capture problem. Adopted
c fom B. Zygelman's code for the 4-channel singlet sigma
c NH4+ calculation. See Makefile for list of required
c subroutines. 
c
c User is required to change nch, ndat, and/or npar
c to suit their problem in all subroutines. These indice
c limits are given in the first parameter statements near
c the top of each subroutine. 
c
c Short range fits and possibly long range fits may need 
c modification in dpot.f to suit the given problem. See
c comments in dpot.f for instructions. Text headers in outdat.f
c for output of partial cross sections (fort.4) and the 
c s- and k-matrics (fort.2) may also need changing.
c No changes are required for spline.f, splint.f, gaussv.f,
c coulfg.f or jwkb.f.
c
c Required input is from input.capture, fort.17,
c and fort.16. See example
c input files and input structure given below for input.capture
c and in input.f for fort.17 and fort.16
c 
c Total and state-selective cross sections are given in 10^-16 cm2
c in fort.30 and fort.20-29
c -----------------------------------------------------------
c nch    number of channels 
c ndat   number of potential data points per channel
c npar   number of partial waves per pass
c npass  number of passes for cycles of npar partial waves
c emu    reduced mass in amu
c rmu    reduced mass in au
c rstart minimum internuclear distance
c rn     maximum internuclear distance
c jmin   minimum partial angular momentum 
c jstep  step size of partial waves
c ne     number of relative energy points to calculate
c h      internuclear distance integration step size
c ev     array of relative energies in eV
c e      relative energy in au 
c delta  absolute energy in au of incoming channel 
c        at infinite separation
c rs     minimum internuclear distance of diabatic
c        potential data, start of short range fit
c rl     maximum internuclear distance of diabatic
c        potential data, start of long range fit
c z      charge of first ion for any channel, must be positive 
c z2     charge of second particule, zero if one specie
c        is neutral
c alpha  static dipole polarizability of neutral specie in au
c ea,ei  array of infinite separation energy defects in au
c        of outgoing channels with respect to lowest channel (1)
c        ea must be ab initio energies, ei can be either
c        ab initio or experimental energies
c gamma
c parsig state-selective partial cross section 
c sig    total state-selective cross section
c sigt   total electron capture cross section
c jv     partial wave angular momentum array
c smat   s-matrix
c fkmat  k-matrix
c
c--------------------------------------------------------------
      subroutine mocc1(nch,ndat,npar,npot)
      implicit real*8 (a-h,o-z)
      integer nch, ndat,npar
     . ,npass,jstep, nprint,ne
     . ,i,j,k,ii,kcount,impot,npot,ici 
czlb      parameter(nch=3,ndat=31,npar=64)
czlb      parameter(npot=nch*(nch+1)/2)
      dimension iclose(nch)
      integer coupindx
      double precision rs,rl,alpha,cinf  
      double precision yn(npar,nch,nch),rjl(npar,nch,nch)
     . ,rnl(npar,nch,nch),rjlp(npar,nch,nch)
     . ,rnlp(npar,nch,nch),rki(nch),eta(nch),gamma(nch)
      double precision ei(nch),ke(nch),sig(nch,nch),sigt(nch)
      double precision ev(100),he,erelk(nch),z2
      double precision z,ea,delta,twomu,scale,e,g(nch),
     . erel,cent,unit,emu,rmu,h,rstart,rn,estart,lam,jin1
      double precision jv(npar),jmin,jin,pi,fkmat,smat,parsig
czlb
      dimension erel(nch),cent(npar,nch)
      dimension unit(nch,nch)
      dimension fkmat(npar,nch,nch),smat(npar,nch,nch),
     1          parsig(npar,nch,nch)
      dimension ea(nch),cinf(npot),impot(nch)
      dimension alpha(nch),z(nch),z2(nch)
      dimension lam(nch),coupindx(nch,nch)
      dimension pot(ndat,npot),ypot(ndat,npot)
      dimension r(ndat)
c
      common/pot/twomu,e,scale
      common/count/lcount
      common/fit/rs,rl,nprint
czlb        06/27/04
      dimension xi(nch),r2av(nch)
czlb 
c
c      common/pot/twomu,e,erel(nch),cent(npar,nch),scale
c      common/une/unit(nch,nch)
c      common/fmtrix/ fkmat(npar,nch,nch),smat(npar,nch,nch)
c     . ,parsig(npar,nch,nch)
c      common/count/lcount
c      common/fit/rs,rl,nprint
c      common/fit1/ea(nch),delta,cinf(npot),impot(nch)
c      common/fit2/alpha(nch),z(nch),z2(nch)
c      common/lamm/lam(nch),coupindx(nch,nch)
czlb
c     set all infinity coupling constants to zero
c
      do 27 i=1,npot
         cinf(i)=0.0d0
27    continue 
c     read in individual run data from input.capture
c
czlb      open(unit=3, file='input.capture')
           write(*,*) 'Point 3'
           read(3,*) rstart,jmin,npass,jstep,nprint
           write(*,*) 'Point 3'
           read(3,*) ne,h,rn,emu
           read(3,*) estart, he
           read(3,*) delta, rs, rl
czlb           read(3,*) (ea(i),alpha(i),z(i),z2(i),lam(i),
czlb     .                     g(i),impot(i),i=1,nch)
           read(3,*) (ea(i),alpha(i),xi(i),r2av(i),z(i),z2(i),lam(i),
     .                     g(i),impot(i),i=1,nch)
           read(3,*) (ei(i),i=1,nch)
           read(3,*) ici
           read(3,*) (j,cinf(j),i=1,ici)
c
        nofg=npar*npass+101
        xmax=nofg-1 
c
c determine gamma and coupling types
        do 21 i=1,nch
           gamma(i)=z(i)*z2(i)
           do 21 j=1,nch
           coupindx(i,j)=0
21      continue
        write(6,998)
        write(6,997)
        write(6,996)
        write(6,995)
        write(6,994)
        do 22 i=1,nch-1
           do 22 j=i+1,nch
             coupindx(i,j)=dabs(dint(lam(i)-lam(j)))+1
             coupindx(j,i)=coupindx(i,j)
             print*,i,j,coupindx(i,j)
22      continue
c
c     read in diabatic potential array from fort.17 
c     and fort.16 and fit to spline
c
      call input(r,pot,ypot,alpha,z,z2,nch,ndat,npar,npot,xi,r2av)
c
c write headers for cross section output files 
c
        do 885 i=1,nch
           write(19+i,883)
           write(19+i,879)
           write(19+i,884) (i,j,j=1,nch)
           write(39+i,878)
           write(39+i,877)
           write(39+i,876) (i,j,j=1,nch)
885     continue
883     format('# state-to-state electron capture cross sections') 
879     format('#    E              sigma_ij (10^-16 cm2)')
c *******
c **************************** fix 3 to be nch *************
c *******
czlb884     format('#  (eV/u)',8x,3(i2,i2,10x))
884     format('#  (eV/u)',8x,10(i2,i2,10x))
        write(30,882)
        write(30,881) (i,i=1,nch)
882     format('# Total electron capture cross sections (10^-16 cm2)') 
c *******
c **************************** fix 3 to be nch *************
c *******
czlb881     format('#  ', 3('eV/u',9x,i2,13x))
881     format('#  ', 10('eV/u',9x,i2,13x))
878     format('# Partial cross sections (a_0^2/g)')
877     format('#    J    sigma_ij(J)')
c *******
c **************************** fix 3 to be nch *************
c *******
czlb876     format('#',8x,3(i2,i2,10x))
876     format('#',8x,10(i2,i2,10x))

c loop over relative energy
c
        do 556 kount=1,ne
           write(*,*) 'E point =', kount
           ev(kount)=estart+dfloat(kount-1)*he
           e=ev(kount)/27.2113962d0


C.....Using At. Dat. Nucl. Dat. Tables 18, 587 (1976) gives 1822.88732

      rmu=1822.88732D0*emu

      pi=4.0d0*datan(1.0d0)
c READ IN GAMMMA VALUES
      DO 130 I=1,NCH
      EREL(I)=E-EI(I)
      ke(i)=erel(i)*27.2113962d0/emu
      erelk(i)=erel(i)
c check for closed channels
       iclose(i)=0
       if(erel(i).lt.0.d0) then
c       print*, i
       iclose(i)=1
       erelk(i)=dabs(erel(i))
       end if
      RKI(I)=DSQRT(2.D0*RMU*(ERELk(I)))
  130 ETA(I)=RMU*GAMMA(I)/RKI(I)
      TWOMU=RMU+RMU

       scale=h*h/6.0d0
      do 201 i=1,nch
         sigt(i)=0.0d0
      do 201 j=1,nch
         sig(i,j)=0.0d0
      unit(i,j)=0.0d0
201   unit(j,j)=1.0d0

c counters are initialized again

c initialize partial wave angular momenta
c array cent will be loaded in common block pot.

      jin=jmin

      do 555 i=1,npass
c      print*, (rki(ii), erel(ii),erelk(ii),ii=1,nch)
      call linita(jin,jstep,jv,erel,cent,nch,ndat,npar,npot)
c get new jmin
      jin1=jin
      jin=jv(npar)+1.0d0
c get asymptotic wavefunctions

      call asymp(rn,rki,eta,ei,jv,rjl,rjlp,rnl,rnlp,
     1           unit,nch,ndat,npar,npot,xmax,nofg)
c      ,nch,ndat,npar,npot)
      
c run the log-derivative integrator

      call ylog(h,rstart,rn,yn,erel,cent,unit,lam,coupindx,
     1          r,pot,ypot,alpha,z,z2,
     2          ea,delta,cinf,impot,nch,ndat,npar,npot,
     3          xi,r2av,jin1)
      
c evaluate k-matrix,s-matrix,cross-sections.

      call kmatrix(yn,rjl,rjlp,rnl,rnlp,jv,rki,nprint,
     1             unit,fkmat,smat,parsig,nch,ndat,npar,npot)

c call outputting routine, but first sum-up all
c     partial cross sections, multiply by degeneracy
c     factor, and convert to 10^-16 cm2
c
c        print*,'enter writeout
          do i99=1,npar
          do j99=1,nch
          do k99=1,nch
            if(parsig(i99,j99,k99).lt.1.0d-99) parsig(i99,j99,k99)=0.0d0
          enddo
          enddo
          enddo
c
          do 403 k=1,npar
           do 400 ii=1,nch
               do 400 jj=1,nch
                   sig(ii,jj)=sig(ii,jj)
     .             +g(ii)*parsig(k,ii,jj)*dfloat(jstep)
     .             *0.280028561d0
c
400        continue
c
c  write out J-partial cross sections
c   
            if (nprint.ge.1) then
               do 404 ii=1,nch
              write(39+ii,999) (k-1)*jstep+(i-1)*npar*jstep,
     .                    (parsig(k,ii,jj),jj=1,nch) 
c     .                    parsig(k,1,4),parsig(k,1,5),
c     .                    parsig(k,2,3),parsig(k,2,4)
404            continue
            end if
403       continue
555    continue
c
c print state-selective cross sections including elastic
c
       do 554 i=1,nch
          write(19+i,455) ke(i), (sig(i,j),j=1,nch)
             do 553 j=1,nch
                if(j.eq.i) goto 553
c sum-up total inelastic cross section
                sigt(i)=sigt(i)+sig(i,j)
553          continue
554    continue
c
c  print total inelastic cross sections   
c
       write(30,455) (ke(i),sigt(i),i=1,nch)
c *********
c ********* special case write for sum inelastic cross sections
c ******
       write(31,455) ke(3),sig(3,1)+sig(2,1)
czlb455    format(20g14.6)
czlb999    format(i6,10g14.6)
455    format(20(1pe14.7))
999    format(i6,10(1pe14.7))
998    format('**** Test of coupling types ****')
997    format('  indx=1, radial coupling')
996    format('  indx=2, rotational coupling')
995    format('  indx>2, other, coupling set to zero')
994    format(' i j indx')

556    continue
       return
       end
c********************************************************************ZLB
c 12-16-03 generalized, pcs
c 11-5-94 modified by pcs
c         subroutine to input diabatic potential matrix
c         from fort.17. Note required serial structure
c         of fort.17. Structure is diagonal elements first
c         (1,1 2,2 3,3 ...) followed by off-diagonal elements
c         (1,2 1,3 ... 2,3 ...). Subroutine also determines
c         spline fit.
c-----------------------------------------------------------
c ndat   number of potential data points per channel
c npot   size of potential array, includes diagonal
c        elements and upper triangle
c pot    array of diabatic potential data
c r      vector of internuclear distances
c ypot   array of spline data
c---------------------------------------------------
        subroutine input(r,pot,ypot,alpha,z,z2,nch,ndat,npar,npot,
     .                   xi,r2av)
        implicit real*8 (a-h,o-z)
        integer ndat, npot,nch,zc,zn
czlb        parameter(ndat=31, nch=3) 
czlb        parameter(npot=nch*(nch+1)/2)
c input reads in x, y values spline nodes
czlb
        dimension pot(ndat,npot),ypot(ndat,npot)
        dimension r(ndat)
        dimension alpha(nch),z(nch),z2(nch)
        dimension potm(nch,nch)
czlb        06/27/04
      dimension xi(nch),r2av(nch)
czlb 
c
c        common/dat/pot(ndat,npot),ypot(ndat,npot)
c        common/grid/r(ndat)
c        common/fit2/alpha(nch),z(nch),z2(nch)
czlb
        double precision pott(ndat),ypott(ndat),
     .     yp1(npot),yp2(npot),z,alpha,z2
     .    ,r,pot,ypot,ypp2,ypp1,dumb,tmp 
c
c set all first and last point derivatives to zero
c 
        do 223 j=1,npot
           yp1(j)=0.0d0
           yp2(j)=0.0d0
223     continue
c
c **** read potential (fort.17) and coupling (fort.16) 
c **** data, may require modification
czlb
czlb read potentials from file <pot.diabatic>
czlb
        do 23 i=1,ndat
          read(16,*)
          read(16,100) r(i)
          do j=1,nch
            read(16,*) (potm(j,k),k=1,j)
          enddo
czlb
czlb make matrix symmatric     
czlb
          do j=1,nch
            do k=1,j
              potm(k,j)=potm(j,k)
            enddo
          enddo
c
          do j=1,nch
            pot(i,j)=potm(j,j)        
          enddo
c
          jmin=nch+1
          jmax=nch+(nch-1)
          do j=1,nch
            do k=j+1,nch
              pot(i,jmin+k-j-1)=potm(j,k)
c              write(*,*) jmin+k-j-1
            enddo
            jmin0=jmin
            jmax0=jmax
            jmin=jmax0+1
            jmax=jmax0+(jmax0-jmin0)
          enddo
c
czlb           read(17,*)r(i),(pot(i,j),j=1,nch)
czlb           read(16,*)r(i),(pot(i,j),j=1+nch,npot)
c        write(51,*) r(i), pot(i,1),pot(i,2),pot(i,3)
23      continue
100     format(40x,E16.8)
c
c determine first point derivatives for all potentials
c and couplings and last point for couplings
c
        do 24 j=1,npot
           yp1(j)=(pot(2,j)-pot(1,j))/(r(2)-r(1))
           yp2(j)=(pot(ndat,j)-pot(ndat-1,j))
     .            /(r(ndat)-r(ndat-1))
24      continue
c
c determine last point derivatives for potentials only
c using asymptotic forms
c 
        write(6,999)
        zc=0
        do 53 j=1,nch
c    check if Coulomb channel
           if (z2(j).ne.0.0d0) then
              yp2(j)=-(z(j)*z2(j))/r(ndat)**2
              zc=zc+1
c              print*,'j=',j,'zc=',zc
              goto 53
           else
              zn=j
              print*,'  ion-neutral channel index =', zn 
           end if
c    for ion-neutral polarization potential 
czlb           yp2(j)=2.0d0*alpha(j)*z(j)**2/r(ndat)**5
           yp2(j)=2.0d0*alpha(j)*z(j)**2/r(ndat)**5
     .           -3.0*xi(j)*r2av(j)/r(ndat)**4  
c           print*,'dev',alpha(j),z(j),yp2(j)
 53     continue
c    number of  ion-ion Coulomb channels 
        print*,'  number of Coulomb channels =',zc       

c223     continue
c
c get spline coefficients
c for pot1
        do 52 j=1,npot
        ypp1=yp1(j)
        ypp2=yp2(j)
        do 50 i=1,ndat
50      pott(i)=pot(i,j)
        call spline(r,pott,ndat,ypp1,ypp2,ypott)
        do 51 i=1,ndat
51      ypot(i,j)=ypott(i)
52      continue
999     format('**** Test of channel types ****')
        return
c pot(i,j), ypot(i,j) contains the spline parameters at node i
c potential j.
        end
c********************************************************************ZLB
c 11-5-94 modified by pcs
c linita initializes channel angular momenta
c for "type-a" parity channels.
c shared arrays
c--------------------------------------------
c nch    number of channels
c npar   number of partial waves per pass
c--------------------------------------------
      subroutine linita(jmin,jstep,jv,erel,cent,nch,ndat,npar,npot)
      implicit real*8 (a-h,o-z)
      integer nch,npar,jstep
czlb      parameter(nch=3,npar=64)
      double precision twomu,e,erel,cent,scale,step,
     . l(npar),jv(npar),jmin,jval
czlb
      dimension erel(nch),cent(npar,nch)
      common/pot/twomu,e,scale
c      common/pot/twomu,e,erel(nch),cent(npar,nch),scale
czlb
c
       jval=jmin
       step=dfloat(jstep)
       do 10 k=1,npar
       jv(k)=jval
10     jval=jval+jstep
       do 20 i=1,nch
       do 20 k=1,npar
          l(k)=jv(k)
20     cent(k,i)=l(k)*(l(k)+1.0d0)
       return
       end
c********************************************************************ZLB
c  2-9-95  modified for SiH2+, two channels, pcs
c 11-5-94  modified for arbitrary number of channels, pcs
c          for high energies increase xmax
c----------------------------------------------------------
c nch    number of channels
c npar   number of partial waves per pass
c---------------------------------------------------------      
      subroutine asymp(rn,rki,eta,ei,jv,rjl,rjlp,rnl,rnlp,
     1                 unit,nch,ndat,npar,npot,xmax,nofg)
      implicit real*8 (a-h,o-z)
      integer nch,npar,l,ln,kfn,mode1,ifail
      double precision rho,deta,xmin,xmax
czlb      parameter(nch=3,npar=64,xmax=2000.0d0)
      double precision jv(1),unit,rn,a,b,t,rtki
c shared arrays
czlb
         dimension unit(nch,nch)
c         common/une/unit(nch,nch)
czlb
      double precision rki(nch),eta(nch),ei(nch),rjl(npar,nch,nch)
     1,rjlp(npar,nch,nch),rnl(npar,nch,nch),rnlp(npar,nch,nch)
c local arrays
       dimension l(nch)
czlb       double precision fc(3001,nch),gc(3001,nch),fcp(3001,nch),
czlb     .   gcp(3001,nch),f1(3001),f2(3001),g1(3001),g2(3001)
       double precision fc(nofg,nch),gc(nofg,nch),fcp(nofg,nch),
     .   gcp(nofg,nch),f1(nofg),f2(nofg),g1(nofg),g2(nofg)

c asymp return regular and irregular asymptotic function matrices.

         xmin=0.d0
         do 112 i=1,nch
         kfn=0
         rho=rn*rki(i)
         mode1=1

         if(eta(i).eq.0.0d0)  then
         kfn=1
         t=rho
         a=1.0d0
         b=-1.0d0
         else
         kfn=0
         t=1.0d0
         b=1.0d0
         a=0.d0
         end if

         deta=eta(i)
 
         call coulfg(rho,deta,xmin,xmax,f1,g1,f2,g2,mode1,kfn,ifail)
         maxl=dint(xmax)

         do 111 k=1,maxl-1
         fc(k,i)=f1(k)*t
         fcp(k,i)=a*f1(k)+t*f2(k)
         gc(k,i)=g1(k)*t*b
         gcp(k,i)=b*(a*g1(k)+t*g2(k))
111      continue
112      continue

30       continue
c construct zero-order asymptotic functions
         do 312 k=1,npar
         do 313 i=1,nch
313      l(i)=idint(jv(k))

         do 412 i=1,nch
         do 412 j=1,nch
         rtki=dsqrt(rki(j))
         ln=l(j)+1
         rjl(k,i,j)=unit(i,j)*fc(ln,j)/rtki
         rjlp(k,i,j)=unit(i,j)*fcp(ln,j)*rtki
         rnl(k,i,j)=unit(i,j)*gc(ln,j)/rtki
         rnlp(k,i,j)=unit(i,j)*gcp(ln,j)*rtki

412      continue
312      continue

c         do k=1,npar
c         write(15,*)rjl(k,1,1),rjlp(k,1,1),rnl(k,1,1),rnlp(k,1,1)
c         end do

         return
         end
c********************************************************************ZLB
c11-5-94  modified by pcs
c ylog is the vectorized log-deritative integrator,
c RSTART is starting value
c YN is the NxN log-deritative at RN. H is the step size.
c local arrays
c--------------------------------------------------------
c nch    number of channels
c npar   number of partial waves per pass
c-------------------------------------------------------
      subroutine ylog(h,rstart,rn,yn,erel,cent,unit,lam,coupindx,
     1                r,pot,ypot,alpha,z,z2, 
     2                ea,delta,cinf,impot,nch,ndat,npar,npot,
     3                xi,r2av,jstart)
      implicit real*8 (a-h,o-z)
      integer nch,npar,nol,n4,n1,m,max
czlb      parameter(nch=3,npar=64,ndat=31)
czlb      parameter(npot=nch*(nch+1)/2)
      double precision work1(npar,nch,nch),work2(npar,nch,nch)
     1,work3(npar,nch,nch),zne(npar,nch,nch),zno(npar,nch,nch),sum(npar)
c shared arrays
      double precision unit,re,ro,rn,rstart,h
      double precision v(npar,nch,nch),yn(npar,nch,nch)
      double precision lam,jstart
      integer coupindx 
      dimension lam(nch),coupindx(nch,nch)
      dimension pot(ndat,npot),ypot(ndat,npot)
      dimension r(ndat)
      dimension alpha(nch),z(nch),z2(nch)
      dimension ea(nch),cinf(npot),impot(nch)
czlb        06/27/04
      dimension xi(nch),r2av(nch)
czlb 
      dimension unit(nch,nch)
      dimension erel(nch),cent(npar,nch)
c      common/une/unit(nch,nch)
czlb      data nol/npar/,n4/npar/
      nol = npar
      n4  = npar
czlb

       n1=nch
       do 150 i=1,nch
       do 150 j=1,nch
       do 149 k=1,npar
149    zne(k,i,j)=unit(i,j)*1.0d20
150    continue
       re=rstart
       max=idnint((rn-rstart)/(2.0d0*h))
       do 200 m=1,max
       ro=re+h
       call potval(ro,h,v,erel,cent,unit,lam,coupindx,
     1             r,pot,ypot,alpha,z,z2,
     2             ea,delta,cinf,impot,nch,ndat,npar,npot,
     3             xi,r2av,jstart)
       do 160 i=1,nch
       do 160 j=1,nch
       do 159 k=1,npar
       work1(k,i,j)=unit(i,j)
       work2(k,i,j)=unit(i,j)
       work3(k,i,j)=unit(i,j)
       zno(k,i,j)=unit(i,j)+zne(k,i,j)
159    work3(k,i,j)=unit(i,j)+v(k,i,j)
160    continue
       call gaussv(zno,work1,sum,nch,nol,n1,n4)
       call gaussv(work3,work2,sum,nch,nol,n1,n4)
       do 171 i=1,nch
       do 171 j=1,nch
       do 168 k=1,npar
168    work3(k,i,j)=0.0d0
       do 170 n=1,nch
       do 169 k=1,npar
169    work3(k,i,j)=work3(k,i,j)+zne(k,i,n)*work1(k,n,j)-v(k,i,n)*
     1work2(k,n,j)*8.0d0
170    continue
       do 172 k=1,npar
172    zno(k,i,j)=work3(k,i,j)
171    continue
       re=ro+h
       call potval(re,h,v,erel,cent,unit,lam,coupindx,
     1             r,pot,ypot,alpha,z,z2,
     2             ea,delta,cinf,impot,nch,ndat,npar,npot,
     3             xi,r2av,jstart)
       do 180 i=1,nch
       do 180 j=1,nch
       do 179 k=1,npar
       work1(k,i,j)=unit(i,j)
       work2(k,i,j)=unit(i,j)
179    zne(k,i,j)=unit(i,j)+zno(k,i,j)
180    continue
       call gaussv(zne,work1,sum,nch,nol,n1,n4)
       do 191 i=1,nch
       do 191 j=1,nch
       do 188 k=1,npar
188    work2(k,i,j)=0.0d0
       do 190 n=1,nch
       do 189 k=1,npar
189    work2(k,i,j)=work2(k,i,j)+zno(k,i,n)*work1(k,n,j)
     1-4.0d0*v(k,i,n)*unit(n,j)
190    continue
       do 192 k=1,npar
192    zne(k,i,j)=work2(k,i,j)
191    continue
200    continue
c last integration point
        do 201 i=1,nch
        do 201 j=1,nch
        do 202 k=1,npar
        zne(k,i,j)=(zne(k,i,j)+2.0d0*v(k,i,j))/h
202     yn(k,i,j)=zne(k,i,j)
201     continue
c       write(21,*) yn(1,1,1),yn(1,2,2)
c       write(21,*) yn(1,1,2),yn(1,2,1)
       return
       end
c********************************************************************ZLB
c 11-5-94  modified for arbitrary number of channels, pcs
c---------------------------------------------------------
c nch    number of channels
c npar   number of partial waves per pass
c---------------------------------------------------------
      subroutine kmatrix(ylog,rjl,rjlp,rnl,rnlp,jv,rki,nprint,
     1                   unit,fkmat,smat,parsig,nch,ndat,npar,npot)
      implicit real*8 (a-h,o-z)
      integer nch,npar,nol,n4,nprint
czlb      parameter(nch=3,npar=64)
c shared arrays
       double precision ylog(npar,nch,nch),rjl(npar,nch,nch)
     1,rnl(npar,nch,nch),rjlp(npar,nch,nch),rnlp(npar,nch,nch)
     1,rki(nch)
      double precision jv(npar),unit,fkmat,smat,parsig
czlb
      dimension unit(nch,nch)
      dimension fkmat(npar,nch,nch),smat(npar,nch,nch),
     1          parsig(npar,nch,nch)
      dimension smr(npar,nch,nch),smi(npar,nch,nch)
c
c      common/une/unit(nch,nch)
c      common/fmtrix/ fkmat(npar,nch,nch),smat(npar,nch,nch)
c     . ,parsig(npar,nch,nch)
c       common/smatrix/ smr(npar,nch,nch),smi(npar,nch,nch)
czlb
c local arrays
       double precision fnum(npar,nch,nch),fden(npar,nch,nch),sum(npar)
     . ,yr(npar,nch,nch),yl(npar,nch,nch),smr,smi
     . ,u(npar,nch,nch)
     . ,term1ij,term2ij,x,sqs,pi
czlb       data nol/npar/,n4/npar/
       nol = npar
       n4  = npar
czlb
       pi=4.0d0*datan(1.0d0)
       n1=nch
       do 350 k=1,npar
       do 340 i=1,nch
       do 340 j=1,nch
       term1ij=0.0d0
       term2ij=0.0d0
       do 329 m=1,nch
      term1ij=term1ij+ylog(k,i,m)*rjl(k,m,j)-rjlp(k,i,m)*unit(m,j)
329   term2ij=term2ij+rnlp(k,i,m)*unit(j,m)-ylog(k,i,m)*rnl(k,m,j)
       fnum(k,i,j)=term1ij
340    fden(k,i,j)=term2ij
350    continue
       call gaussv(fden,fnum,sum,nch,nol,n1,n4)
c*****k-matrix=w.  now calculate the s-matrix.
       do 4 j=1,nch
       do 4 i=1,nch
       x=unit(i,j)
       do 6 k=1,npar
       fkmat(k,i,j)=fnum(k,i,j)
       yr(k,i,j)=x
6      sum(k)=x
       do 7 m=1,nch
       do 7 k=1,npar
7      sum(k)= sum(k) + fnum(k,i,m)*fnum(k,m,j)
       do 5 k=1,npar
5      u(k,i,j)=sum(k)
4      continue
       call  gaussv(u,yr,sum,nch,nol,n1,n4)

       do 12 j=1,nch
       do 12 i=1,nch
       do 13 k=1,npar
13     sum(k)=0.0d0
       do 14 m=1,nch
       do 14 k=1,npar
14     sum(k)= sum(k) + 2.0d0*yr(k,i,m)*fnum(k,m,j)
       do 15 k=1,npar
15     yl(k,i,j)=sum(k)
12     continue

       do 16 j=1,nch
       do 17 i=1,nch
       do 17 k=1,npar
17     yr(k,i,j) = -unit(i,j) +2.0d0*yr(k,i,j)
16     continue
c*****yr=re(s) and yl=im(s)

c...   calculate cross sections
       do 52 i=1,nch
       do 52 j=1,nch
       do 53 k=1,npar
       sqs=(unit(i,j)-yr(k,i,j))**2+yl(k,i,j)**2
       smat(k,i,j)=sqs
       smr(k,i,j)=yr(k,i,j)
       smi(k,i,j)=yl(k,i,j)
53     parsig(k,i,j)=(sqs)*pi*(2.0d0*jv(k)+1.0d0)/(rki(i)**2)
52     continue
c get output routine
       if(nprint.eq.2) call outdat(fkmat,smat,parsig,smr,smi,unit,
     1                             nch,ndat,npar,npot,jv)
       return
       end
c********************************************************************ZLB
c 11-5-94  modified for arbitrary number of channels, pcs
c      output kmatrix and S*S matrix
c      modify format statements 39 and 45 as required.
c---------------------------------------------------------
c nch    number of channels
c npar   number of partial waves per pass
c--------------------------------------------------------
       subroutine outdat(fkmat,smat,parsig,smr,smi,unit,
     1                   nch,ndat,npar,npot,jv)
       implicit real*8 (a-h,o-z)
       integer nch,npar
czlb       parameter(nch=3,npar=64)
      dimension unit(nch,nch)
      dimension fkmat(npar,nch,nch),smat(npar,nch,nch),
     1          parsig(npar,nch,nch)
      dimension smr(npar,nch,nch),smi(npar,nch,nch)
c
c      common/une/unit(nch,nch)
c      common/fmtrix/ fkmat(npar,nch,nch),smat(npar,nch,nch)
c     . ,parsig(npar,nch,nch)
c       common/smatrix/ smr(npar,nch,nch),smi(npar,nch,nch)
czlb
       double precision  jv(1),unit,fkmat,smat,parsig,smr,smi

       do 30 k=1,npar
       write(2,33)jv(k)
       write(2,36)
       do 29 i=1,nch
       write(2,34)(fkmat(k,i,j),j=1,nch)
29       continue
       write(2,35)
       do 28 i=1,nch
      write(2,34)(smat(k,i,j),j=1,nch)
28        continue
c
       write(2,44)
       do 18 i=1,nch
       write(2,34)(smr(k,i,j),j=1,nch)
18        continue
c
       write(2,46)
       do 19 i=1,nch
       write(2,34)(smi(k,i,j),j=1,nch)
19        continue
30     continue
33     format( /,' J= ',f10.1)
34     format(10g12.3)
35     format( ' - - - - - T*T MATRIX - - - - - -')
36      format( ' - - - - - K MATRIX - - - - - - -')
44     format( ' - - - - - Re S MATRIX - - - - - -')
46      format( ' - - - - - Im S MATRIX - - - - - - -')
c output charge transfer cross section
c sigma(3-1),sigma(3-2)
       write(4,39)
39     format( '                        ')

       do 40 k=1,npar

       write(4,33)jv(k)
c********writes may require modification*************
       write(4,41)(smat(k,nch,i),i=1,nch)
c     %,(parsig(k,nch,i),i=nch-1,1,-1)
c output graph data
c       write(7,42)jv(k),(parsig(k,nch,i),i=nch-1,1,-1)
40     continue
c *********formats may require modification **********
41     format( ' S*S ',10g12.4)
c42     format(6g15.5)

c45     format( ' TOTAL ',/, ' SIGMA(2-1)'
c     %,g15.4,
c     %' TOTAL ',g15.4)
       return
       end
c********************************************************************ZLB
c 12-16-03 generalized, pcs 
c  3-4-96 modified for all potentials to go to zero asymptotically
c  2-9-95 modified for SiH2+, two channels, pcs
c 11-5-94 modified by pcs
c         subroutine to return diabatic potential at
c         a given internuclear distance x.
c         Some editing of short range fits is required.
c         See below. If outgoing channels are not coulomb,
c         i.e., for a singly charged system, the long range
c         fits will also need modification.
c----------------------------------------------------------
c nch    number of channels
c ndat   number of potential data points per channel
c npar   number of partial waves per pass
c npot   size of potential array, includes diagonal
c        elements and upper triangle
c-----------------------------------------------------------
        subroutine dpot(x,v,r,pot,ypot,alpha,z,z2,
     1                  ea,delta,cinf,impot,nch,ndat,npar,npot,
     2                  xi,r2av)
      implicit real*8 (a-h,o-z)
        integer nch, ndat,npot,npar,i,j,k,l,impot,
     .        zc,zn,nprint
czlb        parameter(nch=3,ndat=31,npar=64)
czlb        parameter(npot=nch*(nch+1)/2)
czlb
      dimension pot(ndat,npot),ypot(ndat,npot)
      dimension r(ndat)
      dimension ea(nch),cinf(npot),impot(nch)
      dimension alpha(nch),z(nch),z2(nch)
      common/fit/rs,rl,nprint
czlb      06/27/04
      dimension xi(nch),r2av(nch)
czlb 
c
c
c      common/dat/pot(ndat,npot),ypot(ndat,npot)
c      common/grid/r(ndat)
c      common/fit/rs,rl,nprint
c      common/fit1/ea(nch),delta,cinf(npot),impot(nch)
c      common/fit2/alpha(nch),z(nch),z2(nch)
czlb
c
        double precision v(npar,nch,nch),u(npot)
        double precision pott(ndat),ypott(ndat)
     .   ,y,r,pot,ypot,a,b,dem,delta,rs,rl,alpha,
     .    z,x,ea,z2,cinf
c
c rename spline fit data
c
       do 11 j=1,npot
       do 10 i=1,ndat
        ypott(i)=ypot(i,j)
10      pott(i)=pot(i,j)
c
c if internuclear distance is outside of spline range
c skip spline call
c
        if((x.gt.rl).or.(j.gt.nch .and. x.gt.rl)) then
          y=0.0d0
          goto 12
        end if
c
c spline fit all couplings and potentials
c
        call splint(r,pott,ypott,ndat,x,y)
12      u(j)=y
c
c get short range or long range fits
c
c short range potentials of form a*exp(-b*r) + Vmin
c pott(ndat)=Vmin is last data point
c  or asymptotic energy for coulomb channels  and
c  potential at equilibrium distance for incoming
c  channel. Adjust as required.
c
       if(x.lt.rs) then 
c
c   potentials
c
         if(j.le.nch) then
         b=-(pott(2)-pott(1))/(r(2)-r(1))/(pott(1)-pott(impot(j)))
         a=(pott(1)-pott(impot(j)))*dexp(b*rs)
         u(j)=a*dexp(-b*x)+pott(impot(j))
c         write(20, *) 'coulomb',a,b,pott(ndat),u(j)
         end if
c
c
c short range form for diabatic off diagonal potential
c  which goes to zero at R=0. use form a*r^2+b*R.
c  coefficients determined by Cramers' rule. 
         if(j.gt.nch) then
          dem=r(1)**2*r(2)-r(2)**2*r(1)
          a=(r(2)*pott(1)-r(1)*pott(2))/dem
          b=(r(1)**2*pott(2)-r(2)**2*pott(1))/dem
          u(j)=a*x**2+b*x
czlb
c           r0=r(1)*0.95
c           p0=pott(1)*0.95d0
c
c           a=(pott(1)/sin(0.5*r(1))-p0/sin(0.5*r0))/(r(1)-r0)
c           b=p0/sin(0.5*r0)-a*r0
c           u(j)=(a*x+b)*sin(0.5*x)
czlb
               
c        write(20, *) 'off0',a,b,dem, u(j)
         end if
       end if
c        write(21,*) 'x,u',x,u(j)
11      continue
c
c add constants so that all channels go to zero
c at R=infinity. delta is total electronic energy
c of channel 1 at infinity. ea(i) are excitation
c energies of each channel above channel 1 at
c R=infinity. Must use ab initio data. 
c
c Also create v(k,i,j) potential matrix for indices i,j
c and partial wave k
c 
        do 30 k=1,npar
c diagonal diabats
           do 32 i=1,nch
              v(k,i,i)=u(i)+delta-ea(i)
32         continue
c
c    now make long-range fits for the diabatic potentials         
c 
           zc=0
           if(x.gt.rl) then
              do 31 i=1,nch
c    check if Coulomb channel
                 if (z2(i).ne.0.0d0) then
                    v(k,i,i)=(z(i)*z2(i))/x 
                    zc=zc+1
                    goto 31
                 else
                    zn=i
                 end if
c    for ion-neutral polarization potential
czlb              v(k,i,i)=-alpha(i)*z(i)**2/x**4/2.0d0
               v(k,i,i)=-alpha(i)*z(i)**2/x**4/2.0d0
     .                  +xi(i)*r2av(i)/x**3
31            continue
c    for ion-ion Coulomb channels 
              if (zc.gt.0) then
czlb                 print*, 'Coulomb channels,zn=', zn,zc
              end if
c    for couplings
              do 36 j=1+nch,npot
              b=-(pot(ndat,j)-pot(ndat-1,j))
     .          /(r(ndat)-r(ndat-1))/pot(ndat,j)
c               b=0.0
czlb             a=(pot(ndat,j)-cinf(j))*dexp(b*rl)
czlb             u(j)=a*dexp(-b*x)+cinf(j)
c
              a=(pot(ndat,j)-cinf(j))
              u(j)=a*dexp(-b*(x-rl))+cinf(j)
c              write(6, 40) j,a,b,pot(ndat,j),u(j)
36            continue             
           end if
40         format(i3,1x,4e15.5)
c off-diagonal diabatic couplings
           l=nch+1
           do 33 i=1,nch-1
              do 34 j=i+1,nch
                v(k,i,j)=u(l)
                v(k,j,i)=u(l)
              l=l+1
34            continue
33         continue
30      continue
        if(nprint.eq.3) then
           write(60,999) x,(v(1,i,i),i=1,nch)
             do 35 i=1,nch-1
                write(60+i,999) x,(v(1,i,j),j=i+1,nch)
35           continue
c        write(41,999) x,v(1,1,2)
c        write(23,999) x,v(1,3,4)
         do ii=1,nch
         do jj=1,nch
           if(abs(v(1,ii,jj)).lt.1.0d-99) v(1,ii,jj) =0.0d0
         enddo
         enddo
         write(400,999) x,v(1,1,1),v(1,2,2),v(1,3,3),v(1,4,4),v(1,5,5),
     4                                               v(1,6,6),v(1,7,7) 
         write(100,999) x,v(1,2,1),
     1                    v(1,3,1),v(1,3,2),
     2                    v(1,4,1),v(1,4,2),v(1,4,3)
         write(200,999) x,v(1,5,1),v(1,5,2),v(1,5,3),v(1,5,4),
     3                    v(1,6,1),v(1,6,2),v(1,6,3),v(1,6,4),v(1,6,5)
         write(300,999) x,v(1,7,1),v(1,7,2),v(1,7,3),v(1,7,4),v(1,7,5),
     4                                               v(1,7,6)   
        end if
999     format(11g14.6)
        return
        end
c********************************************************************ZLB
c 11-5-94 modified by pcs for n-channel problem
c potval returns the potentials u(npar,nch,nch)
c at radius r, for an nch channel problem. 
c h is the step size.
c all other parameters are in common block POT.
c------------------------------------------------
c nch    number of channels
c npar   number of partial waves per pass
c------------------------------------------------
      subroutine potval(r,h,u,erel,cent,unit,lam,coupindx,
     1                  ra,pot,ypot,alpha,z,z2,
     2                  ea,delta,cinf,impot,nch,ndat,npar,npot,
     3                  xi,r2av,jstart)
      implicit real*8 (a-h,o-z)
      integer nch, npar, coupindx
czlb      parameter(nch=3,npar=64,ndat=31)
czlb        parameter(npot=nch*(nch+1)/2)
      double precision v(npar,nch,nch)
c shared arrays
      double precision u(npar,nch,nch),
     . twomu,e,erel,cent,scale,unit,r,lam,lamm,jstart
czlb
      dimension erel(nch),cent(npar,nch)
      dimension unit(nch,nch)
      dimension lam(nch),coupindx(nch,nch)
      dimension pot(ndat,npot),ypot(ndat,npot)
      dimension ra(ndat)
      dimension alpha(nch),z(nch),z2(nch)
      dimension ea(nch),cinf(npot),impot(nch)
c
      common/pot/twomu,e,scale
      common/fit/rs,rl,nprint
czlb        06/27/04
      dimension xi(nch),r2av(nch)
czlb 
c
c      common/pot/twomu,e,erel(nch),cent(npar,nch),scale
c      common/une/unit(nch,nch)
c      common/lamm/lam(nch),coupindx(nch,nch)
czlb
c potential is rescaled by alpha parameter.
       call dpot(r,v,ra,pot,ypot,alpha,z,z2,
     1           ea,delta,cinf,impot,nch,ndat,npar,npot,
     2           xi,r2av)
C  CALCULATE POTENTIAL AT R = T AND SCALE BY -H**2/6
c note: v already contains value of threshold energies
      do 10 k=1,npar
c rotational couplings
c       v(k,1,2)=(2.0d0/twomu/r**2*(cent(k,1))**0.5)*v(k,1,2)
c       v(k,2,1)=v(k,1,2)
c       v(k,2,3)=(2.0d0/twomu/r**2*(cent(k,1))**0.5)*v(k,2,3)
c       v(k,3,2)=v(k,2,3)
       do 10 i=1,nch
       do 10 j=1,nch
          if(coupindx(i,j).eq.2) then
             lamm=min(lam(i),lam(j))
             if((dfloat(k-1)+jstart).ge.lamm) then
                v(k,i,j)=2.0d0/twomu/r**2*v(k,i,j)
     .           *((dfloat(k-1)+jstart-lamm)
     .           *(dfloat(k-1)+jstart+lamm+1.0d0))**0.5       
c           if(r.lt.0.055) print*,dfloat(k-1)+jstart,i,j,lamm,
c     .            ((dfloat(k-1)+jstart-lamm)
c     .          *(dfloat(k-1)+jstart+lamm+1.0d0))**0.5,
c     .           v(k,i,j)
             else
               v(k,i,j)=0.0d0
             end if        
          end if
          if(coupindx(i,j).gt.2)  then
             v(k,i,j)=0.0d0
c          print*,k,i,j
          end if
cnew 12/22/2004
c pcs, to ensure J(J+1)-lambda^2 >0, i.e. no contributions from these
c      partial waves due to radial coupling
          if((coupindx(i,j).eq.1).and.(cent(k,i).lt.lam(i)**2)) 
     .         v(k,i,j)=0.0d0
c
c all potentials go to zero NOT THRESHOLD Energy. Put that into erel(i)
      u(k,i,j)=-scale*(twomu*(v(k,i,j)-erel(i)*unit(i,j))+
     1unit(i,j)*(cent(k,i)-lam(i)**2)/(r*r))
      u(k,j,i)=u(k,i,j)
cnew 12/22/2004
10    continue
c      write(23,*)r,v(1,5,5)+420.0d0/(r*r)/twomu,
c     .   v(1,5,5)+650.0d0/(r*r)/twomu
c      write(20,*)r,v(1,1,1)-erel(1)+e,v(1,2,2)-erel(2)+e
c      write(21,*)r,v(1,3,3)-erel(3)+e,v(1,4,4)-erel(4)+e
c      write(22,*)r,v(1,3,1),v(1,3,2)
c      write(26,*)r,v(1,3,4),v(1,3,5)
c      write(27,*)r,v(1,4,5)
      return
      end
c********************************************************************ZLB
      subroutine gaussv (a,b,s, n,nop, n1,n4)
      implicit real*8 (a-h,o-z)
c*****solves n by n matrix equation  by gaussian elimination:
c...      a(p,i,j)x(p,j,k)=b(p,i,k)
c...   with vectorisation on passive index p<=nol
c...      s(nop) is scratch space
c*****answer is returned in b, a is destroyed.
       integer  p,n,n1,n4,nop
       double precision a(n4,n1,n1),b(n4,n1,n1), s(1)
      if(n.eq.1) go to 280
 2030  format("  gauss n nop n1 n4 ",5i5)
      ns = n - 1
      do 200 i=1,ns
      ig = i + 1
       do 1 p=1,nop
    1  s(p)=1.0d0/a(p,i,i)
      do 120 j=ig,n
       do 120 p=1,nop
  120  a(p,i,j)=a(p,i,j)*s(p)
      do 130 j=1,n
       do 130 p=1,nop
  130  b(p,i,j)=b(p,i,j)*s(p)
      do 150 k=1,n
      if(k-i) 151,150,151
  151  do 2 p=1,nop
    2  s(p)=a(p,k,i)
      do 140 l=ig,n
       do 140 p=1,nop
  140  a(p,k,l)=a(p,k,l) - a(p,i,l)*s(p)
      do 145 l=1,n
       do 145 p=1,nop
  145  b(p,k,l)=b(p,k,l) - b(p,i,l)*s(p)
150   continue
200   continue
       do 3 p=1,nop
    3  s(p)=1.0d0/a(p,n,n)
      do 230 j=1,n
       do 230 p=1,nop
  230  b(p,n,j)=b(p,n,j)*s(p)
      do 250 k=1,ns
       do 4 p=1,nop
    4  s(p)=a(p,k,n)
      do 250 l=1,n
       do 250 p=1,nop
  250  b(p,k,l)=b(p,k,l) - b(p,n,l)*s(p)
      return
  280  do 5 p=1,nop
    5  b(p,1,1)=b(p,1,1)/a(p,1,1)
      return
      end
c********************************************************************ZLB
      SUBROUTINE COULFG(XX,ETA1,XLMIN,XLMAX, FC,GC,FCP,GCP,
     X                  MODE1,KFN,IFAIL)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C                                                                      C
C  Revised Coulomb wavefunction program using Steed's method           C
C                                                                      C
C  A. R. Barnett           Manchester  March   1981                    C
C                                                                      C
C  original program 'RCWFN'      in    CPC  8 (1974) 377-395           C
C                 + 'RCWFF'      in    CPC 11 (1976) 141-142           C
C  full description of algorithm in    CPC 21 (1981) 297-314           C
C  this version written up       in    CPC 27 (1982) 147-zzz           C
C                                                                      C
C  COULFG returns F,G,F',G', for real xx.GT.0,real eta1 (including 0), C
C   and real lambda(xlmin) .GT. -1 for integer-spaced lambda values    C
C   thus giving positive-energy solutions to the Coulomb Schrodinger   C
C   equation,to the Klein-Gordon equation and to suitable forms of     C
C   the Dirac equation ,also spherical & cylindrical Bessel equations  C
C                                                                      C
C  for a range of lambda values (xlmax - xlmin) must be an integer,    C
C  starting array element is m1 = max0(idint(xlmin+accur),0) + 1       C
C      see text for modifications for integer l-values                 C
C                                                                      C
C  if 'MODE' = 1  get F,G,F',G'   for integer-spaced lambda values     C
C            = 2      F,G      unused arrays must be dimensioned in    C
C            = 3      F               call to at least length (1)      C
C  if 'KFN'  = 0 real        Coulomb functions are returned            C
C            = 1 spherical   Bessel      "      "     "                C
C            = 2 cylindrical Bessel      "      "     "                C
C  the use of 'MODE' and 'KFN' is independent                          C
C                                                                      C
C  Precision:  results to within 2-3 decimals of 'machine accuracy'    C
C   in oscillating region x .GE. eta1 + sqrt(eta1**2 + xlm(xlm+1))     C
C   COULFG is coded for REAL*8 on IBM or equivalent  accur = 10**-16   C
C   use AUTODBL + extended precision on HX compiler  accur = 10**-33   C
C   for mantissas of 56 & 112 bits. For single precision CDC (48 bits) C
C   reassign SQRT=SQRT etc.  See text for complex arithmetic version  C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION    FC(1),GC(1),FCP(1),GCP(1)
      LOGICAL      ETANE0,XLTURN
C***  COULFG has calls to: SQRT,ABS,MOD,INT,SIGN,DFLOAT,MIN1
      DATA ZERO,ONE,TWO,TEN2,ABORT /0.0D0, 1.0D0, 2.0D0, 1.0D2, 2.0D4/
      DATA HALF,TM30 / 0.5D0, 1.0D-30 /
      DATA RT2DPI /0.79788 45608 02865 35587 98921 19868 76373 D0/
C *** This constant is  SQRT(TWO/PI):  use Q0 for IBM REAL*16: D0 for
C ***  REAL*8 & CDC double p:  E0 for CDC single p; and truncate value.
C
                        ACCUR = 1.0D-16
C ***            change accur to suit machine and precision required
      MODE  = 1
      IF(MODE1 .EQ. 2 .OR. MODE1 .EQ. 3 ) MODE = MODE1
      IFAIL = 0
      IEXP  = 1
      NPQ   = 0
      ETA   = ETA1
      GJWKB = ZERO
      PACCQ = ONE
      IF(KFN .NE. 0) ETA = ZERO
                 ETANE0  = ETA .NE. ZERO
      ACC   = ACCUR
      ACC4  = ACC*TEN2*TEN2
      ACCH  = dSQRT(ACC)
C ***    test range of XX, exit if.LE.SQRT(ACCUR) or if negative
C
      IF(XX .LE. ACCH)                          GO TO 120
      X     = XX
      XLM   = XLMIN
      IF(KFN .EQ. 2)  XLM = XLM - HALF
      IF(XLM .LE. -ONE .OR. XLMAX .LT. XLMIN)   GO TO 130
      E2MM1 = ETA*ETA + XLM*XLM + XLM
      XLTURN= X*(X - TWO*ETA) .LT. XLM*XLM + XLM
      DELL  = XLMAX - XLMIN + ACC
      IF(dABS(dMOD(DELL,ONE)) .GT. ACC) WRITE(6,1060)XLMAX,XLMIN,DELL
      LXTRA = idINT(DELL)
      XLL   = XLM + DFLOAT(LXTRA)
C ***       LXTRA is number of additional lambda values to be computed
C ***       XLL  is max lambda value, or 0.5 smaller for J,Y Bessels
C ***         determine starting array element (M1) from XLMIN
      M1  = MAX0(idINT(XLMIN + ACC),0) + 1
      L1  = M1 + LXTRA
C
C ***    evaluate cf1  =  f   =  fprime(xl,eta,x)/f(xl,eta,x)
C
      XI  = ONE/X
      FCL = ONE
      PK  = XLL + ONE
      PX  = PK  + ABORT
   10 EK  = ETA / PK
      F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
      PK1 =  PK + ONE
C ***   test ensures b1 .ne. zero for negative eta; fixup is exact.
             IF(dABS(ETA*X + PK*PK1) .GT. ACCH)  GO TO 20
             FCL  = (ONE + EK*EK)/(ONE + (ETA/PK1)**2)
             PK   =  TWO + PK
      GO TO 10
   20 D   =  ONE/((PK + PK1)*(XI + EK/PK1))
      DF  = -FCL*(ONE + EK*EK)*D
            IF(FCL .NE. ONE )  FCL = -ONE
            IF(D   .LT. ZERO)  FCL = -FCL
      F   =  F  + DF
C
C ***   begin cf1 loop on pk = k = lambda + 1
C
      P     = ONE
   30 PK    = PK1
        PK1 = PK1 + ONE
        EK  = ETA / PK
        TK  = (PK + PK1)*(XI + EK/PK1)
        D   =  TK - D*(ONE + EK*EK)
              IF(dABS(D) .GT. ACCH)             GO TO 40
              WRITE (6,1000) D,DF,ACCH,PK,EK,ETA,X
              P = P  +   ONE
              IF( P .GT. TWO )                  GO TO 140
   40 D     = ONE/D
              IF (D .LT. ZERO) FCL = -FCL
        DF  = DF*(D*TK - ONE)
        F   = F  + DF
              IF(PK .GT. PX)                    GO TO 140
      IF(dABS(DF) .GE. dABS(F)*ACC)             GO TO 30
                  NFP = PK - XLL - 1
      IF(LXTRA .EQ. 0)                          GO TO 60
C
C *** downward recurrence to lambda = XLM. Array GC,if present,stores RL
C
      FCL = FCL*TM30
      FPL = FCL*F
      IF(MODE .EQ. 1) FCP(L1) = FPL
                      FC (L1) = FCL
      XL  = XLL
      RL  = ONE
      EL  = ZERO
      DO 50  LP = 1,LXTRA
         IF(ETANE0) EL = ETA/XL
         IF(ETANE0) RL = dSQRT(ONE + EL*EL)
         SL    =  EL  + XL*XI
         L     =  L1  - LP
         FCL1  = (FCL *SL + FPL)/RL
         FPL   =  FCL1*SL - FCL *RL
         FCL   =  FCL1
         FC(L) =  FCL
         IF(MODE .EQ. 1) FCP(L)  = FPL
         IF(MODE .NE. 3 .AND. ETANE0) GC(L+1) = RL
   50 XL = XL - ONE
      IF(FCL .EQ. ZERO) FCL = ACC
      F  = FPL/FCL
C ***    now we have reached lambda = XLMIN = XLM
C ***    evaluate CF2 = P + I.Q  again using Steed's algorithm
C ***    see text for compact complex code for sp CDC or non-ANSI IBM
C
   60 IF( XLTURN ) CALL JWKB(X,ETA,dMAX1(XLM,ZERO),FJWKB,GJWKB,IEXP)
      IF( IEXP .GT. 1 .OR. GJWKB .GT. ONE/(ACCH*TEN2))  GO TO 80
          XLTURN = .FALSE.
      TA =  TWO*ABORT
      PK =  ZERO
      WI =  ETA + ETA
      P  =  ZERO
      Q  =  ONE - ETA*XI
      AR = -E2MM1
      AI =  ETA
      BR =  TWO*(X - ETA)
      BI =  TWO
      DR =  BR/(BR*BR + BI*BI)
      DI = -BI/(BR*BR + BI*BI)
      DP = -XI*(AR*DI + AI*DR)
      DQ =  XI*(AR*DR - AI*DI)
   70 P     = P  + DP
         Q  = Q  + DQ
         PK = PK + TWO
         AR = AR + PK
         AI = AI + WI
         BI = BI + TWO
         D  = AR*DR - AI*DI + BR
         DI = AI*DR + AR*DI + BI
         C  = ONE/(D*D + DI*DI)
         DR =  C*D
         DI = -C*DI
         A  = BR*DR - BI*DI - ONE
         B  = BI*DR + BR*DI
         C  = DP*A  - DQ*B
         DQ = DP*B  + DQ*A
         DP = C
         IF(PK .GT. TA)                         GO TO 150
      IF(dABS(DP)+dABS(DQ).GE.(dABS(P)+dABS(Q))*ACC)   GO TO 70
                      NPQ   = PK/TWO
                      PACCQ = HALF*ACC/dMIN1(dABS(Q),ONE)
                      IF(dABS(P) .GT. dABS(Q)) PACCQ = PACCQ*dABS(P)
C
C *** solve for FCM = F at lambda = XLM,then find norm factor W=W/FCM
C
      GAM = (F - P)/Q
            IF(Q .LE. ACC4*dABS(P))             GO TO 160
      W   = ONE/dSQRT((F - P)*GAM + Q)
            GO TO 90
C *** arrive here if G(XLM) .GT. 10**6 or IEXP .GT. 70 & XLTURN = .TRUE.
   80 W   = FJWKB
      GAM = GJWKB*W
      P   = F
      Q   = ONE
C
C *** normalise for spherical or cylindrical Bessel functions
C
   90                     ALPHA = ZERO
          IF(KFN  .EQ. 1) ALPHA = XI
          IF(KFN  .EQ. 2) ALPHA = XI*HALF
                          BETA  = ONE
          IF(KFN  .EQ. 1) BETA  = XI
          IF(KFN  .EQ. 2) BETA  = dSQRT(XI)*RT2DPI
      FCM  = dSIGN(W,FCL)*BETA
           FC(M1)  = FCM
                      IF(MODE .EQ. 3)           GO TO 100
           IF(.NOT. XLTURN)   GCL =  FCM*GAM
           IF(      XLTURN)   GCL =  GJWKB*BETA
           IF( KFN .NE. 0 )   GCL = -GCL
           GC(M1)  = GCL
           GPL =  GCL*(P - Q/GAM) - ALPHA*GCL
                      IF(MODE .EQ. 2)           GO TO 100
           GCP(M1) = GPL
           FCP(M1) = FCM*(F - ALPHA)
  100 IF(LXTRA .EQ. 0 ) RETURN
C *** upward recurrence from GC(M1),GCP(M1)  stored value is RL
C *** renormalise FC,FCP at each lambda and correct regular derivative
C ***    XL   = XLM here  and RL = ONE , EL = ZERO for Bessels
         W    = BETA*W/dABS(FCL)
         MAXL = L1 - 1
      DO 110 L = M1,MAXL
                      IF(MODE .EQ. 3)           GO TO 110
                      XL = XL + ONE
         IF(ETANE0)   EL = ETA/XL
         IF(ETANE0)   RL = GC(L+1)
                      SL = EL + XL*XI
         GCL1     = ((SL - ALPHA)*GCL - GPL)/RL
         GPL      =   RL*GCL -  (SL + ALPHA)*GCL1
         GCL      = GCL1
         GC(L+1)  = GCL1
                      IF(MODE .EQ. 2)           GO TO 110
         GCP(L+1) = GPL
         FCP(L+1) = W*(FCP(L+1) - ALPHA*FC(L+1))
  110 FC(L+1)     = W* FC(L+1)
      RETURN
 1000 FORMAT(/,' CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ',1P,
     X       7D9.2,/)
C
C ***    error messages
C
  120 IFAIL = -1
      WRITE(6,1010) XX,ACCH
 1010 FORMAT(' FOR XX = ',1P,D12.3,' TRY SMALL-X  SOLUTIONS',
     X' OR X NEGATIVE',/ ,' SQUARE ROOT ACCURACY PARAMETER =  ',D12.3,/)
      RETURN
  130 IFAIL = -2
      WRITE (6,1020) XLMAX,XLMIN,XLM
 1020 FORMAT(/,' PROBLEM WITH INPUT ORDER VALUES:XLMAX,XLMIN,XLM = ',
     X1P,3D15.6,/)
      RETURN
  140 IFAIL =  1
      WRITE (6,1030) ABORT,F ,DF,PK,PX,ACC
 1030 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',F10.0,' ITERATIONS',/
     X,' F,DF,PK,PX,ACCUR =  ',1P,5D12.3,//)
      RETURN
  150 IFAIL =  2
      WRITE (6,1040) ABORT,P,Q,DP,DQ,ACC
 1040 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',F7.0,' ITERATIONS',/
     X,' P,Q,DP,DQ,ACCUR =  ',1P,4D17.7,D12.3,//)
      RETURN
  160 IFAIL =  3
      WRITE (6,1050) P,Q,ACC,DELL,LXTRA,M1
 1050 FORMAT(' FINAL Q.LE.ABS(P)*ACC*10**4 , P,Q,ACC = ',1P,3D12.3,4X,
     X' DELL,LXTRA,M1 = ',D12.3,2I5,/)
      RETURN
 1060 FORMAT(' XLMAX - XLMIN = DELL NOT AN INTEGER ',1P,3D20.10,/)
      END
c********************************************************************ZLB
      SUBROUTINE JWKB(XX,ETA1,XL,FJWKB,GJWKB,IEXP)
      implicit real*8 (a-h,o-z)
      double precision    XX,ETA1,XL,FJWKB,GJWKB,DZERO,zero
     . ,half,one,six,ten,rl35,aloge,x,eta,gh2,xll1,hll,hl,sl,
czlb     .  rl2,gh,phi,phi10,iexp
     .  rl2,gh,phi,phi10
C *** Computes JWKB approximations to Coulomb functions    for DL.GE. 0
C *** as modified by Biedenharn et al. Phys Rev 97 (1955) 542-554
C *** calls MAX1,SQRT,ALOG,EXO,ATAN2,FLOAT,INT        Barnett Feb 1981
      DATA   ZERO,HALF,ONE,SIX,TEN/ 0.0d0, 0.5d0, 1.0d0, 6.0d0, 10.0d0 /
      DATA  DZERO, RL35, ALOGE  /0.0D0, 35.0d0, 0.43429 45 d0 /
      X     = XX
      ETA   = ETA1
      GH2   = X*(ETA + ETA - X)
      XLL1  = dMAX1(XL*XL + XL,DZERO)
      IF(GH2 + XLL1 .LE. ZERO) RETURN
       HLL  = XLL1 + SIX/RL35
       HL   = dSQRT(HLL)
       SL   = ETA/HL + HL/X
       RL2  = ONE + ETA*ETA/HLL
       GH   = dSQRT(GH2 + HLL)/X
       PHI  = X*GH - HALF*( HL*dLOG((GH + SL)**2/RL2) - dLOG(GH) )
          IF(ETA .NE. ZERO) PHI = PHI - ETA*dATAN2(X*GH,X - ETA)
      PHI10 = -PHI*ALOGE
      IEXP  =  idINT(PHI10)
      IF(IEXP .GT. 70) GJWKB = TEN**(PHI10 - dble(IEXP))
      IF(IEXP .LE. 70) GJWKB = dEXP(-PHI)
      IF(IEXP .LE. 70) IEXP  = 0
      FJWKB = HALF/(GH*GJWKB)
      RETURN
      END
c********************************************************************ZLB
c P.C. Stancil 10-8-92
c Determines parameter array for spline
c From Numerical Recipes p.88 (1986)
c
      subroutine spline(x,y,n,yp1,ypn,y2)
      implicit real*8 (a-h,o-z)
c
      integer i,n,k,nmax
czlb      parameter (nmax=100)
      parameter (nmax=5000)
      double precision yp1,ypn,sig,p,qn,un
      double precision x(n),y(n),y2(n),u(nmax)
c
      if (yp1 .gt. 0.99d30) then
c      if (yp1 .gt. 0.99d60) then
         y2(1)=0.0d0
         u(1)=0.0d0
      else
         y2(1)=-0.5d0
         u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
c
      do 10 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.d0
         y2(i)=(sig-1.0d0)/p
         u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     .         /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
 10   continue
c
c      if (ypn .gt. 0.99d30) then
      if (ypn .gt. 0.99d60) then
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
c      do i=1,n
c      write(11,999) x(i),y(i),y2(i)
c      end do
999   format(3g12.4)
      return
      end      
c********************************************************************ZLB
c P.C. Stancil 10-8-92
c Evaluates a function using spline array from spline.f
c From Numerical Recipes p.89 (1986)
c
      subroutine splint(xa,ya,y2a,n,x,y)
      implicit real*8 (a-h,o-z)
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
c********************************************************************ZLB
