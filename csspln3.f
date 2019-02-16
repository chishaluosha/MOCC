c P.C. Stancil 10-11-92
c creates spline for adiabatic d/dr matrix element 
c
c________________________________________
c
c
c $Id: csspln3.f,v 1.1 2014/02/15 18:57:13 lyonsdav Exp $
c
c SUMMARY
c
c REVISION HISTORY
c $Log: csspln3.f,v $
c Revision 1.1  2014/02/15 18:57:13  lyonsdav
c Initial revision
c
c

      subroutine csspln (n,e,cs,cs2,l12,l122,l23,l232)
c
      integer n,ie
      double precision csp1,cspn,dumb
      double precision e(n),cs(n),cs2(n),l23(n),l232(n)
     .    ,l12(n),l122(n)
      external spline
c
      open(unit=18, file ='adiab.new')
      open(unit=19, file ='ddr1.new')
c
c reads in data
      do 10 ie=1,n
         read (18,*) e(ie), cs(ie),l12(ie),dumb,dumb,dumb
         read (19,*) e(ie), l23(ie),dumb,dumb,dumb,dumb,dumb
 10   continue
c
      csp1=0.0d0
      cspn=0.0d0
c
      call spline(e,cs,n,csp1,cspn,cs2)
      call spline(e,l23,n,csp1,cspn,l232)
      call spline(e,l12,n,csp1,cspn,l122)
c
      return
      end

