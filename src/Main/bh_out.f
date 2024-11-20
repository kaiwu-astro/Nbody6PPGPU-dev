       subroutine bh_out(xs,vs,bodys)
*
*      Additional output of bh neighbour particles 
*
       include "common6.h"
       integer nb, cnt 
     &, ind_sort(nmax),var_sort(nmax)
*
*    use single precision
      real*4 tmp_r(nmax),tmp_v,tmp_a,tmp_adot,fs(3,nmax),
     &       fdots(3,nmax)
      REAL*4  XS(3,NMAX),VS(3,NMAX),BODYS(NMAX),RHOS(NMAX),AS(20)
*
      do i=ifirst,ntot 
*    create array of distances:
           tmp_r(i-ifirst+1) = sqrt(xs(1,i)**2+xs(2,i)**2+xs(3,i)**2)
*    Convert force & derivative to single precision:
        do j=1,3
           fs(j,i) = real( f(j,i) )
           fdots(j,i) = real( fdot(j,i) )        
        enddo
*
      enddo
*
*    Make sure time is updated
      ttot = time+toff
*    Number of BH neighbour particles
      nb = 100
*
*    Sort distances and create array of indexes
*    using indexx routine from num receipts library
      call indexx(ntot-ifirst+1, tmp_r, ind_sort)
*
       OPEN (UNIT=99,STATUS='UNKNOWN',FORM='FORMATTED',
     &   FILE='bh_nb.dat', ACCESS='APPEND')
*
*    put option advance=no to avoid newline character for each particle
          write (99,97,advance="no") ttot
   97     format (e9.3)     
*
      do i=1,nb
         ii = ind_sort(i)
         tmp_v = sqrt(vs(1,ii)**2+vs(2,ii)**2+vs(3,ii)**2)
         tmp_a = sqrt(fs(1,ii)**2+fs(2,ii)**2+fs(3,ii)**2)
         tmp_adot = sqrt(fdots(1,ii)**2+fdots(2,ii)**2+fdots(3,ii)**2)  
*          
       write (99, 98, advance="no") 
     &          i,name(ii), cmbh,bodys(ii), (xs(j,ii), j=1,3),tmp_r(ii),
     &          (vs(j,ii), j=1,3), tmp_v, (fs(j,ii), j=1,3), tmp_a,
     &          (fdots(j,ii), j=1,3), tmp_adot, step(ii), stepr(ii)
   98  format (i7.4, i7.5, 1x, 20(e15.6), 1x)   
      enddo 
*    Put newline character after time interval
          write (99,99) 
   99     format ()     
*
      call flush(99)
*
      return
*
      end
