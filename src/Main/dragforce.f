       SUBROUTINE DRAGFORCE(I,FIRR,FD,XI,VI)
*
*
*       Drag force & first derivative.
*       -------------------------
*
      INCLUDE 'common6.h'
      REAL*8  FIRR(3),FD(3),FDRAG(3),FDDRAG(3),XI(3),VI(3),
     &        R2DISC,R2DENS,RDISC,DDENS,RHODOT,VXRE,VYRE,VZRE,VRE,
     &        CKFDRAG,DVXRE,DVYRE,DVZRE,DVRE,KFDOT,EXCENTR1(NMAX),
     &        AORBIT(NMAX),BORBIT(NMAX),CORBIT(NMAX),INCLIN1(NMAX),
     &        DISCMOM,STARMOM,ANGLE1,INCLIN2(NMAX),RR0_S,
     &        RDOT,Rzero,sigma,R2_DELTA,R2zero,H2Rzero,S,BETA_S,
     &        RRzero, RR2zero, RDISC_3over,RDISC_5over, SQRT_MBH,
     &        alpha, ZETA, RRCRIT_ZETA,R_ZETA,R0_ZETA,
     &        RCRIT_ZETA,R_,Z2, R_ZETA_R, RRCRIT, VXDISC,VYDISC,
     &        r_lim2, h_lim2
*     &        ,xold(3), vold(3), vv_new, vv_old
*         xold and vold are used to keep track of x  and xdot
*
      INTEGER I, K, N0
*
*      LOGICAL LPR
*
*     Limiting factors for the accretion disk
      r_lim2 = 0.5*0.5
      h_lim2 = 5.0d-03**2
*
      ttot = time + toff
*
      t_diss_on = 0.125
      if (ttot.lt.t_diss_on) goto 10

*      if (ttot.eq.0.0) then
*      LPR=NAME(I).EQ.250
*
      Rzero=0.22
      S=4
      BETA_S=0.7
      alpha = 0.75
*
*      R_CRIT = 0.0257314
      R2zero=Rzero*Rzero
      H2Rzero=HZ*HZ*R2zero
*
*     Qzero = 0.01
      N0 = 4000
*
*      end if
*     
*      if (kz(19).gt.0) then
*     in case of stellar evolution:
*          st_cross_sec = 0.0
*          do k = ifirst, ntot
*             st_cross_sec = st_cross_sec + radius(k)**2
*          enddo
*          Q_DRAG = Qzero*st_cross_sec/R2zero
*      else 
*     no stellar evolution:    
*     Qdrag is calculated automatically:
           Q_DRAG = Qzero*(float(N0)/float(N))*LOG(0.4*N)/LOG(0.4*N0) 
*      endif
*
*
*              Auxiliary Calculations
      R2DISC=XI(1)**2+XI(2)**2
      Z2=XI(3)**2
*
      IF (r2disc.lt.r_lim2.and.z2.lt.h_lim2) then
*
      RDISC=DSQRT(R2DISC)
      R_=DSQRT(R2DISC+Z2)

      sigma=(2.0-ALPHA)*Mdisc/(TWOPI*dsqrt(TWOPI)*HZ*Rzero)
      RDOT=(XI(1)*VI(1)+XI(2)*VI(2))/RDISC
*	
*      New density profile 
*
      IF (RDISC.LT.R_CRIT) THEN 
*
      RRCRIT=RDISC/R_CRIT
      ZETA=1.0
      RRCRIT_ZETA=RRCRIT
      R_ZETA=1.0/R2DISC
      R0_ZETA=1.0/R2zero
      RCRIT_ZETA=1.0/(R_CRIT*R_CRIT)
*
      ELSE
*
      ZETA=0.0
      RRCRIT_ZETA=1.0
      R_ZETA=1.0
      R0_ZETA=1.0
      RCRIT_ZETA=1.0
*      
      END IF
*
      R_ZETA_R=R_ZETA/R_
*
      RRzero=RDISC/Rzero
      RR2zero=R2DISC/R2zero
      RR0_S=RRzero**S
      R2_DELTA=RRzero**(-alpha)*EXP(-BETA_S*RR0_S)/RRCRIT_ZETA
*
*   Disk density
      DDENS=sigma*R2_DELTA*DEXP(-Z2/(2.0*H2Rzero*RRCRIT_ZETA**2))
*
*   The density derivative      
        RHODOT=(-RDOT)*(ZETA+alpha+BETA_S*S*RR0_S)/RDISC
     &       -(VI(3)*XI(3)*R_ZETA-ZETA*(R_ZETA_R)*Z2*RDOT)
     &        /(H2Rzero*RCRIT_ZETA)          
*
      RDISC_3over=RDISC**(1.5)
      RDISC_5over=RDISC**(2.5)
      SQRT_MBH=SQRT(CMBH)
*
*   Velocity of the disc
      VXDISC=-SQRT_MBH*XI(2)/RDISC_3over
      VYDISC=SQRT_MBH*XI(1)/RDISC_3over
*
*   Relative velocity
      VXRE=VXDISC-VI(1)
      VYRE=VYDISC-VI(2)
      VZRE=-VI(3)
      VRE=DSQRT(VXRE*VXRE + VYRE*VYRE + VZRE*VZRE)
      CKFDRAG=Q_DRAG*DDENS*VRE
*
*              Mean of Drag Force
      FDRAG(1)=CKFDRAG*VXRE
      FDRAG(2)=CKFDRAG*VYRE
      FDRAG(3)=CKFDRAG*VZRE
*      	
*     The dissipation energy:
*      EDISS1=EDISS1+BODY(I)*(FDRAG(1)*VI(1)+FDRAG(2)*VI(2)+
*     &       FDRAG(3)*VI(3))*STEP(I) 
*     
      ediss1=ediss1+0.5*body(i)*((a_drag(1,i)+fdrag(1))*(x(1,i)-x0(1,i))
     &                        +(a_drag(2,i)+fdrag(2))*(x(2,i)-x0(2,i)) 
     &                        +(a_drag(3,i)+fdrag(3))*(x(3,i)-x0(3,i)) )
* 
*
*      do ii=1,3
*         adrag(ii) = fdrag(ii)
*      end do
*      
        DTR = TIME - T0R(I)
*
      DVXRE=-SQRT_MBH*VI(2)/RDISC_3over
     &      +SQRT_MBH*1.5*XI(2)*RDOT/RDISC_5over
     &      -F(1,I)
      DVYRE=-SQRT_MBH*VI(1)/RDISC_3over
     &      +SQRT_MBH*1.5*XI(1)*RDOT/RDISC_5over
     &      -F(2,I)
      DVZRE=-F(3,I)
*     
      DVRE=(VXRE*DVXRE+VYRE*DVYRE+VZRE*DVZRE)/VRE
*
      KFDOT=Q_DRAG*DDENS
*   
*              Mean of Drag Force Derivative
      FDDRAG(1)=KFDOT*(RHODOT*VXRE*VRE+DVRE*VXRE+DVXRE*VRE)
      FDDRAG(2)=KFDOT*(RHODOT*VYRE*VRE+DVRE*VYRE+DVYRE*VRE)
      FDDRAG(3)=KFDOT*(RHODOT*VZRE*VRE+DVRE*VZRE+DVZRE*VRE)
*
*      PRINT *, "FDDRAG 1, 2,3 =", FDDRAG(1), FDDRAG(2), FDDRAG(3)
*      FDDRAG(1)=KFDOT*(DDOTDENS*VXRE+DVRE*VXRE+DVXRE)
*      FDDRAG(2)=KFDOT*(DDOTDENS*VYRE+DVRE*VYRE+DVYRE)
*      FDDRAG(3)=KFDOT*(DDOTDENS*VZRE+DVRE*VZRE+DVZRE)
*     PRINT*,"fddragx, fddragy, fddragz", FDDRAG(1),FDDRAG(2),FDDRAG(3)
*   
*     END IF
*    
*      fdx = dsqrt(fdrag(1)**2+fdrag(2)**2+fdrag(3)**2)
*      fddx = dsqrt(fddrag(1)**2+fddrag(2)**2+fddrag(3)**2)
*      stepdd = fdx/fddx
*      if(stepdd.lt.step(i))print*,' i,dt,dtd=',i,step(i),stepdd
*
        do k=1,3
           a_drag(k,i) = fdrag(k)
        end do
*
      ELSE 
        do k=1,3
           a_drag(k,i) = 0.0
           fddrag(k) = 0.0
        end do
* end if (r_lim2 & h_lim2):
      END IF
*
*              Total Force Acting on a Star
      DO K=1,3
         FIRR(K)=FIRR(K)+A_DRAG(K,I)
         FD(K)=FD(K)+FDDRAG(K)
*
*   !!  Try experimental ediss calculation:
*     We need to add drag component to the total force and it's
*     derivative in order to use XVPRED
*         F(K,I) = F(K,I) + FDRAG(K)
*         FDOT(K,I) = FDOT(K,I) + FDDRAG(K)  
*         xold(k) = x(k,i)
*         vold(k) = xdot(k,i)
      END DO 
*
*   Predict coordinates and velocities
*   experimental feature
*      if (kdr.eq.0) then
*      call xvpred(i,nnb)
*      vv_old2 = vold(1)**2+vold(2)**2+vold(3)**2
*      vv_new2 = xdot(1,i)**2+xdot(2,i)**2+xdot(3,i)**2
*      ediss1 = ediss1 + 0.5*body(i)*(vv_new2-vv_old2)      
*         F(K,I) = F(K,I) - FDRAG(K)
*         FDOT(K,I) = FDOT(K,I) - FDDRAG(K)
*      end if
  10    RETURN
*
      END
      
