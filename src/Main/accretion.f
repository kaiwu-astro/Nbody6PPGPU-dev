       SUBROUTINE ACCRETION(I)
*
*
*       Accretion of stars onto central massive black hole.
*       Based on "star_destr_ext.c"
*       -------------------------
*       by Taras Panamarev
*
      INCLUDE 'common6.h'
      REAL*8 R_TIDAL0, R_TIDAL, E_kin, E_pot_bh, E_pot, E_corr_tmp,
     &       EPS, EPS_bh, EPS_bh2, RSB, RKB2, RKS2, R2, RR, NNB,
     &       X_IJ, Y_IJ, Z_IJ, VX_IJ, VY_IJ, VZ_IJ, V, TMP,
     &       A_EXT_I(3), ADOT_EXT_I(3), TEMP1, EPS_SS, EPS_SS2,
     &       DTR,DTR2HALF,DT_TMP,DT2HALF, DT3OVER6,RKS,XP(3),EPOT1
     &       ,t_diss_on, raccr, semi, E_tot, SEMIMAJOR, ECC, JSQR
     &       , rtmp, beta, m_ratio, rho_ratio, XVEC(3), VVEC(3)
     &       , v_kick, C_g, cos_theta, r_ratio
      real*4 b_rs1,b_rs2,b_l1,b_l2,b_te1,b_te2,b_rc1,b_rc2,b_m1,b_m2,
     &       b_acc, b_a, b_p,
     &       s_rs, s_l, s_te, s_rc, s_m, s_mc
*     &       ,FIRR(3) 
      INTEGER ii, L, KW, PERI_FLAG
*     & ,NXTLEN, NXTLST(NMAX)
*      LOGICAL is_active
*
*     Time when accretion gets 'switched on' 
      t_diss_on = 0.125
      t_diss_on = 5.000
*
      kw=-1
*     Update total time:
      TTOT = TIME + TOFF
*     
      if (ttot.lt.t_diss_on) goto 10     
*
      EPS_bh2 = EPS_bh**2
*
      E_pot = 0.0D0
      E_kin = 0.0D0
      E_pot_bh = 0.0D0
*
      R_TIDAL0 = 0.22*RACCR
*
      iaccr = 0
      K = I
*     Compute the real tidal radius according to current mass
      if(BODY(K)==0.0) return
*
      m_ratio = BODY(K)/(1.0/65536.0)
      if(m_ratio < 0.08) return
*     M-R relation of Adiabatic Model
*      call M_Rho_intplt(m_ratio,rho_ratio)
*      R_TIDAL = R_TIDAL0 * rho_ratio**(-1./3.)

*     M-R relation of Main-Sequence Model
      r_ratio = m_ratio**(0.8)
      R_TIDAL = R_TIDAL0 * r_ratio * m_ratio**(-1./3.)

*    Change x to x0, because x is predicted value
      RSB = SQRT(X0(1,K)**2+X0(2,K)**2+X0(3,K)**2)
*
      X_IJ = X0(1,K)
      Y_IJ = X0(2,K)
      Z_IJ = X0(3,K)
*     
      VX_IJ = X0DOT(1,K)
      VY_IJ = X0DOT(2,K)
      VZ_IJ = X0DOT(3,K)
*
      XVEC(1:3) = X0(1:3,K)
      VVEC(1:3) = X0DOT(1:3,K)
      call GET_JSQR(XVEC,VVEC,JSQR)
*
      V = SQRT(VX_IJ**2+VY_IJ**2+VZ_IJ**2)
*
      RV_IJ = VX_IJ*X_IJ+VY_IJ*Y_IJ+VZ_IJ*Z_IJ
      E_tot = 0.5*V*V - CMBH/RSB
      if( E_tot < 0.0 ) then
       SEMIMAJOR = -CMBH / (2.0*E_tot)
       ECC = sqrt(1.0-JSQR/(CMBH*SEMIMAJOR))
       beta = R_TIDAL/(SEMIMAJOR*(1.-ECC))
      else
       SEMIMAJOR = CMBH / (2.0*E_tot)
       ECC = sqrt(1.0+JSQR/(CMBH*SEMIMAJOR))
       beta = R_TIDAL/(SEMIMAJOR*(ECC-1.))
      end if

*      IF (RSB.LT.R_TIDAL.and.RSB.ne.0.0) THEN
*    Modified by Zhong
*    gamma=4/3 star case, full disruption
*    TDE only happens when star is approaching to BH
      IF((RSB.LT.(R_TIDAL/1.85)).and.(RSB.ne.0.0)
     &   .and.(RV_IJ/RSB/V.le.0.0)) THEN

            IF (K.GT.N) THEN 
*             binary is being accreted:
             KSPAIR = K - N
             
             call KSRES_OP(KSPAIR,J1,J2,XIREL,VIREL,0)

             if (rank.eq.0) PRINT*, 'A BINARY IS ACCRETED', TTOT,
     &        K-N,NAME(N), name(j1), name(j2)

*     Stellar evolution of two components
            IF (KZ(12).GT.0) THEN
               call sev_one_star_new(J1,kstar(j1),B_RS1,B_L1,
     &              B_TE1,B_MC1,B_RC1,B_M1)
               call sev_one_star_new(J2,kstar(j2),B_RS2,B_L2,
     &              B_TE2,B_MC2,B_RC2,B_M2)
            END IF
*     Binary parameters
            SEMI = -0.5*BODY(K)/H(KSPAIR)
            B_ECC = REAL(SQRT((1.0D0 - R(KSPAIR)/SEMI)**2+
     &           TDOT2(KSPAIR)**2/(BODY(K)*SEMI)))
            B_P = REAL(DAYS*SEMI*SQRT(ABS(SEMI)/BODY(I)))
            B_A = REAL(SEMI*RAU)
*     CONTINUE
      if (rank.eq.0) then
             OPEN (UNIT=28,STATUS='UNKNOWN',FORM='FORMATTED',
     &       FILE='accr_data_bin.dat', ACCESS='APPEND')
*
         WRITE (28, 96) TTOT,
     &    CMBH,NAME(J1),NAME(J2),NAME(K),BODY(J1),BODY(J2),
     &             (X0(J,K), J=1,3),(X0DOT(J,K), J=1,3),B_A,B_ECC,B_P,
     &             B_RS1, B_RS2, B_L1,B_L2,B_TE1,B_TE2,
     &             kstar(J1), kstar(J2), kstar(K) 
   96 FORMAT(f10.5,1x,e18.6,1x,3(I7,1x),17(1X,e18.6),1x,3(1x,I2))
         CALL FLUSH(28)
      endif
              IPHASE = 2
*
*            Check for rare case of merged binary component
             IF (NAME(K).LT.0) IPHASE = 7
*      
             GOTO 10
*         
            ENDIF
*       
         X_IJ = X0(1,K)
         Y_IJ = X0(2,K)
         Z_IJ = X0(3,K)
*         
         VX_IJ = X0DOT(1,K)
         VY_IJ = X0DOT(2,K)
         VZ_IJ = X0DOT(3,K)
*
         V = SQRT(VX_IJ**2+VY_IJ**2+VZ_IJ**2)         
*
         R2 = X_IJ**2+Y_IJ**2+Z_IJ**2
         RR = SQRT(R2)
*         
         RV_IJ = VX_IJ*X_IJ+VY_IJ*Y_IJ+VZ_IJ*Z_IJ
         TEMP1 = CMBH/(RR*R2)
*
*	 External force from BH acting on a particle
*     
         A_EXT_I(1) = -TeMP1*X_IJ
         A_EXT_I(2) = -TeMP1*Y_IJ
         A_EXT_I(3) = -TeMP1*Z_IJ
*
         ADOT_EXT_I(1) = -TeMP1*(VX_IJ-3.0*RV_IJ*X_IJ/R2)
         ADOT_EXT_I(2) = -TeMP1*(VY_IJ-3.0*RV_IJ*Y_IJ/R2)
         ADOT_EXT_I(3) = -TeMP1*(VZ_IJ-3.0*RV_IJ*Z_IJ/R2)
*            print *, "before sev_one_star:", s_rs, s_l, s_te, kw
                     call sev_one_star_new(K,kw,S_RS,S_L,
     &                    S_TE,S_MC,S_RC,S_M)
*            print *, "after sev_one_star:", s_rs, s_l, s_te, kw
      if (rank.eq.0) then
      OPEN (UNIT=28,STATUS='UNKNOWN',FORM='FORMATTED',
     &   FILE='accr_data.dat', ACCESS='APPEND')
*
         WRITE (28, 99) NAME(K), BODY(K), TTOT, CMBH,(X0(J,K), J=1,3),
     &             (X0DOT(J,K), J=1,3), ECC,beta,SEMIMAJOR, R_TIDAL
*     &             , S_RS, S_L, S_TE 
   99 FORMAT (I7, 13(1X, e24.16))
         CALL FLUSH(28)
      endif
                     
* !!!     PUT HERE STARDISC ON/OFF CRITERIUM
*
         E_KIN = 0.5*BODY(K)*(X0DOT(1,K)**2+X0DOT(2,K)**2+X0DOT(3,K)**2)
         E_POT_BH = -CMBH*BODY(K)/SQRT(RSB**2+EPS_bh2)
         NNB = LIST(1,K)
*      EPOT1 = 0.0D0
*      DO 20 J = 1,N
*      IF (J.EQ.K .OR. J.EQ.2*IPAIR-1 .OR. J.EQ.2*IPAIR .OR.
*     *    BODY(J).EQ.0.0D0 .OR. BODY(K).EQ.0.0D0)  GO TO 20
*          A1 = X(1,K) - X(1,J)
*          A2 = X(2,K) - X(2,J)
*          A3 = X(3,K) - X(3,J)
*      A4 = BODY(J)*BODY(K)/DSQRT(A1*A1 + A2*A2 + A3*A3)
*      EPOT1 = EPOT1 - A4
*
         TMP = 0.0D0
*
* !!! The body of do loop below needs to be corrected due to reg and ireg
* !!!      timestep - Done!
* !!!		
        DO 20 ii = IFIRST, NTOT
        IF  (ii.NE.K) THEN
*         RKS = SQRT(RKS2)
*	  IF  (RKS.LE.RS(K)) THEN
*      DO L = 1,NXTLEN
*         IF (ii.eq.NXTLST(L)) is_active = .true.
*      END DO   
*
*      IF (is_active) THEN
* 
            RKB2 = X0(1,ii)**2+X0(2,ii)**2+X0(3,ii)**2
            TMP = TMP - BODY(ii)/SQRT(RKB2+EPS_bh2)
*
      RKS2 = (X0(1,K)-X0(1,ii))**2+(X0(2,K)-X0(2,ii))**2+
     &   (X0(3,K)-X0(3,ii))**2
      TMP = TMP + BODY(ii)/SQRT(RKS2+EPS_ss2)
* !!! DON'T FORGET TO INCLUDE EPS_ss2 DEFINITION! - Done
*          ELSE
*            DT_TMP = TIME - T0(ii)
**            DTR = TIME - T0R(ii)
*            DT2HALF = 0.5*DT_TMP**2
**           DTR2HALF = 0.5*DTR**2
*            DT3OVER6 = ONE3*DT_TMP*DT2HALF
*            XP(1)=X(1,ii)+XDOT(1,ii)*DT_TMP
**     & +FIRR(1)*DT2HALF+FR(1,ii)*DTR2HALF+FDOT(1,ii)*DT3OVER6
*            XP(2)=X(2,ii)+XDOT(2,ii)*DT_TMP
**     & +FIRR(2)*DT2HALF+FR(2,ii)*DTR2HALF+ FDOT(2,ii)*DT3OVER6
*            XP(3)=X(3,ii)+XDOT(3,ii)*DT_TMP
**     & +FIRR(3)*DT2HALF+FR(3,ii)*DTR2HALF+FDOT(3,ii)*DT3OVER6

*          RKB2=XP(1)**2+XP(2)**2+XP(3)**2
*           TMP = TMP - BODY(ii)/SQRT(RKB2+EPS_bh2)
*
*           RKS2=(XP(1)-X(1,ii))**2+(XP(2)-X(2,ii))**2+
*     &          (XP(3)-X(3,ii))**2
*           TMP = TMP + BODY(ii)/SQRT(RKS2+EPS_ss2)          
*      End If is_active:
*          END IF
*!      End If ii:
          END IF
 20      CONTINUE 
*
*       !Accrete the mass of a star:
*       (Comment this out for Zhong Shiyan)
*        CMBH = CMBH + BODY(K)

***       NUM_ACCR = NUM_ACCR+1
*
*		!Create zero-mass particle and put the accreted star very far away from the cluster:
*        E_POT = EPOT1
       E_POT = -BODY(K)*TMP
        E_CORR_TMP = E_KIN + E_POT + E_POT_BH
        E_CORR = E_CORR + E_CORR_TMP
*
       if(rank.eq.0)PRINT*,'ACCR: RSB,NAME,TIME,EPOT,EKIN,EBH,ECORR,',
     & RSB, NAME(K), TTOT, E_POT, E_KIN, E_POT_BH, E_CORR  
* 
       BODY(K) = 0.0
*
*       t0(k) = time
*        T0(K) = TADJ + DTADJ
*        STEP(K) = 1.0D+06
*        STEPR(K) = 1.0D+06
*
        DO 30 L = 1,3
              X0(L,K) = 10000.0*RSCALE*X(L,K)/RR
              X(L,K) = X0(L,K)
              X0DOT(L,K) = SQRT(0.004*ZMASS/RSCALE)*XDOT(L,K)/V
              F(L,K) = 0.0D0
              FDOT(L,K) = 0.0D0
              D0(L,K) = 0.0
              D1(L,K) = 0.0
              D2(L,K) = 0.0
              D3(L,K) = 0.0
              D0R(L,K) = 0.0
              D1R(L,K) = 0.0
              D2R(L,K) = 0.0
              D3R(L,K) = 0.0
 30    CONTINUE
*
*      
        call delay_remove_tlist(K,step,dtk)
        step(k)=2.0*DTK(1)
        TEV(k)=1E6

*     !Remove accreted star from neighbour lists
*	NNB = LIST(1,K)		
        CALL NBREM(K,1,NNB)
        LIST(1,K) = 0
*
*		Set IPHASE < 0 to ensure updating of time-step sequence
*
        IPHASE = -1
        IACCR = 1
*
       if(rank.eq.0) PRINT*, 'BHDATA :',  TTOT, CMBH
*      50 FORMAT ('BHDATA: ',F8.2,E15.6)

*
*      GOTO 10
*
* 40   CONTINUE
*  ! end if r_tidal:
       END IF
*      !End DO K=IFIRST,NTOT

*      IF(RSB>R_TIDAL.and.RSB<(2.0*R_TIDAL)) THEN
*    Modified by Zhong
*    gamma=4/3 star, PTDE zone
*      IF(RSB>(0.4*R_TIDAL0).and.RSB<(2.0*R_TIDAL0)) THEN
*        if (rank.eq.0) then
*          OPEN (UNIT=48,STATUS='UNKNOWN',FORM='FORMATTED',
*     &    FILE='pTDE_data.dat', ACCESS='APPEND')
*
*          WRITE(48, 97) NAME(K), BODY(K), TTOT, CMBH,(X0(J,K), J=1,3),
*     &             (X0DOT(J,K), J=1,3), (F(J,K),J=1,3), kstar(K)
*   97   FORMAT(I7,1X,e18.6,1X,e24.16,10(1X, e18.6),1x,i2)
*          CALL FLUSH(48)
*          CLOSE(48)
*        endif
*      END IF
*
      IF(RSB.ge.(R_TIDAL/1.85).and.RSB<(R_TIDAL/0.6)) THEN
*     -- 200916 add by Zhong
*     -- Compute inner product of R & V, compare with RVSIGN
*        and determine PERI_FLAG
        rtmp = RV_IJ / RSB / V
        PERI_FLAG = 0
        if( rtmp.ge.0.0 .and. RVSIGN(NAME(K)) < 0.0 ) PERI_FLAG=1
        RVSIGN(NAME(K)) = rtmp
*     -- if PERI_FLAG==1, write to file
        if (rank.eq.0 .and. PERI_FLAG==1) then
          OPEN (UNIT=48,STATUS='UNKNOWN',FORM='FORMATTED',
     &    FILE='pTDE_data_peri.dat', ACCESS='APPEND')
*
          WRITE(48, 98) NAME(K), BODY(K), TTOT, CMBH,(X0(J,K), J=1,3),
     &             (X0DOT(J,K), J=1,3), ECC,beta,SEMIMAJOR, R_TIDAL
   98   FORMAT(I7,1X,e18.6,1X,e24.16,11(1X, e24.16))
          CALL FLUSH(48)
          CLOSE(48)
        endif

*     -- 200916 end by Zhong
      END IF

*     do the kick due to PTDE, 200922 add by Zhong
      IF(RSB.ge.(R_TIDAL/1.85).and.RSB<(R_TIDAL/0.6)
     &  .and.PERI_FLAG==1)THEN
*        call GET_V_KICK_Dai(v_kick,m_ratio,rho_ratio,beta)
        call GET_V_KICK_MS(v_kick,m_ratio,r_ratio,beta)
*     Because this position might not be the real pericenter,
*     we need to adjust the value of v_kick according to current
*     position and velocity
        cos_theta=RVSIGN(NAME(K))
        v_kick=-V*cos_theta+sqrt((V*cos_theta)**2.+v_kick**2.)
        X0DOT(1,K)=X0DOT(1,K)+v_kick*X0(1,K)/RSB
        X0DOT(2,K)=X0DOT(2,K)+v_kick*X0(2,K)/RSB
        X0DOT(3,K)=X0DOT(3,K)+v_kick*X0(3,K)/RSB
*     Modify stellar mass
        call C_4_3(C_g,beta)
        if(BODY(K)*(1.-C_g)>0.0) then
          BODY(K)=BODY(K)*(1.-C_g)
        else
          BODY(K)=0.0
        end if

*     Compute the potential generated by star cluster at the position
*     of the disruptee
        TMP = 0.0
        DO ii = IFIRST, NTOT
          IF  (ii.NE.K) THEN
            RKS2 = (X0(1,K)-X0(1,ii))**2+(X0(2,K)-X0(2,ii))**2
     &       +(X0(3,K)-X0(3,ii))**2
            TMP = TMP - BODY(ii)/SQRT(RKS2)
          END IF
        END DO

        if (rank.eq.0) then
*          print*,NAME(K),ECC,beta,SEMIMAJOR,JSQR/(CMBH*SEMIMAJOR)
          OPEN (UNIT=48,STATUS='UNKNOWN',FORM='FORMATTED',
     &    FILE='pTDE_data_kick.dat', ACCESS='APPEND')
*
          WRITE(48, 90) NAME(K), BODY(K), TTOT, CMBH,(X0(J,K), J=1,3),
     &             (X0DOT(J,K), J=1,3), ECC,beta,SEMIMAJOR, TMP
   90   FORMAT(I7,1X,e18.6,1X,e24.16,11(1X, e24.16))
          CALL FLUSH(48)
          CLOSE(48)
        endif

      END IF
*     End of kick stuff, 200922 add by Zhong

*
10      CONTINUE

      RETURN
*
      END
*----------------------------------------------------------------------------
      subroutine sev_one_star_new(I,KW,RM,LUM,TE,MC,RCC,M1)
*
*
*     Get stellar evolution parameters for one star
*     ---------------------------------------------
*
      include 'common6.h'
      INTEGER I,KW
      REAL*8  LUMS(10),TSCLS(20),GB(10)
      REAL*8  RM8,LUM8,MC8,RCC8,M18,M0
      REAL*4  RM,LUM,TE,MC,RCC,M1
      REAL*8 RSCALE_OUT,MSCALE_OUT,VSCALE_OUT,RAU_OUT,TSCALE_OUT

      IF (KZ(19).EQ.0.AND.KZ(12).EQ.-1) THEN
         RSCALE_OUT=1.0
         MSCALE_OUT=1.0
         VSCALE_OUT=1.0
         RAU_OUT=1.0
         TSCALE_OUT=1.0
      else
         RSCALE_OUT=RBAR
         MSCALE_OUT=ZMBAR
         VSCALE_OUT=VSTAR
         RAU_OUT=RAU
         TSCALE_OUT=TSTAR
      END IF

*
      KW = KSTAR(I)
      AGE = MAX(TIME,TEV0(I))*TSCALE_OUT - EPOCH(I)
      M0 = BODY0(I)*MSCALE_OUT
      M18 = BODY(I)*MSCALE_OUT
      CALL STAR(KW,M0,M18,TM,TN,TSCLS,LUMS,GB,ZPARS)
      CALL HRDIAG(M0,AGE,M18,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &     RM8,LUM8,KW,MC8,RCC8,ME,RE,K2)
*     Temperature (Solar temperature got from Williams, D. R. (1 July 2013). "Sun Fact Sheet". NASA. Retrieved 12 August 2013.)
      TE = REAL(5778*(LUM8/(RM8*RM8))**0.25)

      RM = REAL(RM8)
      LUM = REAL(LUM8)
      MC = REAL(MC8)
      RCC = REAL(RCC8)
      M1 = REAL(M18)

      RETURN

      END
*     200922 Add by Zhong
*     *********************************************
*     Compute the coefficient C_4_3
      subroutine C_4_3(C_g,beta)
      REAL*8  C_g, beta, a, b

      if(beta < 0.6) then
        C_g = 0.0
      else if(beta<=1.85) then
        a = 12.996 - 31.149*beta + 12.865*beta**2.0
        b = 1.0000 - 5.3232*beta + 6.4262*beta**2.0
        C_g = exp(a/b)
      else
        C_g = 1.00
      end if
      RETURN

      END

*     *********************************************
*     Compute the mean density of star with mass M
*     Rho is used to compute tidal radius
*     Note: M is actually M/M0
*     Note: Rho is actually Rho/Rho0
      subroutine M_rho_intplt(M,rho)
      INTEGER I
      real*8 x, y, k, rho, M, tab(2,11)
      tab(1,1) = 0.1200
      tab(1,2) = 0.3200
      tab(1,3) = 0.4288
      tab(1,4) = 0.5136
      tab(1,5) = 0.5840
      tab(1,6) = 0.6496
      tab(1,7) = 0.7136
      tab(1,8) = 0.7760
      tab(1,9) = 0.8384
      tab(1,10) = 0.9120
      tab(1,11) = 1.0000

      tab(2,1) = 0.6229
      tab(2,2) = 2.2546
      tab(2,3) = 2.9798
      tab(2,4) = 3.4590
      tab(2,5) = 3.7698
      tab(2,6) = 3.8863
      tab(2,7) = 3.7309
      tab(2,8) = 3.4849
      tab(2,9) = 3.0187
      tab(2,10) = 2.1640
      tab(2,11) = 1.0000

      if( M == 1.0000) then
        Rho = 1.0000
        return
      end if
      if( M < tab(1,1) ) then
        i = 1
        k = (tab(2,i+1) - tab(2,i))/(tab(1,i+1) - tab(1,i))
      else
        do i=1,10
         if( M>=tab(1,i) .and. M<tab(1,i+1)) exit
        end do
        k = (tab(2,i+1) - tab(2,i))/(tab(1,i+1) - tab(1,i))
      end if
      rho = tab(2,i)+k*(M-tab(1,i))
      RETURN
      END

*     *********************************************
*     Compute the kick velocity for Adiabatic Model
*     Note: M is actually M/M0
*     Note: Rho is actually Rho/Rho0
      subroutine GET_V_KICK_Dai(v_kick,M,Rho,beta)
      real*8 beta, v_kick, M, Rho
      real*8 v_esc, v_esc0
*     convert from km/s to model unit
      v_esc0 = 617.0 / 106.9
*     The real v_esc varies with stellar mass and radius
      v_esc = v_esc0 * M**(1./3.) * Rho**(1./6.)
      v_kick = v_esc * (0.0745 + 0.0571*beta**4.539)
      return
      end

*     *********************************************
*     Compute the kick velocity for MS Model
*     Note: M is actually M/M0
*     Note: Rho is actually Rho/Rho0
      subroutine GET_V_KICK_MS(v_kick,M,r,beta)
      real*8 beta, v_kick, M, r
      real*8 v_esc, v_esc0
*     convert from km/s to model unit
      v_esc0 = 617.0 / 106.9
*     The real v_esc varies with stellar mass and radius
      v_esc = v_esc0 * M**(1./2.) * r**(-1./2.)
      v_kick = v_esc * (0.0745 + 0.0571*beta**4.539)
      return
      end

*     *********************************************
*     Compute the square of angular momentum
      subroutine GET_JSQR(x,v,JSQR)
      real*8 x(3), v(3), JSQR, JVEC(3)
      JVEC(1) = x(2)*v(3)-x(3)*v(2)
      JVEC(2) = x(3)*v(1)-x(1)*v(3)
      JVEC(3) = x(1)*v(2)-x(2)*v(1)
      JSQR=JVEC(1)**2.0+JVEC(2)**2.0+JVEC(3)**2.0
      return
      end
