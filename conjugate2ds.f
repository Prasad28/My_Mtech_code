C#############################################################
      PROGRAM CONJUGATE2DS_MULTI
C#############################################################
      INCLUDE 'com.inc'
      LOGICAL LWRITE,LREAD,LCAL,LAXIS,LTIME
      CHARACTER*20 FILRES,FILIN,FILGR
C--------------------------------------------------------------
C
      FILIN='conjugate2ds.inp'
      FILGR='filgr'
      FILRES='res'  
C
      OPEN (UNIT=5,FILE=FILIN)
      OPEN (UNIT=1,FILE=FILGR,FORM='UNFORMATTED')
      OPEN (UNIT=3,FILE=FILRES,FORM='UNFORMATTED')
      REWIND 3
      REWIND 5
      REWIND 1
C
C.....INPUT AND BOUNDARY DATA, INITIALIZATION, OUTPUT TITLE, ETC.
C
      CALL MODINP
C
      IF(LREAD)THEN
        OPEN(UNIT=7,FILE='residual1')
        REWIND 7
       ELSE
        OPEN(UNIT=7,FILE='residual')
        REWIND 7  
      ENDIF
C
C==============================================
C.....TIME LOOP
C==============================================
C
 160  ITIM=0
      TIME=0.0D0
      ITIMS=ITIM+1
      ITIME=ITIM+ITST
C
C
      DO ITIM=ITIMS,ITIME
        TIME=TIME+DT
C
C.....SHIFT SOLUTIONS IN TIME
C
      IF(LTIME) THEN
        DO IJ=1,NIJS
          TOO(IJ)=TO(IJ)
          TO(IJ)=T(IJ)
          UOO(IJ)=UO(IJ)
          UO(IJ)=U(IJ)
          VOO(IJ)=VO(IJ)
          VO(IJ)=V(IJ)
        END DO
      ENDIF
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++
C.....OUTER ITERATIONS (SIMPLE RELAXATIONS)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DO ITER=1,MAXIT
        IF(LCAL(IU))CALL CALCUV
        IF(LCAL(IP))CALL CALCP
        IF(LCAL(IEN))CALL CALCSC
C
C.....CHECK CONVERGENCE OF OUTER ITERATIONS
        SOURCE=MAX(RESOR(IU),RESOR(IV),RESOR(IP),RESOR(IEN))
        PRINT*,'ITIM =',ITIM
        PRINT*,'ITERATION =',ITER
        PRINT*,'RESOR(IU) =',RESOR(IU)
        PRINT*,'RESOR(IV) =',RESOR(IV)
        PRINT*,'RESOR(IP) =',RESOR(IP)
        PRINT*,'RESOR(IEN) =',RESOR(IEN)
c
        IF(SOURCE.GT.SLARGE) GO TO 510
        IF(SOURCE.LT.SORMAX) GO TO 250
      ENDDO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
  250 CONTINUE
C
      END DO
C
      WRITE(3)(F1(IJ),IJ=1,NIJS),(F2(IJ),IJ=1,NIJS),(U(IJ),IJ=1,NIJS),
     *        (V(IJ),IJ=1,NIJS),(P(IJ),IJ=1,NIJS),(T(IJ),IJ=1,NIJS),
     *       (UO(IJ),IJ=1,NIJS),(VO(IJ),IJ=1,NIJS),(TO(IJ),IJ=1,NIJS)
     *       ,(DTX(IJ),IJ=1,NIJS),(DTY(IJ),IJ=1,NIJS)
      CLOSE(3)
C

 150  STOP
C
C==============================================================
C......MESSAGE FOR DIVERGENCE 
C==============================================================
C
  510 PRINT *,'  *** TERMINATED - OUTER ITERATIONS DIVERGING ***'
      WRITE(11,*)'  *** TERMINATED - OUTER ITERATIONS DIVERGING ***'
      STOP 
C
  710 FORMAT(1X,I5,1X,I5,2X,1P4E12.4)
C
      END
C
C###########################################################
      SUBROUTINE CALCUV
C###########################################################
      INCLUDE 'com.inc'

      LOGICAL LWRITE,LREAD,LAXIS,LCAL,LTIME

C----------------------------------------------------------
C
C.....RECIPROCAL VALUES OF UNDER-RELAXATION FACTORS FOR U AND V
C
      URFU=1.0D0/URF(IU)
      URFV=1.0D0/URF(IV)
C
      CALL GRADFIUVP(U,DUX,DUY)
      CALL GRADFIUVP(V,DVX,DVY)
      CALL GRADFIUVP(P,DPX,DPY)
C.....INITIALIZE TEMPORARILY STORED VARIABLES
C
      DO IJ=1,NIJS
        SU(IJ)=0.0D0
        SV(IJ)=0.0D0
        APU(IJ)=0.0D0
        APV(IJ)=0.0D0
      END DO
C
C==========================================================
C.....FLUXES THROUGH INTERNAL EAST CV FACES 
C==========================================================
      GDS1=GDS(IU)
      DO IB=1,1
      DO I=2,NIM(IB)-1
         DO J=2,NJM(IB)
         IJ=LB(IB)+LI(IB,I)+J
         CALL FLUXUV(IJ,IJ+NJ(IB),IJ,IJ-1,F1(IJ),AW(IJ+NJ(IB)),AE(IJ),
     *FX(IJ),GDS(IU),SEX(IJ),SEY(IJ))
      END DO
      END DO
      ENDDO  
C=========================================================
C.....FLUXES THROUGH INTERNAL NORTH CV FACES 
C=========================================================
      DO IB=1,1  
      DO I=2,NIM(IB)
         DO J=2,NJM(IB)-1
           IJ=LB(IB)+LI(IB,I)+J
      CALL FLUXUV(IJ,IJ+1,IJ-NJ(IB),IJ,F2(IJ),AS(IJ+1),AN(IJ),
     *FY(IJ),GDS(IV),SNX(IJ),SNY(IJ))
      END DO
      END DO
      ENDDO
C
C=============================================================
C.....VOLUME INTEGRALS (SOURCE TERMS)
C=============================================================
      DO IB=1,1
      DO I=2,NIM(IB)
C
        DO J=2,NJM(IB)
          IJ=LB(IB)+LI(IB,I)+J
          SU(IJ)=SU(IJ)-DPX(IJ)*VOL(IJ)
          SV(IJ)=SV(IJ)-DPY(IJ)*VOL(IJ)
C
C..... BUOYANCY SOURCE CONTRIBUTION
C
          IF(LCAL(IEN)) THEN
            SBB=-DEN(IJ)*BETA*VOL(IJ)*(T(IJ)-TREF)
            SU(IJ)=SU(IJ)+GRAVX*SBB
            SV(IJ)=SV(IJ)+GRAVY*SBB
          ENDIF
C
C..... UNSTEADY TERM CONTRIBUTION TO AP AND SU
C
          IF(LTIME) THEN
            APT=DEN(IJ)*VOL(IJ)*DTR
       SU(IJ)=SU(IJ)+(1.0D0+GAMT)*APT*UO(IJ)-0.5D0*GAMT*APT*UOO(IJ)
       SV(IJ)=SV(IJ)+(1.0D0+GAMT)*APT*VO(IJ)-0.5D0*GAMT*APT*VOO(IJ)
            APV(IJ)=APV(IJ)+(1.0D0+0.5D0*GAMT)*APT
            APU(IJ)=APU(IJ)+(1.0D0+0.5D0*GAMT)*APT
          ENDIF
C
        END DO
      END DO
      ENDDO
C
C.....AXISYMMETRIC CONTRIBUTION
C
      IF(LAXIS) THEN
       DO IB=1,1
        DO I=2,NIM(IB)
        DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)
          RCR=4.0D0/(X(IJ)+X(IJ-1)+X(IJ-NJ(IB))+X(IJ-NJ(IB)-1))
          APV(IJ)=APV(IJ)+2.0D0*VIS(IJ)*VOL(IJ)*RCR**2
        END DO
        END DO
       ENDDO 
      ENDIF
C
      CALL BCUV
C
C=============================================================
C.....UNDER-RELAXATION, SOLVING EQUATION SYSTEM FOR U-VELOCITY
C=============================================================
C
      DO IB=1,1
      DO I=2,NIM(IB)
      DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)
        AP(IJ)=(-AE(IJ)-AW(IJ)-AN(IJ)-AS(IJ)+APU(IJ))*URFU
        SU(IJ)=SU(IJ)+(1.0D0-URF(IU))*AP(IJ)*U(IJ)
        APU(IJ)=1.0D0/(AP(IJ)+1.0E-20)
      END DO
      END DO
      ENDDO
C
      CALL SIPSOLUVP(U,IU,1)
C
C=============================================================
C.....UNDER-RELAXATION, SOLVING EQUATION SYSTEM FOR V-VELOCITY
C=============================================================
C
      DO IB=1,1
      DO I=2,NIM(IB)
      DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)
        AP(IJ)=(-AE(IJ)-AW(IJ)-AN(IJ)-AS(IJ)+APV(IJ))*URFV
      SU(IJ)=SV(IJ)+(1.0D0-URF(IV))*AP(IJ)*V(IJ)
        APV(IJ)=1.0D0/(AP(IJ)+1.0E-20)
      END DO
      END DO
      ENDDO
C
      CALL SIPSOLUVP(V,IV,1)
C
      RETURN
      END 
C
C
C###########################################################
      SUBROUTINE BCUV
C###########################################################
      INCLUDE 'com.inc'

      DO IB=1,1
C
C.....BOUNDARY CONDITION OF BLOCK 1
C
      IF(IB.EQ.1)THEN
C.....SOUTH (WALL)
C
      DO I=2,NIM(IB)
         IJ=LB(IB)+LI(IB,I)+2
         IJB=IJ-1
         SX=SSX(IJ)
         SY=SSY(IJ)
         ARE=SQRT(SX**2+SY**2)
         XNN=SX/(ARE+1.0D-20)
         YNN=SY/(ARE+1.0D-20)
         XTW=-YNN
         YTW=XNN
         DN=(XC(IJB)-XC(IJ))*XNN+(YC(IJB)-YC(IJ))*YNN
         SRDW=ARE/(DN+1.0E-20)
         COEF=VIS(IJ)*SRDW
         APU(IJ) =APU(IJ) +COEF*XTW**2
         APV(IJ) =APV(IJ)+COEF*YTW**2
         SU(IJ)=SU(IJ)+COEF*(U(IJB)*XTW**2-
     *         (V(IJ)-V(IJB))*XTW*YTW)
         SV(IJ)=SV(IJ)+COEF*(V(IJB)*YTW**2-
     *         (U(IJ)-U(IJB))*XTW*YTW)
      END DO
C
C.....NORTH BOUNDARY (WALL)
C
      DO I=2,NIM(IB)
         IJ=LB(IB)+LI(IB,I)+NJM(IB)
         IJB=IJ+1
         SX=SNX(IJ)
         SY=SNY(IJ)
         ARE=SQRT(SX**2+SY**2)
         XNN=SX/(ARE+1.0D-20)
         YNN=SY/(ARE+1.0D-20)
         XTW=-YNN
         YTW=XNN
         DN=(XC(IJB)-XC(IJ))*XNN+(YC(IJB)-YC(IJ))*YNN
         SRDW=ARE/(DN+1.0E-20)
         COEF=VIS(IJ)*SRDW
         APU(IJ) =APU(IJ) +COEF*XTW**2
           APV(IJ) =APV(IJ)+COEF*YTW**2
           SU(IJ)=SU(IJ)+COEF*(U(IJB)*XTW**2-
     *         (V(IJ)-V(IJB))*XTW*YTW)
           SV(IJ)=SV(IJ)+COEF*(V(IJB)*YTW**2-
     *         (U(IJ)-U(IJB))*XTW*YTW)
      END DO
C
C.....WEST BOUNDARY (WALL)
C
      DO J=2,NJM(IB)
         IJ=LB(IB)+LI(IB,2)+J
         IJB=IJ-NJ(IB)
         SX=SWX(IJ)
         SY=SWY(IJ)
         ARE=SQRT(SX**2+SY**2)
         XNN=SX/(ARE+1.0D-20)
         YNN=SY/(ARE+1.0D-20)
         XTW=-YNN
         YTW=XNN
         DN=(XC(IJB)-XC(IJ))*XNN+(YC(IJB)-YC(IJ))*YNN
         SRDW=ARE/(DN+1.E-20)
         COEF=VIS(IJ)*SRDW
         APU(IJ) =APU(IJ) +COEF*XTW**2
         APV(IJ) =APV(IJ)+COEF*YTW**2
         SU(IJ)=SU(IJ)+COEF*(U(IJB)*XTW**2-
     *         (V(IJ)-V(IJB))*XTW*YTW)
         SV(IJ)=SV(IJ)+COEF*(V(IJB)*YTW**2-
     *         (U(IJ)-U(IJB))*XTW*YTW)
      END DO
C
C.....EAST BOUNDARY (WALL)
C
      DO J=2,NJM(IB)
         IJ=LB(IB)+LI(IB,NIM(IB))+J
         IJB=IJ+NJ(IB)
         SX=SEX(IJ)
         SY=SEY(IJ)
         ARE=SQRT(SX**2+SY**2)
         XNN=SX/(ARE+1.0D-20)
         YNN=SY/(ARE+1.0D-20)
         XTW=-YNN
         YTW=XNN
         DN=(XC(IJB)-XC(IJ))*XNN+(YC(IJB)-YC(IJ))*YNN
         SRDW=ARE/(DN+1.E-20)
         COEF=VIS(IJ)*SRDW
         APU(IJ) =APU(IJ) +COEF*XTW**2
         APV(IJ) =APV(IJ)+COEF*YTW**2
         SU(IJ)=SU(IJ)+COEF*(U(IJB)*XTW**2-
     *         (V(IJ)-V(IJB))*XTW*YTW)
         SV(IJ)=SV(IJ)+COEF*(V(IJB)*YTW**2-
     *         (U(IJ)-U(IJB))*XTW*YTW)
      END DO
C
      ENDIF  
C.....END OF BLOCK 1
C
      ENDDO
CXXXXXXXXX
      RETURN
      END
C
C################################################################
      SUBROUTINE FLUXUV(IJP,IJN,IJ1,IJ2,FM,CAP,CAN,FAC,G,SX,SY)
C################################################################
C==============================================================
      INCLUDE 'com.inc'
C
C.....INTERPOLATED CELL FACE VALUES
C
      FACP=1.0D0-FAC
C   
      XI=XC(IJN)*FAC+XC(IJP)*FACP
      YI=YC(IJN)*FAC+YC(IJP)*FACP
      XF=0.5D0*(X(IJ1)+X(IJ2))
      YF=0.5D0*(Y(IJ1)+Y(IJ2))
C
      DUXI=DUX(IJN)*FAC+DUX(IJP)*FACP
      DVXI=DVX(IJN)*FAC+DVX(IJP)*FACP
      DUYI=DUY(IJN)*FAC+DUY(IJP)*FACP
      DVYI=DVY(IJN)*FAC+DVY(IJP)*FACP
      UI=U(IJN)*FAC+U(IJP)*FACP+DUXI*(XF-XI)+DUYI*(YF-YI)
      VI=V(IJN)*FAC+V(IJP)*FACP+DVXI*(XF-XI)+DVYI*(YF-YI)
C
C.....SURFACE AND DISTANCE VECTOR COMPONENTS, DIFFUSION COEFFICIENT
C
      VISI=VIS(IJN)*FAC+VIS(IJP)*FACP
      XPN=XC(IJN)-XC(IJP)
      YPN=YC(IJN)-YC(IJP)
      VSOL=VISI*SQRT((SX**2+SY**2)/(XPN**2+YPN**2))
C
C.....EXPLICIT CONVECTIVE AND DIFFUSIVE FLUXES
C
      FCUE=FM*UI
      FCVE=FM*VI
      FDUE=VISI*(2.0D0*DUXI*SX+(DUYI+DVXI)*SY)
      FDVE=VISI*((DUYI+DVXI)*SX+2.0D0*DVYI*SY)
C
C.....IMPLICIT CONVECTIVE AND DIFFUSIVE FLUXES
C
      FCUI=MIN(FM,0.0D0)*U(IJN)+MAX(FM,0.0D0)*U(IJP)
      FCVI=MIN(FM,0.0D0)*V(IJN)+MAX(FM,0.0D0)*V(IJP)
      FDUI=VSOL*(DUXI*XPN+DUYI*YPN)
      FDVI=VSOL*(DVXI*XPN+DVYI*YPN)
C
C.....COEFFICIENTS, DEFERRED CORRECTION, SOURCE TERMS
C
      CAN=-VSOL+MIN(FM,0.0D0)
      CAP=-VSOL-MAX(FM,0.0D0)
      FUC=G*(FCUE-FCUI)
      FVC=G*(FCVE-FCVI)
C
      SU(IJP)=SU(IJP)-FUC+(FDUE-FDUI)
      SU(IJN)=SU(IJN)+FUC+(-FDUE+FDUI)
      SV(IJP)=SV(IJP)-FVC+(FDVE-FDVI)
      SV(IJN)=SV(IJN)+FVC+(-FDVE+FDVI)
C
      RETURN
      END
C
C############################################################## 
      SUBROUTINE CALCP 
C##############################################################
C-------------------------------------------------------------- 
      INCLUDE 'com.inc'

      LOGICAL LWRITE,LREAD,LAXIS,LCAL,LTIME
C
C-------------------------------------------------------------- 
C
      DO IB=1,1
      DO I=2,NIM(IB)
         DO J=2,NJM(IB)
            IJ=LB(IB)+LI(IB,I)+J
            SU(IJ)=0.0D0
            AP(IJ)=0.0D0
         ENDDO
      ENDDO
      ENDDO
C============================================================
C.....EAST CV FACES (S - AREA, VOLE - VOLUME BETWEEN P AND E)
C============================================================
C
      DO IB=1,1
      DO I=2,NIM(IB)-1
      DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)
        CALL FLUXME(IJ,IJ+NJ(IB),IJ,IJ-1,F1(IJ),AW(IJ+NJ(IB)),
     *AE(IJ),FX(IJ),SEX(IJ),SEY(IJ))
      END DO
      END DO
      ENDDO
C
C=============================================================
C.....NORTH CV FACES (S - AREA, VOLN - VOLUME BETWEEN P AND N)
C=============================================================
C
      DO IB=1,1
      DO I=2,NIM(IB)
         DO J=2,NJM(IB)-1
        IJ=LB(IB)+LI(IB,I)+J
        CALL FLUXMN(IJ,IJ+1,IJ-NJ(IB),IJ,F2(IJ),AS(IJ+1),AN(IJ),FY(IJ),
     *SNX(IJ),SNY(IJ))
         ENDDO
      END DO
      END DO
C
C===============================================================
C..... SORCE TERM AND COEFFICIENT OF NODE P
C===============================================================
C
      SUM=0.0D0
      DO IB=1,1  
      DO I=2,NIM(IB)
      DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)
        PP(IJ)=0.0D0
        SU(IJ)=SU(IJ)+F1(IJ-NJ(IB))-F1(IJ)+F2(IJ-1)-F2(IJ)
        AP(IJ)=AP(IJ)-(AE(IJ)+AW(IJ)+AN(IJ)+AS(IJ))
        SUM=SUM+SU(IJ)
      END DO
      END DO
      ENDDO
C
C===============================================================
C.....SOLVE EQUATIONS SYSTEM FOR P' AND APPLY CORRECTIONS
C===============================================================
C
      CALL SIPSOLUVP(PP,IP,1)
C
C.....CALCULATE PRESSURE CORRECTION AT BOUNDARIES
C
      CALL PBOUND(PP)

C.....VALUE OF P' AT REFERENCE LOCATION TO BE SUBTRACTED FROM ALL P'
C
      CALL GRADFIUVP(PP,DPX,DPY)
      IJPREF=LB(1)+LI(1,IPR)+JPR
      PPO=PP(IJPREF)
C
C.....CORRECT EAST MASS FLUXES 
C
      DO IB=1,1
      DO I=2,NIM(IB)-1
      DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)
        F1(IJ)=F1(IJ)+AE(IJ)*(PP(IJ+NJ(IB))-PP(IJ))
      END DO
      END DO
      ENDDO
C
C.....CORRECT NORTH MASS FLUXES 
C
      DO IB=1,1
      DO I=2,NIM(IB)
      DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)-1
        F2(IJ)=F2(IJ)+AN(IJ)*(PP(IJ+1)-PP(IJ))
      END DO
      END DO
      ENDDO  
C
C.....CORRECT PRESSURE AND VELOCITIES AT CELL CENTER
C
      DO IB=1,1
      DO I=2,NIM(IB)
        DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)
          U(IJ)=U(IJ)-DPX(IJ)*VOL(IJ)*APU(IJ)
          V(IJ)=V(IJ)-DPY(IJ)*VOL(IJ)*APV(IJ)
          P(IJ)=P(IJ)+URF(IP)*(PP(IJ)-PPO)
          SU(IJ)=0.0D0
        END DO
      END DO
      ENDDO
C
      CALL PBOUND(P)
C
      RETURN
      END
C
C##############################################################
      SUBROUTINE FLUXME(IJP,IJN,IJ1,IJ2,FM,CAP,CAN,FAC,SX,SY)
C##############################################################
      INCLUDE 'com.inc'
C
C.....INTERPOLATED CELL FACE VALUES
C
      FACP=1.0D0-FAC
      XI=XC(IJN)*FAC+XC(IJP)*FACP
      YI=YC(IJN)*FAC+YC(IJP)*FACP
      XF=0.5D0*(X(IJ1)+X(IJ2))
      YF=0.5D0*(Y(IJ1)+Y(IJ2))
C
      DUXI=DUX(IJN)*FAC+DUX(IJP)*FACP
      DVXI=DVX(IJN)*FAC+DVX(IJP)*FACP
      DUYI=DUY(IJN)*FAC+DUY(IJP)*FACP
      DVYI=DVY(IJN)*FAC+DVY(IJP)*FACP
      UI=U(IJN)*FAC+U(IJP)*FACP+DUXI*(XF-XI)+DUYI*(YF-YI)
      VI=V(IJN)*FAC+V(IJP)*FACP+DVXI*(XF-XI)+DVYI*(YF-YI)
C
      DENI=DEN(IJN)*FAC+DEN(IJP)*FACP
C
C.....SURFACE AND DISTANCE VECTOR COMPONENTS
C
      XPN=XC(IJN)-XC(IJP)
      YPN=YC(IJN)-YC(IJP)
      SMDPN=(SX**2+SY**2)/(SX*XPN+SY*YPN)   
C
C.....MASS FLUX, COEFFICIENTS FOR THE P'-EQUATION
C
      CAP=-0.5D0*(APU(IJP)*DEN(IJP)+APU(IJN)*DEN(IJN))*(SX**2+SY**2)
      CAN=CAP

      DPXI=0.5D0*(DPX(IJN)+DPX(IJP))*XPN
      DPYI=0.5D0*(DPY(IJN)+DPY(IJP))*YPN
      FM=DENI*(UI*SX+VI*SY)+CAP*(P(IJN)-P(IJP)-DPXI-DPYI)
C
      RETURN
      END
C
C##############################################################
      SUBROUTINE FLUXMN(IJP,IJN,IJ1,IJ2,FM,CAP,CAN,FAC,SX,SY)
C##############################################################
      INCLUDE 'com.inc'
C
C.....INTERPOLATED CELL FACE VALUES
C
      FACP=1.0D0-FAC
      XI=XC(IJN)*FAC+XC(IJP)*FACP
      YI=YC(IJN)*FAC+YC(IJP)*FACP
      XF=0.5D0*(X(IJ1)+X(IJ2))
      YF=0.5D0*(Y(IJ1)+Y(IJ2))
C
      DUXI=DUX(IJN)*FAC+DUX(IJP)*FACP
      DVXI=DVX(IJN)*FAC+DVX(IJP)*FACP
      DUYI=DUY(IJN)*FAC+DUY(IJP)*FACP
      DVYI=DVY(IJN)*FAC+DVY(IJP)*FACP
      UI=U(IJN)*FAC+U(IJP)*FACP+DUXI*(XF-XI)+DUYI*(YF-YI)
      VI=V(IJN)*FAC+V(IJP)*FACP+DVXI*(XF-XI)+DVYI*(YF-YI)
C
      DENI=DEN(IJN)*FAC+DEN(IJP)*FACP
C
C.....SURFACE AND DISTANCE VECTOR COMPONENTS
C
      XPN=XC(IJN)-XC(IJP)
      YPN=YC(IJN)-YC(IJP)
      SMDPN=(SX**2+SY**2)/(SX*XPN+SY*YPN)   
C
C.....MASS FLUX, COEFFICIENTS FOR THE P'-EQUATION
C
      CAP=-0.5D0*(APV(IJP)*DEN(IJP)+APV(IJN)*DEN(IJN))*(SX**2+SY**2)
      CAN=CAP

      DPXI=0.5D0*(DPX(IJN)+DPX(IJP))*XPN
      DPYI=0.5D0*(DPY(IJN)+DPY(IJP))*YPN
      FM=DENI*(UI*SX+VI*SY)+CAP*(P(IJN)-P(IJP)-DPXI-DPYI)
C
      RETURN
      END
C
C###############################################################
      SUBROUTINE GRADFIUVP(FI,DFX,DFY)
C###############################################################
      INCLUDE 'com.inc'

      LOGICAL LWRITE,LREAD,LCAL,LAXIS,LTIME

      DIMENSION FI(NXY),DFX(NXY),DFY(NXY),DFXO(NXY),DFYO(NXY)
C
      DO IJ=1,NIJS
        DFXO(IJ)=0.0D0
        DFYO(IJ)=0.0D0
      END DO
C
C.....SET INDICES, INITIALIZE FIELDS
C
      DO IJ=1,NIJS
        DFX(IJ)=0.0D0
        DFY(IJ)=0.0D0
      END DO
C
C.....CONTRIBUTION FROM INNER EAST SIDES
C
      DO IB=1,1
      DO I=2,NIM(IB)-1
         DO J=2,NJM(IB)
         IJ=LB(IB)+LI(IB,I)+J
         IJP=IJ
         IJN=IJ+NJ(IB)
         IJ1=IJP
         IJ2=IJP-1
         FAC=FX(IJ)
         FACP=1.0D0-FAC
         XI=XC(IJN)*FAC+XC(IJP)*FACP
         YI=YC(IJN)*FAC+YC(IJP)*FACP
         XF=0.5D0*(X(IJ1)+X(IJ2))
         YF=0.5D0*(Y(IJ1)+Y(IJ2))
         DFXI=DFXO(IJN)*FAC+DFXO(IJP)*FACP
         DFYI=DFYO(IJN)*FAC+DFYO(IJP)*FACP
C
C.....COORDINATES OF THE CELL-FACE CENTER, VARIABLE VALUE THERE 
C
         FIE=FI(IJN)*FAC+FI(IJP)*FACP+DFXI*(XF-XI)+DFYI*(YF-YI)
          SX=SEX(IJ)
          SY=SEY(IJ)
          DFXE=FIE*SX
          DFYE=FIE*SY
C
          DFX(IJ)=DFX(IJ)+DFXE
          DFY(IJ)=DFY(IJ)+DFYE
          DFX(IJ+NJ(IB))=DFX(IJ+NJ(IB))-DFXE
          DFY(IJ+NJ(IB))=DFY(IJ+NJ(IB))-DFYE
        END DO
      END DO
      ENDDO  
C
C.....CONTRIBUTION FROM INNER NORTH SIDES
C
      DO IB=1,1
      DO I=2,NIM(IB)
         DO J=2,NJM(IB)-1
         IJ=LB(IB)+LI(IB,I)+J
         IJP=IJ
         IJN=IJ+1
         IJ1=IJP
         IJ2=IJP-NJ(IB)
         FAC=FY(IJ)
         FACP=1.0D0-FAC
         XI=XC(IJN)*FAC+XC(IJP)*FACP
         YI=YC(IJN)*FAC+YC(IJP)*FACP
         XF=0.5D0*(X(IJ1)+X(IJ2))
         YF=0.5D0*(Y(IJ1)+Y(IJ2))
         DFXI=DFXO(IJN)*FAC+DFXO(IJP)*FACP
         DFYI=DFYO(IJN)*FAC+DFYO(IJP)*FACP
C
C.....COORDINATES OF THE CELL-FACE CENTER, VARIABLE VALUE THERE 
C
         FIN=FI(IJN)*FAC+FI(IJP)*FACP+DFXI*(XF-XI)+DFYI*(YF-YI)
          SX =SNX(IJ)
          SY =SNY(IJ)
          DFXN=FIN*SX
          DFYN=FIN*SY
C
          DFX(IJ)=DFX(IJ)+DFXN
          DFY(IJ)=DFY(IJ)+DFYN
          DFX(IJ+1)=DFX(IJ+1)-DFXN
          DFY(IJ+1)=DFY(IJ+1)-DFYN
        END DO
      END DO
      ENDDO
C
C.....CONTRIBUTION FROM WALL BOUNDARIES
C
       DO IB=1,1
C
       IF(IB.EQ.1)THEN 
          DO I=2,NIM(IB)
             IJ=LB(IB)+LI(IB,I)+2
             IJB=IJ-1
             SX=SSX(IJ)
             SY=SSY(IJ)
             DFX(IJ)=DFX(IJ)+FI(IJB)*SX
             DFY(IJ)=DFY(IJ)+FI(IJB)*SY
          ENDDO
C
          DO I=2,NIM(IB)
             IJ=LB(IB)+LI(IB,I)+NJM(IB)
             IJB=IJ+1
             SX=SNX(IJ)
             SY=SNY(IJ)
             DFX(IJ)=DFX(IJ)+FI(IJB)*SX
             DFY(IJ)=DFY(IJ)+FI(IJB)*SY
          ENDDO
C
          DO J=2,NJM(IB)
             IJ=LB(IB)+LI(IB,2)+J
             IJB=IJ-NJ(IB)
             SX=SWX(IJ)
             SY=SWY(IJ)
             DFX(IJ)=DFX(IJ)+FI(IJB)*SX
             DFY(IJ)=DFY(IJ)+FI(IJB)*SY
          ENDDO
C
          DO J=2,NJM(IB)
             IJ=LB(IB)+LI(IB,NIM(IB))+J
             IJB=IJ+NJ(IB)
             SX=SEX(IJ)
             SY=SEY(IJ)
             DFX(IJ)=DFX(IJ)+FI(IJB)*SX
             DFY(IJ)=DFY(IJ)+FI(IJB)*SY
          ENDDO
       ENDIF
C
       ENDDO 
C
C.....CALCULATE GRADIENT COMPONENTS AT CV-CENTERS
C
      DO IB=1,1
      DO I=2,NIM(IB)
        DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)
          DFX(IJ)=DFX(IJ)/VOL(IJ)
          DFY(IJ)=DFY(IJ)/VOL(IJ)
        END DO
      END DO
      ENDDO
C
      RETURN
      END
C
C###############################################################
      SUBROUTINE GRADFI(FI,DFX,DFY)
C###############################################################
      INCLUDE 'com.inc'

      LOGICAL LWRITE,LREAD,LCAL,LAXIS,LTIME

      DIMENSION FI(NXY),DFX(NXY),DFY(NXY),DFXO(NXY),DFYO(NXY)
C
      DO IJ=1,NIJS
        DFXO(IJ)=0.0D0
        DFYO(IJ)=0.0D0
      END DO
C.....SET INDICES, INITIALIZE FIELDS
C
      DO IJ=1,NIJS
        DFX(IJ)=0.0D0
        DFY(IJ)=0.0D0
      END DO
C
C.....CONTRIBUTION FROM INNER EAST SIDES
C
      DO IB=1,MIF
      DO I=2,NIM(IB)-1
         DO J=2,NJM(IB)
         IJ=LB(IB)+LI(IB,I)+J
         IJP=IJ
         IJN=IJ+NJ(IB)
         IJ1=IJP
         IJ2=IJP-1
         FAC=FX(IJ)
         FACP=1.0D0-FAC
         XI=XC(IJN)*FAC+XC(IJP)*FACP
         YI=YC(IJN)*FAC+YC(IJP)*FACP
         XF=0.5D0*(X(IJ1)+X(IJ2))
         YF=0.5D0*(Y(IJ1)+Y(IJ2))
         DFXI=DFXO(IJN)*FAC+DFXO(IJP)*FACP
         DFYI=DFYO(IJN)*FAC+DFYO(IJP)*FACP
C
C.....COORDINATES OF THE CELL-FACE CENTER, VARIABLE VALUE THERE 
C
         FIE=FI(IJN)*FAC+FI(IJP)*FACP+DFXI*(XF-XI)+DFYI*(YF-YI)
          SX=SEX(IJ)
          SY=SEY(IJ)
          DFXE=FIE*SX
          DFYE=FIE*SY
C
          DFX(IJ)=DFX(IJ)+DFXE
          DFY(IJ)=DFY(IJ)+DFYE
          DFX(IJ+NJ(IB))=DFX(IJ+NJ(IB))-DFXE
          DFY(IJ+NJ(IB))=DFY(IJ+NJ(IB))-DFYE
        END DO
      END DO
      ENDDO  
C
C.....CONTRIBUTION FROM INNER NORTH SIDES
C
      DO IB=1,MIF
      DO I=2,NIM(IB)
         DO J=2,NJM(IB)-1
         IJ=LB(IB)+LI(IB,I)+J
         IJP=IJ
         IJN=IJ+1
         IJ1=IJP
         IJ2=IJP-NJ(IB)
         FAC=FY(IJ)
         FACP=1.0D0-FAC
         XI=XC(IJN)*FAC+XC(IJP)*FACP
         YI=YC(IJN)*FAC+YC(IJP)*FACP
         XF=0.5D0*(X(IJ1)+X(IJ2))
         YF=0.5D0*(Y(IJ1)+Y(IJ2))
         DFXI=DFXO(IJN)*FAC+DFXO(IJP)*FACP
         DFYI=DFYO(IJN)*FAC+DFYO(IJP)*FACP
C
C.....COORDINATES OF THE CELL-FACE CENTER, VARIABLE VALUE THERE 
C
         FIN=FI(IJN)*FAC+FI(IJP)*FACP+DFXI*(XF-XI)+DFYI*(YF-YI)
          SX =SNX(IJ)
          SY =SNY(IJ)
          DFXN=FIN*SX
          DFYN=FIN*SY
C
          DFX(IJ)=DFX(IJ)+DFXN
          DFY(IJ)=DFY(IJ)+DFYN
          DFX(IJ+1)=DFX(IJ+1)-DFXN
          DFY(IJ+1)=DFY(IJ+1)-DFYN
        END DO
      END DO
      ENDDO
C
C.....CONTRIBUTION FROM INTERFACE
C
      DO I=1,NPIF
        IJP=IJL(I)
        IJN=IJR(I)
        IJ1=IJIF1(I)
        IJ2=IJIF2(I)
        FAC=FB(I)
        FACP=1.0D0-FAC
        XI=XC(IJN)*FAC+XC(IJP)*FACP
        YI=YC(IJN)*FAC+YC(IJP)*FACP
        XF=0.5D0*(X(IJ1)+X(IJ2))
        YF=0.5D0*(Y(IJ1)+Y(IJ2))
        DFXI=DFXO(IJN)*FAC+DFXO(IJP)*FACP
        DFYI=DFYO(IJN)*FAC+DFYO(IJP)*FACP
C
C.....COORDINATES OF THE CELL-FACE CENTER, VARIABLE VALUE THERE 
C
        FII=FI(IJN)*FAC+FI(IJP)*FACP+DFXI*(XF-XI)+DFYI*(YF-YI)
        SX=SBIX(I)
        SY=SBIY(I)
        DFXE=FII*SX
        DFYE=FII*SY
        DFX(IJP)=DFX(IJP)+DFXE
        DFY(IJP)=DFY(IJP)+DFYE
        DFX(IJN)=DFX(IJN)-DFXE
        DFY(IJN)=DFY(IJN)-DFYE
      END DO
C
C.....CONTRIBUTION FROM WALL BOUNDARIES
C
       DO IB=1,MIF
C
       IF(IB.EQ.1)THEN 
          DO I=2,NIM(IB)
             IJ=LB(IB)+LI(IB,I)+2
             IJB=IJ-1
             SX=SSX(IJ)
             SY=SSY(IJ)
             DFX(IJ)=DFX(IJ)+FI(IJB)*SX
             DFY(IJ)=DFY(IJ)+FI(IJB)*SY
          ENDDO
C
          DO J=2,NJM(IB)
             IJ=LB(IB)+LI(IB,2)+J
             IJB=IJ-NJ(IB)
             SX=SWX(IJ)
             SY=SWY(IJ)
             DFX(IJ)=DFX(IJ)+FI(IJB)*SX
             DFY(IJ)=DFY(IJ)+FI(IJB)*SY
          ENDDO
C
       ENDIF
C
       IF(IB.EQ.2)THEN
          DO I=2,NIM(IB)
             IJ=LB(IB)+LI(IB,I)+2
             IJB=IJ-1
             SX=SSX(IJ)
             SY=SSY(IJ)
             DFX(IJ)=DFX(IJ)+FI(IJB)*SX
             DFY(IJ)=DFY(IJ)+FI(IJB)*SY
          ENDDO
C 
          DO J=2,NJM(IB)
             IJ=LB(IB)+LI(IB,NIM(IB))+J
             IJB=IJ+NJ(IB)
             SX=SEX(IJ)
             SY=SEY(IJ)
             DFX(IJ)=DFX(IJ)+FI(IJB)*SX
             DFY(IJ)=DFY(IJ)+FI(IJB)*SY
          ENDDO
       ENDIF
C
       IF(IB.EQ.3)THEN
          DO I=2,NIM(IB)
             IJ=LB(IB)+LI(IB,I)+NJM(IB)
             IJB=IJ+1
             SX=SNX(IJ)
             SY=SNY(IJ)
             DFX(IJ)=DFX(IJ)+FI(IJB)*SX
             DFY(IJ)=DFY(IJ)+FI(IJB)*SY
          ENDDO
          DO J=2,NJM(IB)
             IJ=LB(IB)+LI(IB,2)+J
             IJB=IJ-NJ(IB)
             SX=SWX(IJ)
             SY=SWY(IJ)
             DFX(IJ)=DFX(IJ)+FI(IJB)*SX
             DFY(IJ)=DFY(IJ)+FI(IJB)*SY
          ENDDO
C 
          DO J=2,NJM(IB)
             IJ=LB(IB)+LI(IB,NIM(IB))+J
             IJB=IJ+NJ(IB)
             SX=SEX(IJ)
             SY=SEY(IJ)
             DFX(IJ)=DFX(IJ)+FI(IJB)*SX
             DFY(IJ)=DFY(IJ)+FI(IJB)*SY
          ENDDO
       ENDIF
C
       ENDDO 
C
C.....CALCULATE GRADIENT COMPONENTS AT CV-CENTERS
C
      DO IB=1,MIF
      DO I=2,NIM(IB)
        DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)
          DFX(IJ)=DFX(IJ)/VOL(IJ)
          DFY(IJ)=DFY(IJ)/VOL(IJ)
        END DO
      END DO
      ENDDO
C
      RETURN
      END
C
C
C#############################################################
      SUBROUTINE CALCSC
C#############################################################
      INCLUDE 'com.inc'
      LOGICAL LWRITE,LREAD,LCAL,LAXIS,LTIME
C
C.....CALCULATE GRADIENTS OF FI
C
      CALL GRADFI(T,DTX,DTY)
C
C.....INITIALIZE ARRAYS, SET BLENDING AND UNDER_RELAXATION COEFF.
C
      DO IJ=1,NIJS
        SU(IJ)=0.0D0
        AP(IJ)=0.0D0
      END DO
C
      GFI=GDS(IEN)
      URFFI=1.0D0/URF(IEN)
C
C.....CALCULATE FLUXES THROUGH INNER CV-FACES (EAST & NORTH)
C
      DO IB=1,MIF
      DO I=2,NIM(IB)-1
      DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)
       CALL FLUXSC(IJ,IJ+NJ(IB),IJ,IJ-1,AW(IJ+NJ(IB)),AE(IJ),
     *              FX(IJ),GFI,SEX(IJ),SEY(IJ),F1(IJ))
      END DO
      END DO
      ENDDO
C
      DO IB=1,MIF
      DO I=2,NIM(IB)
      DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)-1
        CALL FLUXSC(IJ,IJ+1,IJ-NJ(IB),IJ,AS(IJ+1),AN(IJ),
     *              FY(IJ),GFI,SNX(IJ),SNY(IJ),F2(IJ))
      END DO
      END DO
      ENDDO
C
C.....CONTRIBUTION FROM INTERFACE
C
      DO I=1,NPIF
        IJP=IJL(I)
        IJN=IJR(I)
        CALL FLUXSC(IJP,IJN,IJIF1(I),IJIF2(I),
     *              AL(I),AR(I),FB(I),GFI,SBIX(I),SBIY(I),0.0D0)
        AP(IJP)=AP(IJP)-AR(I)
        AP(IJN)=AP(IJN)-AL(I)
      END DO
C
C.....WALL BOUNDARY CONDITIONS AND SOURCES FOR TEMPERATURE
C
       CALL TEMP
C
C.....FINAL COEFFICIENT AND SOURCE MATRIX FOR FI-EQUATION
C
      DO IB=1,MIF
      DO I=2,NIM(IB)
      DO J=2,NJM(IB)
         IJ=LB(IB)+LI(IB,I)+J
        AP(IJ)=(AP(IJ)-AE(IJ)-AW(IJ)-AN(IJ)-AS(IJ))*URFFI
        SU(IJ)=SU(IJ)+(1.0D0-URF(IEN))*AP(IJ)*T(IJ)
      END DO
      END DO
      ENDDO
C
C.....SOLVING EQUATION SYSTEM FOR FI-EQUATION
C
      CALL SIPSOL(T,IEN,1)
C
c
      RETURN
      END 
C
C
C################################################################
      SUBROUTINE FLUXSC(IJP,IJN,IJ1,IJ2,CAP,CAN,FAC,G,SX,SY,FM)
C################################################################
      INCLUDE 'com.inc'
C
C.....INTERPOLATED CELL FACE VALUES
C
      FACP=1.0D0-FAC

      XI=XC(IJN)*FAC+XC(IJP)*FACP
      YI=YC(IJN)*FAC+YC(IJP)*FACP
      XF=0.5D0*(X(IJ1)+X(IJ2))
      YF=0.5D0*(Y(IJ1)+Y(IJ2))
C
      DTXI=DTX(IJN)*FAC+DTX(IJP)*FACP
      DTYI=DTY(IJN)*FAC+DTY(IJP)*FACP
      TI=T(IJN)*FAC+T(IJP)*FACP+DTXI*(XF-XI)+DTYI*(YF-YI)
C
C.....SURFACE AND DISTANCE VECTOR COMPONENTS, DIFFUSION COEFF.
C
      TERM1=FACP/TCON(IJP)
      TERM2=FAC/TCON(IJN)
      TCONI=1.0D0/(TERM1+TERM2)
      XPN=XC(IJN)-XC(IJP)
      YPN=YC(IJN)-YC(IJP)
      VSOL=TCONI*SQRT((SX**2+SY**2+1.0D-20)/(XPN**2+YPN**2))
C.....EXPLICIT CONVECTIVE AND DIFFUSIVE FLUXES
C
      FCFIE=FM*TI
      FDFIE=TCONI*(DTXI*SX+DTYI*SY)
C
C.....IMPLICIT CONVECTIVE AND DIFFUSIVE FLUXES
C
      FCFII=MIN(FM,0.0D0)*T(IJN)+MAX(FM,0.0D0)*T(IJP)
      FDFII=VSOL*(DTXI*XPN+DTYI*YPN)
C
C.....COEFFICIENTS, DEFERRED CORRECTION, SOURCE TERMS
C
      CAN=-VSOL+MIN(FM,0.0D0)
      CAP=-VSOL-MAX(FM,0.0D0)
      FFIC=G*(FCFIE-FCFII)

      SU(IJP)=SU(IJP)-FFIC+(FDFIE-FDFII)
      SU(IJN)=SU(IJN)+FFIC+(-FDFIE+FDFII)
C
      RETURN
      END
C
C
C############################################################### 
      SUBROUTINE TEMP 
C###############################################################
      INCLUDE 'com.inc'
      LOGICAL LWRITE,LREAD,LCAL,LAXIS,LTIME
C
      DO IB=1,MIF
CCCCCCCCCCCCCCCC
         IF(IB.EQ.1)THEN
C........WEST (GIVEN BOUNDARY CONDITION)
         DO J=2,NJM(IB)
            IJ=LB(IB)+LI(IB,2)+J     
            IJB=IJ-NJ(IB)
            SX=SWX(IJ)
            SY=SWY(IJ)
            ARR=SQRT(SX**2+SY**2+1.0D-20)
            XNN=SX/(ARR+1.0D-20)
            YNN=SY/(ARR+1.0D-20)
            DNI=(XC(IJB)-XC(IJ))*XNN+(YC(IJB)-YC(IJ))*YNN
            SRDW=ARR/(DNI+1.0D-20)
            COEF=TCON(IJB)*SRDW
            AP(IJ)=AP(IJ)+COEF
            SU(IJ)=SU(IJ)+COEF*T(IJB)
         ENDDO
C
C........SOUTH (ADIABATIC BOUNDARY CONDITION)
         DO I=2,NIM(IB)
            IJ=LB(IB)+LI(IB,I)+2
            IJB=IJ-1
            SX=SSX(IJ)
            SY=SSY(IJ)
            ARR=SQRT(SX**2+SY**2+1.0D-20)
            XNN=SX/(ARR+1.0D-20)
            YNN=SY/(ARR+1.0D-20)
            DNI=(XC(IJB)-XC(IJ))*XNN+(YC(IJB)-YC(IJ))*YNN
            T(IJB)=T(IJ)+DTX(IJ)*(XC(IJB)-XC(IJ)-DNI*XNN)
     *                     +DTY(IJ)*(YC(IJB)-YC(IJ)-DNI*YNN)
         ENDDO
C
         ENDIF
CCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCC
         IF(IB.EQ.2)THEN
C........SOUTH (ADIABATIC BOUNDARY CONDITION)
         DO I=2,NIM(IB)
            IJ=LB(IB)+LI(IB,I)+2
            IJB=IJ-1
            SX=SSX(IJ)
            SY=SSY(IJ)
            ARR=SQRT(SX**2+SY**2+1.0D-20)
            XNN=SX/(ARR+1.0D-20)
            YNN=SY/(ARR+1.0D-20)
            DNI=(XC(IJB)-XC(IJ))*XNN+(YC(IJB)-YC(IJ))*YNN
            T(IJB)=T(IJ)+DTX(IJ)*(XC(IJB)-XC(IJ)-DNI*XNN)
     *                     +DTY(IJ)*(YC(IJB)-YC(IJ)-DNI*YNN)
         ENDDO
C
C........EAST (GIVEN BOUNDARY CONDITION)
         DO J=2,NJM(IB)
            IJ=LB(IB)+LI(IB,NIM(IB))+J     
            IJB=IJ+NJ(IB)
            SX=SEX(IJ)
            SY=SEY(IJ)
            ARR=SQRT(SX**2+SY**2+1.0D-20)
            XNN=SX/(ARR+1.0D-20)
            YNN=SY/(ARR+1.0D-20)
            DNI=(XC(IJB)-XC(IJ))*XNN+(YC(IJB)-YC(IJ))*YNN
            SRDW=ARR/(DNI+1.0D-20)
            COEF=TCON(IJB)*SRDW
            COEFF=1.0D0/(1.0D0/COEF+1.0D0/(ARR*BI))
            AP(IJ)=AP(IJ)+COEFF
            SU(IJ)=SU(IJ)+COEFF*TREFATM
            TERM=BI*(XC(IJB)-XC(IJ))/TCON(IJB)
            T(IJB)=(T(IJ)+TERM*TREFATM)/(1.0D0+TERM)
         ENDDO
C
         ENDIF
CCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCC
         IF(IB.EQ.3)THEN
C........NORTH (GIVEN BOUNDARY CONDITION)
         DO I=2,NIM(IB)
            IJ=LB(IB)+LI(IB,I)+NJM(IB)
            IJB=IJ+1
            SX=SNX(IJ)
            SY=SNY(IJ)
            ARR=SQRT(SX**2+SY**2+1.0D-20)
            XNN=SX/(ARR+1.0D-20)
            YNN=SY/(ARR+1.0D-20)
            DNI=(XC(IJB)-XC(IJ))*XNN+(YC(IJB)-YC(IJ))*YNN
            SRDW=ARR/(DNI+1.0D-20)
            COEF=TCON(IJB)*SRDW
            AP(IJ)=AP(IJ)+COEF
            SU(IJ)=SU(IJ)+COEF*T(IJB)
         ENDDO
C
C........WEST (GIVEN BOUNDARY CONDITION)
         DO J=2,NJM(IB)
            IJ=LB(IB)+LI(IB,2)+J     
            IJB=IJ-NJ(IB)
            SX=SWX(IJ)
            SY=SWY(IJ)
            ARR=SQRT(SX**2+SY**2+1.0D-20)
            XNN=SX/(ARR+1.0D-20)
            YNN=SY/(ARR+1.0D-20)
            DNI=(XC(IJB)-XC(IJ))*XNN+(YC(IJB)-YC(IJ))*YNN
            SRDW=ARR/(DNI+1.0D-20)
            COEF=TCON(IJB)*SRDW
            AP(IJ)=AP(IJ)+COEF
            SU(IJ)=SU(IJ)+COEF*T(IJB)
         ENDDO
C
C........EAST (GIVEN BOUNDARY CONDITION)
         DO J=2,NJM(IB)
            IJ=LB(IB)+LI(IB,NIM(IB))+J     
            IJB=IJ+NJ(IB)
            SX=SEX(IJ)
            SY=SEY(IJ)
            ARR=SQRT(SX**2+SY**2+1.0D-20)
            XNN=SX/(ARR+1.0D-20)
            YNN=SY/(ARR+1.0D-20)
            DNI=(XC(IJB)-XC(IJ))*XNN+(YC(IJB)-YC(IJ))*YNN
            SRDW=ARR/(DNI+1.0D-20)
            COEF=TCON(IJB)*SRDW
            COEFF=1.0D0/(1.0D0/COEF+1.0D0/(ARR*BI))
            AP(IJ)=AP(IJ)+COEFF
            SU(IJ)=SU(IJ)+COEFF*TREFATM
            TERM=BI*(XC(IJB)-XC(IJ))/TCON(IJB)
            T(IJB)=(T(IJ)+TERM*TREFATM)/(1.0D0+TERM)
         ENDDO
C
         ENDIF
CCCCCCCCCCCC
      ENDDO
C
      RETURN
      END
C
C#############################################################
      SUBROUTINE SIPSOL(FI,IFI,LISS)
C#############################################################
      INCLUDE 'com.inc'
      DIMENSION FI(NXY),UE(NXY),UN(NXY),RES(NXY)
      REAL BLW(NXY),BLS(NXY),BLPR(NXY)
      LOGICAL LWRITE,LREAD,LCAL,LAXIS,LTIME
C-------------------------------------------------------------
      DATA UE,UN,RES /NXY*0.0D0,NXY*0.0D0,NXY*0.0D0/
C
      
C.....COEFFICIENTS OF UPPER AND LOWER TRIANGULAR MATRICES
C
      IF(LISS.EQ.1)THEN
      DO IB=1,MIF
      DO I=2,NIM(IB)
        DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)
          BLW(IJ)=AW(IJ)/(1.0D0+ALFA*UN(IJ-NJ(IB)))
          BLS(IJ)=AS(IJ)/(1.0D0+ALFA*UE(IJ-1))
          P1=ALFA*BLW(IJ)*UN(IJ-NJ(IB))
          P2=ALFA*BLS(IJ)*UE(IJ-1)
          BLPR(IJ)=1.0D0/(AP(IJ)+P1+P2-BLW(IJ)*UE(IJ-NJ(IB))-
     *BLS(IJ)*UN(IJ-1)+1.0E-20)
          UN(IJ)=(AN(IJ)-P1)*BLPR(IJ)
          UE(IJ)=(AE(IJ)-P2)*BLPR(IJ)
C
        END DO
      END DO
      ENDDO
      ENDIF  
C
C==============================================================
C.....INNER ITERATIONS LOOP
C==============================================================
C
      DO L=1,NSW(IFI)
        RESL=1.0D-20
C      
C.....CALCULATE RESIDUAL AND OVERWRITE IT BY INTERMEDIATE VECTOR
C
        DO IB=1,MIF
        DO I=2,NIM(IB)
        DO J=2,NJM(IB)
           IJ=LB(IB)+LI(IB,I)+J
          RES(IJ)=SU(IJ)-AN(IJ)*FI(IJ+1)-AS(IJ)*FI(IJ-1)-
     *          AE(IJ)*FI(IJ+NJ(IB))-AW(IJ)*FI(IJ-NJ(IB))-AP(IJ)*FI(IJ)
        END DO
        END DO
        ENDDO
C
C
        DO I=1,NPIF
          RES(IJL(I))=RES(IJL(I))-AR(I)*FI(IJR(I))
          RES(IJR(I))=RES(IJR(I))-AL(I)*FI(IJL(I))
        END DO
C
C.....FORWARD SUBSTITUTION
C
        DO IB=1,MIF
        DO I=2,NIM(IB)
        DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)
          RESL=RESL+ABS(RES(IJ))
          RES(IJ)=(RES(IJ)-BLS(IJ)*RES(IJ-1)-BLW(IJ)*
     *            RES(IJ-NJ(IB)))*BLPR(IJ)
        END DO
        END DO
        ENDDO 
C
C.....STORE INITIAL RESIDUAL SUM FOR CHECKING CONV. OF OUTER ITER.
C
        IF(L.EQ.1) RESOR(IFI)=RESL
        RSM=RESL/(RESOR(IFI)+1.0E-20)
C
C.....BACK SUBSTITUTION AND CORRECTION
C
        DO IB=1,MIF
        DO I=NIM(IB),2,-1
        DO IJ=LB(IB)+LI(IB,I)+NJM(IB),LB(IB)+LI(IB,I)+2,-1
          RES(IJ)=RES(IJ)-UN(IJ)*RES(IJ+1)-UE(IJ)*RES(IJ+NJ(IB))
          FI(IJ)=FI(IJ)+RES(IJ)
        END DO
        END DO
        ENDDO 
C
C.....CHECK CONVERGENCE OF INNER ITERATIONS
C
        IF(RSM.LT.SOR(IFI)) RETURN
C
      END DO
C
      RETURN
      END
C
C#############################################################
      SUBROUTINE SIPSOLUVP(FI,IFI,LISS)
C#############################################################
      INCLUDE 'com.inc'
      DIMENSION FI(NXY),UE(NXY),UN(NXY),RES(NXY)
      REAL BLW(NXY),BLS(NXY),BLPR(NXY)
      LOGICAL LWRITE,LREAD,LCAL,LAXIS,LTIME
C-------------------------------------------------------------
      DATA UE,UN,RES /NXY*0.0D0,NXY*0.0D0,NXY*0.0D0/
C
      
C.....COEFFICIENTS OF UPPER AND LOWER TRIANGULAR MATRICES
C
      IF(LISS.EQ.1)THEN
      DO IB=1,1
      DO I=2,NIM(IB)
        DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)
          BLW(IJ)=AW(IJ)/(1.0D0+ALFA*UN(IJ-NJ(IB)))
          BLS(IJ)=AS(IJ)/(1.0D0+ALFA*UE(IJ-1))
          P1=ALFA*BLW(IJ)*UN(IJ-NJ(IB))
          P2=ALFA*BLS(IJ)*UE(IJ-1)
          BLPR(IJ)=1.0D0/(AP(IJ)+P1+P2-BLW(IJ)*UE(IJ-NJ(IB))-
     *BLS(IJ)*UN(IJ-1)+1.0E-20)
          UN(IJ)=(AN(IJ)-P1)*BLPR(IJ)
          UE(IJ)=(AE(IJ)-P2)*BLPR(IJ)
C
        END DO
      END DO
      ENDDO
      ENDIF  
C
C==============================================================
C.....INNER ITERATIONS LOOP
C==============================================================
C
      DO L=1,NSW(IFI)
        RESL=1.0D-20
C      
C.....CALCULATE RESIDUAL AND OVERWRITE IT BY INTERMEDIATE VECTOR
C
        DO IB=1,1
        DO I=2,NIM(IB)
        DO J=2,NJM(IB)
           IJ=LB(IB)+LI(IB,I)+J
          RES(IJ)=SU(IJ)-AN(IJ)*FI(IJ+1)-AS(IJ)*FI(IJ-1)-
     *          AE(IJ)*FI(IJ+NJ(IB))-AW(IJ)*FI(IJ-NJ(IB))-AP(IJ)*FI(IJ)
        END DO
        END DO
        ENDDO
C
C.....FORWARD SUBSTITUTION
C
        DO IB=1,1
        DO I=2,NIM(IB)
        DO IJ=LB(IB)+LI(IB,I)+2,LB(IB)+LI(IB,I)+NJM(IB)
          RESL=RESL+ABS(RES(IJ))
          RES(IJ)=(RES(IJ)-BLS(IJ)*RES(IJ-1)-BLW(IJ)*
     *            RES(IJ-NJ(IB)))*BLPR(IJ)
        END DO
        END DO
        ENDDO 
C
C.....STORE INITIAL RESIDUAL SUM FOR CHECKING CONV. OF OUTER ITER.
C
        IF(L.EQ.1) RESOR(IFI)=RESL
        RSM=RESL/(RESOR(IFI)+1.0E-20)
C
C.....BACK SUBSTITUTION AND CORRECTION
C
        DO IB=1,1
        DO I=NIM(IB),2,-1
        DO IJ=LB(IB)+LI(IB,I)+NJM(IB),LB(IB)+LI(IB,I)+2,-1
          RES(IJ)=RES(IJ)-UN(IJ)*RES(IJ+1)-UE(IJ)*RES(IJ+NJ(IB))
          FI(IJ)=FI(IJ)+RES(IJ)
        END DO
        END DO
        ENDDO 
C
C.....CHECK CONVERGENCE OF INNER ITERATIONS
C
        IF(RSM.LT.SOR(IFI)) RETURN
C
      END DO
C
      RETURN
      END
C
C#############################################################
      SUBROUTINE PBOUND(FI)
C#############################################################
      INCLUDE'com.inc'

      DIMENSION FI(NXY)
C--------------------------------------------------------------
C
C.....EXTRAPOLATE TO SOUTH AN NORTH BOUNDARIES
C
      DO IB=1,1
C
      DO I=2,NIM(IB)
        IJ=LB(IB)+LI(IB,I)+1
        FI(IJ)=FI(IJ+1)+(FI(IJ+1)-FI(IJ+2))*FY(IJ+1)
        IJ=LB(IB)+LI(IB,I)+NJ(IB)
        FI(IJ)=FI(IJ-1)+(FI(IJ-1)-FI(IJ-2))*(1.0D0-FY(IJ-2))
      END DO
C
C.....EXTRAPOLATE TO WEST AND EAST BOUNDARIES
C
      DO J=2,NJM(IB)
        IJ=LB(IB)+LI(IB,1)+J
        FI(IJ)=FI(IJ+NJ(IB))+(FI(IJ+NJ(IB))-FI(IJ+2*NJ(IB)))*
     *FX(IJ+NJ(IB)) 
      ENDDO
C
      DO J=2,NJM(IB)
        IJ=LB(IB)+LI(IB,NI(IB))+J
        FI(IJ)=FI(IJ-NJ(IB))+(FI(IJ-NJ(IB))-FI(IJ-2*NJ(IB)))*
     *(1.0D0-FX(IJ-2*NJ(IB)))
      END DO
C
      ENDDO  
C
C
      RETURN
      END
C
C############################################################
      SUBROUTINE MODINP
C############################################################
      INCLUDE'com.inc'
      LOGICAL LWRITE,LREAD,LCAL,LAXIS,LTIME
C
C-----------------------------------------------------------
C.....READ INPUT DATA FROM UNIT 5
C-----------------------------------------------------------
C
      READ(5,*) LREAD,LWRITE,LAXIS,LTIME
      READ(5,*) MAXIT,IMON,JMON,IPR,JPR,SORMAX,SLARGE,ALFA
      READ(5,*) DENSIT,VISC,BETA,GRAVX,GRAVY,TREF,TREFATM
      READ(5,*) UIN,VIN,PIN,TIN        
      READ(5,*) ITST,DT,GAMT
      READ(5,*) (LCAL(I),I=1,NPHI)
      READ(5,*) (URF(I),I=1,NPHI)
      READ(5,*) (SOR(I),I=1,NPHI)
      READ(5,*) (NSW(I),I=1,NPHI)
      READ(5,*) (GDS(I),I=1,NPHI)
      READ(5,*) TKON,HKAP,BI
C
      IU=1
      IV=2
      IP=3
      IEN=4
      DTR=1.0D0/DT
      PR=0.7D0      
C
C-----------------------------------------------------------
C.....READ GRID DATA (GENERATED USING GRID GENERATOR FOR MG.)
C-----------------------------------------------------------
C
      READ(1)MIF   
      DO IB=1,MIF
         READ(1)NI(IB)
         READ(1)NJ(IB)
         READ(1)NIM(IB)
         READ(1)NJM(IB)
         READ(1)NIJ(IB)
      ENDDO
      READ(1)NIJS
      READ(1)(X(I),I=1,NIJS),(XC(I),I=1,NIJS)
     *       ,(Y(I),I=1,NIJS),(YC(I),I=1,NIJS)
     *       ,(FX(I),I=1,NIJS),(FY(I),I=1,NIJS)
     *       ,(VOL(I),I=1,NIJS)
     *       ,(SEX(I),I=1,NIJS),(SEY(I),I=1,NIJS)
     *       ,(SWX(I),I=1,NIJS),(SWY(I),I=1,NIJS)
     *       ,(SNX(I),I=1,NIJS),(SNY(I),I=1,NIJS)
     *       ,(SSX(I),I=1,NIJS),(SSY(I),I=1,NIJS)
C
C.....READING DATA FOR INTERFACE
C
      READ(1)NPIF,NP1,NP2,NP3
C
      DO I=1,NPIF
         READ(1)SBIX(I),SBIY(I),FB(I)
      ENDDO
C
      DO I=1,NPIF
         READ(1)IJIF1(I),IJIF2(I),IJL(I),IJR(I)
      ENDDO
C
      DO IB=1,MIF
      DO I=1,NI(IB)
         LI(IB,I)=(I-1)*NJ(IB) 
      ENDDO
      ENDDO
C
      LB(1)=0
      DO IB=2,MIF
         LB(IB)=LB(IB-1)+NI(IB-1)*NJ(IB-1) 
      ENDDO
C
C
      DO IB=1,1
      DO I=1,NI(IB)
         DO J=1,NJ(IB)
            IJ=LB(IB)+LI(IB,I)+J
            T(IJ)=TIN
            TO(IJ)=TIN
            VIS(IJ)=VISC
            F1(IJ)=0.0D0
            F2(IJ)=0.0D0
            U(IJ)=UIN
            V(IJ)=VIN
            P(IJ)=PIN
            DEN(IJ)=DENSIT
            TCON(IJ)=VIS(IJ)/PR
            HCAP(IJ)=HKAP  
         ENDDO
      ENDDO
      ENDDO
C
      DO IB=2,3
      DO I=1,NI(IB)
         DO J=1,NJ(IB)
            IJ=LB(IB)+LI(IB,I)+J
            T(IJ)=TIN
            TO(IJ)=TIN
            VIS(IJ)=VISC
            F1(IJ)=0.0D0
            F2(IJ)=0.0D0
            U(IJ)=UIN
            V(IJ)=VIN
            P(IJ)=PIN
            DEN(IJ)=DENSIT
            TCON(IJ)=VIS(IJ)*100.0D0/PR
            HCAP(IJ)=HKAP
         ENDDO
      ENDDO
      ENDDO
C
      DO IB=2,2
      DO I=1,NI(IB)
         DO J=1,NJ(IB)
            IJ=LB(IB)+LI(IB,I)+J
            TCON(IJ)=VIS(IJ)*100.0D0/PR
         ENDDO
      ENDDO
      ENDDO
C
C
C---------------------------------------------------
C.....BOUNDARY CONDITIONS 
C---------------------------------------------------
C
       DO IB=1,MIF
          IF(IB.EQ.1)THEN
            DO J=1,NJ(IB)
               IJ=LB(IB)+LI(IB,1)+J
               T(IJ)=TREF
            ENDDO
          ENDIF
C
          IF(IB.EQ.3)THEN
            DO J=1,NJ(IB)
               IJ=LB(IB)+LI(IB,1)+J
               T(IJ)=TREF
            ENDDO
C
            DO I=2,NI(IB)
               IJ=LB(IB)+LI(IB,I)+NJ(IB)
               T(IJ)=1.0D0
            ENDDO
          ENDIF
       ENDDO
C
      RETURN
      END
C
