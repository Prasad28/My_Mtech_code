      PROGRAM NONGRID
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NX=502,NY=502,NXY=NX*NY,MB=3)
      DIMENSION X(NXY),Y(NXY),XC(NXY),YC(NXY),FX(NXY),FY(NXY),FB(NXY)
      DIMENSION SEX(NXY),SEY(NXY),SWX(NXY),SWY(NXY),SNX(NXY),SNY(NXY)
     *         ,SSX(NXY),SSY(NXY),SBIX(NXY),SBIY(NXY),VOL(NXY) 
      DIMENSION IJ1(NXY),IJ2(NXY),IJL(NXY),IJR(NXY) 
      DIMENSION LB(MB),LI(MB,NX),NI(MB),NJ(MB),NIM(MB),NJM(MB),NIJ(MB)
C
      CHARACTER FILOUT*20  
C
      OPEN(UNIT=6,FILE='2dgrid',FORM='UNFORMATTED')  
      READ(6)MIF,NIJS
      DO IB=1,MIF
         READ(6)NI(IB),NJ(IB),NIM(IB),NJM(IB),NIJ(IB)
      ENDDO
      READ(6)(X(IJ),IJ=1,NIJS),(Y(IJ),IJ=1,NIJS)
      CLOSE(6)
C
      NB=MIF
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
C.....CALCULATION FOR SURFACE AREA COMPONENTS
         DO IB=1,NB
         DO I=2,NIM(IB)
            DO J=2,NJM(IB)
               IJ=LB(IB)+LI(IB,I)+J
               SEX(IJ)=(Y(IJ)-Y(IJ-1))
               SEY(IJ)=(X(IJ-1)-X(IJ))
               SWX(IJ)=(Y(IJ-NJ(IB)-1)-Y(IJ-NJ(IB)))
               SWY(IJ)=(X(IJ-NJ(IB))-X(IJ-NJ(IB)-1))
               SNX(IJ)=(Y(IJ-NJ(IB))-Y(IJ))
               SNY(IJ)=(X(IJ)-X(IJ-NJ(IB)))
               SSX(IJ)=(Y(IJ-1)-Y(IJ-NJ(IB)-1))
               SSY(IJ)=(X(IJ-NJ(IB)-1)-X(IJ-1))
            ENDDO
         ENDDO
         ENDDO
C
C.....CALCULATION OF CELL VOLUMES FOR PLANE GEOMETRY
C
          DO IB=1,NB
          DO I=2,NIM(IB)
          DO J=2,NJM(IB)
            IJ=LB(IB)+LI(IB,I)+J
            DXNESW=X(IJ)-X(IJ-NJ(IB)-1)
            DYNESW=Y(IJ)-Y(IJ-NJ(IB)-1)
            DXNWSE=X(IJ-NJ(IB))-X(IJ-1)
            DYNWSE=Y(IJ-NJ(IB))-Y(IJ-1)
            VOL(IJ)=0.5*ABS(DXNESW*DYNWSE-DXNWSE*DYNESW)
          END DO
          END DO
          ENDDO
C
          VOLUME=0.0D0
          DO IB=1,NB
          DO I=2,NIM(IB)
          DO J=2,NJM(IB)
            IJ=LB(IB)+LI(IB,I)+J
            VOLUME=VOLUME+VOL(IJ)
          END DO
          END DO
          ENDDO
C
      PRINT*,'VOLUME =',VOLUME
C.....CALCULATION OF NODE COORDINATES: CORNER (DUMMY) NODES
C
       DO IB=1,NB
        IJ=LB(IB)+LI(IB,1)+1
        XC(IJ)=X(IJ)
        YC(IJ)=Y(IJ)
C
        IJ=LB(IB)+LI(IB,NIM(IB))+1
        XC(IJ+NJ(IB))=X(IJ)
        YC(IJ+NJ(IB))=Y(IJ)
C
        IJ=LB(IB)+LI(IB,1)+NJM(IB)
        XC(IJ+1)=X(IJ)
        YC(IJ+1)=Y(IJ)
C
        IJ=LB(IB)+LI(IB,NIM(IB))+NJM(IB)
        XC(IJ+NJ(IB)+1)=X(IJ)
        YC(IJ+NJ(IB)+1)=Y(IJ)
C
C.....CALCULATION OF NODE COORDINATES: BOUNDARY NODES
C
        DO J=2,NJM(IB)
          IJ=LB(IB)+LI(IB,NIM(IB))+J
          XC(IJ+NJ(IB))=0.5D0*(X(IJ)+X(IJ-1))
          YC(IJ+NJ(IB))=0.5D0*(Y(IJ)+Y(IJ-1))
          IJ=LB(IB)+LI(IB,1)+J
          XC(IJ)=0.5D0*(X(IJ)+X(IJ-1))
          YC(IJ)=0.5D0*(Y(IJ)+Y(IJ-1))
        END DO
C
        DO I=2,NIM(IB)
          IJ=LB(IB)+LI(IB,I)+1
          XC(IJ)=0.5D0*(X(IJ)+X(IJ-NJ(IB)))
          YC(IJ)=0.5D0*(Y(IJ)+Y(IJ-NJ(IB)))
          IJ=LB(IB)+LI(IB,I)+NJM(IB)
          XC(IJ+1)=0.5D0*(X(IJ)+X(IJ-NJ(IB)))
          YC(IJ+1)=0.5D0*(Y(IJ)+Y(IJ-NJ(IB)))
        END DO
C
C.....CALCULATION OF NODE COORDINATES: CELL CENTERS
C
        DO I=2,NIM(IB)
        DO J=2,NJM(IB)
          IJ=LB(IB)+LI(IB,I)+J
          XC(IJ)=0.25D0*(X(IJ)+X(IJ-1)+X(IJ-NJ(IB))+X(IJ-NJ(IB)-1))
          YC(IJ)=0.25D0*(Y(IJ)+Y(IJ-1)+Y(IJ-NJ(IB))+Y(IJ-NJ(IB)-1))
        END DO
        END DO
C
      ENDDO
C==========================================================     
C..... CALCULATION OF INTERPOLATION FACTORS
C==========================================================     
C
      DO IB=1,NB
        DO I=2,NIM(IB)
        DO J=2,NJM(IB)
          IJ=LB(IB)+LI(IB,I)+J
C
C.....INTERPOLATION IN I-DIRECTION: FX = Pe/PE
C
          XE=0.5D0*(X(IJ)+X(IJ-1))
          YE=0.5D0*(Y(IJ)+Y(IJ-1))
          DLPE=SQRT((XE-XC(IJ))**2+(YE-YC(IJ))**2)
          DLEE=SQRT((XC(IJ+NJ(IB))-XE)**2+(YC(IJ+NJ(IB))-YE)**2)
          FX(IJ)=DLPE/(DLPE+DLEE+1.E-20)
C
C.....INTERPOLATION IN J-DIRECTION: FY = Pn/PN
C
          XN=0.5D0*(X(IJ)+X(IJ-NJ(IB)))
          YN=0.5D0*(Y(IJ)+Y(IJ-NJ(IB)))
          DLPN=SQRT((XN-XC(IJ))**2+(YN-YC(IJ))**2)
          DLNN=SQRT((XC(IJ+1)-XN)**2+(YC(IJ+1)-YN)**2)
          FY(IJ)=DLPN/(DLPN+DLNN+1.E-20)
        END DO
        END DO
      ENDDO
C
C.....PREPARATION OF INTERFACE DATA
C
      NP1=0
      NP2=0
      NP3=0
C
C.....BLOCK 1 AND BLOCK 2
      DO IB=1,1
         JK=2
  30     DO J=JK,NJM(IB)
            IJP=LB(IB)+LI(IB,NIM(IB))+J
            IJP1=LB(IB)+LI(IB,NIM(IB))+J
            IJP2=LB(IB)+LI(IB,NIM(IB))+J-1
            DO KK=2,NJM(IB+1)
               IJN=LB(IB+1)+LI(IB+1,2)+KK
               IJN1=LB(IB+1)+LI(IB+1,1)+KK
               IJN2=LB(IB+1)+LI(IB+1,1)+KK-1
               IF(ABS(Y(IJP1)-Y(IJN1)).LT.1.0D-10)THEN
                 NP1=NP1+1
                 IJ1(NP1)=IJP1
                 IJ2(NP1)=IJP2
                 IJL(NP1)=IJP
                 IJR(NP1)=IJN
                 JK=JK+1
                 XB=0.5D0*(X(IJ1(NP1))+X(IJ2(NP1)))
                 YB=0.5D0*(Y(IJ1(NP1))+Y(IJ2(NP1)))
                 DLPB=SQRT((XB-XC(IJL(NP1)))**2+(YB-YC(IJL(NP1)))**2)
                 DLNB=SQRT((XC(IJR(NP1))-XB)**2+(YB-YC(IJR(NP1)))**2)
                 FB(NP1)=DLPB/(DLPB+DLNB+1.0D-20)
                 SBIX(NP1)=SEX(IJP)
                 SBIY(NP1)=SEY(IJP)
                 GOTO 30
               ENDIF
            ENDDO
         ENDDO 
      ENDDO
C
      PRINT*,'NP1 =',NP1 
C
      NP2=NP2+NP1
C.....BLOCK 1 AND BLOCK 3
      DO IB=1,1
         IK=2
  40     DO I=IK,NIM(IB)
            IJP=LB(IB)+LI(IB,I)+NJM(IB)
            IJP1=LB(IB)+LI(IB,I)+NJM(IB)
            IJP2=LB(IB)+LI(IB,I-1)+NJM(IB)
            DO KK=2,NIM(IB+2)
               IJN=LB(IB+2)+LI(IB+2,KK)+2
               IJN1=LB(IB+2)+LI(IB+2,KK)+1
               IJN2=LB(IB+2)+LI(IB+2,KK-1)+1
               IF(ABS(X(IJP1)-X(IJN1)).LT.1.0D-10)THEN
                 NP2=NP2+1
                 IJ1(NP2)=IJP1
                 IJ2(NP2)=IJP2
                 IJL(NP2)=IJP
                 IJR(NP2)=IJN
                 IK=IK+1
                 XB=0.5D0*(X(IJ1(NP2))+X(IJ2(NP2)))
                 YB=0.5D0*(Y(IJ1(NP2))+Y(IJ2(NP2)))
                 DLPB=SQRT((XB-XC(IJL(NP2)))**2+(YB-YC(IJL(NP2)))**2)
                 DLNB=SQRT((XC(IJR(NP2))-XB)**2+(YB-YC(IJR(NP2)))**2)
                 FB(NP2)=DLPB/(DLPB+DLNB+1.0D-20)
                 SBIX(NP2)=SNX(IJP)
                 SBIY(NP2)=SNY(IJP)
                 GOTO 40
               ENDIF
            ENDDO
         ENDDO 
      ENDDO
C
      PRINT*,'NP2 =',NP2-NP1 
C
      NP3=NP3+NP2
C.....BLOCK 2 AND BLOCK 3
      DO IB=2,2
         IK=2
  50     DO I=IK,NIM(IB)
            IJP=LB(IB)+LI(IB,I)+NJM(IB)
            IJP1=LB(IB)+LI(IB,I)+NJM(IB)
            IJP2=LB(IB)+LI(IB,I-1)+NJM(IB)
            DO KK=2,NIM(IB+1)
               IJN=LB(IB+1)+LI(IB+1,KK)+2
               IJN1=LB(IB+1)+LI(IB+1,KK)+1
               IJN2=LB(IB+1)+LI(IB+1,KK-1)+1
               IF(ABS(X(IJP1)-X(IJN1)).LT.1.0D-10)THEN
                 PRINT*,'IK KK IJP IJN =',IK,KK,IJP,IJN
                 NP3=NP3+1
                 IJ1(NP3)=IJP1
                 IJ2(NP3)=IJP2
                 IJL(NP3)=IJP
                 IJR(NP3)=IJN
                 IK=IK+1
                 XB=0.5D0*(X(IJ1(NP3))+X(IJ2(NP3)))
                 YB=0.5D0*(Y(IJ1(NP3))+Y(IJ2(NP3)))
                 DLPB=SQRT((XB-XC(IJL(NP3)))**2+(YB-YC(IJL(NP3)))**2)
                 DLNB=SQRT((XC(IJR(NP3))-XB)**2+(YB-YC(IJR(NP3)))**2)
                 FB(NP3)=DLPB/(DLPB+DLNB+1.0D-20)
                 SBIX(NP3)=SNX(IJP)
                 SBIY(NP3)=SNY(IJP)
                 GOTO 50
               ENDIF
            ENDDO
         ENDDO 
      ENDDO
C
      PRINT*,'NP3 =',NP3-NP2 
C
      NP=NP3
      NPIF=NP
      PRINT*,'NPIF =',NPIF
C
C
      OPEN(UNIT=4,FILE='filgr',FORM='UNFORMATTED')  
C
      WRITE(4)MIF   
      DO IB=1,MIF
         WRITE(4)NI(IB)
         WRITE(4)NJ(IB)
         WRITE(4)NIM(IB)
         WRITE(4)NJM(IB)
         WRITE(4)NIJ(IB)
      ENDDO
      WRITE(4)NIJS
      WRITE(4)(X(I),I=1,NIJS),(XC(I),I=1,NIJS)
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
      WRITE(4)NPIF,NP1,NP2,NP3
C
      DO I=1,NPIF
         WRITE(4)SBIX(I),SBIY(I),FB(I)
      ENDDO
C
      DO I=1,NPIF
         WRITE(4)IJ1(I),IJ2(I),IJL(I),IJR(I)
      ENDDO
C
      CLOSE(4)
C
      OPEN(UNIT=2,FILE='blockc.dat')
      WRITE(2,*) 'TITLE = "FALCOL COORDINATE DATA"'
      WRITE(2,*)'VARIABLES = "X", "Y","XC","YC"'
      DO IB=1,NB
      WRITE(FILOUT,'(a5,I1)') 'BLOCK',IB
      WRITE(2,*)'ZONE T=', FILOUT,' ','I=',NJ(IB),' ','J=',NI(IB),' '
     *,'DATAPACKING=POINT'
      DO I=1,NI(IB)
         DO J=1,NJ(IB)
            IJ=LB(IB)+LI(IB,I)+J
           WRITE(2,*)X(IJ),Y(IJ),XC(IJ),YC(IJ)
         ENDDO
      ENDDO 
      ENDDO
      CLOSE(2)
C
      STOP
      END


