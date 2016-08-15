      INCLUDE 'com.inc'
      DIMENSION PSI(NXY),TEMPGRADS(NXY),TEMPGRADT(NXY),TGRADS(NXY)
     *,TGRADT(NXY)
      CHARACTER*20 FILOUT,FILGR,FILRES

      PRINT *,' ENTER RATIO OF THERMAL CONDUCTIVITY:  '
      READ *,RATIO
      FILOUT='res'
      FILRES='tec.dat'
      FILGR='filgr' 
      OPEN (UNIT=2,FILE=FILOUT,FORM='UNFORMATTED')
      OPEN (UNIT=7,FILE=FILRES)
      OPEN (UNIT=9,FILE=FILGR,FORM='UNFORMATTED')
C
      READ(9)MIF   
      DO IB=1,MIF
         READ(9)NI(IB)
         READ(9)NJ(IB)
         READ(9)NIM(IB)
         READ(9)NJM(IB)
         READ(9)NIJ(IB)
      ENDDO
      READ(9)NIJS
      READ(9)(X(I),I=1,NIJS),(XC(I),I=1,NIJS)
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
      READ(9)NPIF,NP1,NP2,NP3
C
      DO I=1,NPIF
         READ(9)SBIX(I),SBIY(I),FB(I)
      ENDDO
C
      DO I=1,NPIF
         READ(9)IJIF1(I),IJIF2(I),IJL(I),IJR(I)
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
      CLOSE(9)
C
      READ(2)(F1(IJ),IJ=1,NIJS),(F2(IJ),IJ=1,NIJS),(U(IJ),IJ=1,NIJS),
     *        (V(IJ),IJ=1,NIJS),(P(IJ),IJ=1,NIJS),(T(IJ),IJ=1,NIJS),
     *       (UO(IJ),IJ=1,NIJS),(VO(IJ),IJ=1,NIJS),(TO(IJ),IJ=1,NIJS)
     *       ,(DTX(IJ),IJ=1,NIJS),(DTY(IJ),IJ=1,NIJS)
      CLOSE(2)
C
      IB=1
      DO J=1,NJ(IB)
         IJ1=LB(IB)+LI(IB,NIM(IB))+J
         IJB1=IJ1+NJ(IB)
         IJ2=LB(IB+1)+LI(IB+1,2)+J
         IJB2=IJ2-NJ(IB+1) 
         FACTOR=RATIO*(XC(IJB1)-XC(IJ1))/(XC(IJ2)-XC(IJB2))
         ANUM=T(IJ1)+FACTOR*T(IJ2)
         DENOM=1.0D0+FACTOR
         T(IJB1)=ANUM/(DENOM+1.0D-20)
         T(IJB2)=T(IJB1)
      ENDDO   
C
      IB=1
      DO I=1,NI(IB)
         IJ1=LB(IB)+LI(IB,I)+NJM(IB)
         IJB1=IJ1+1
         IJ2=LB(IB+2)+LI(IB+2,I)+2
         IJB2=IJ2-1 
         FACTOR=RATIO*(YC(IJB1)-YC(IJ1))/(YC(IJ2)-YC(IJB2))
         ANUM=T(IJ1)+FACTOR*T(IJ2)
         DENOM=1.0D0+FACTOR
         T(IJB1)=ANUM/(DENOM+1.0D-20)
         T(IJB2)=T(IJB1)
      ENDDO   
C
      T(LB(2)+LI(2,1)+NJ(2))=T(LB(1)+LI(1,NI(1))+NJ(1))
C
      IB=2
      DO I=2,NI(IB)
         IJ1=LB(IB)+LI(IB,I)+NJM(IB)
         IJB1=IJ1+1
         IJ2=LB(IB+1)+LI(IB+1,NIM(IB-1)-1+I)+2
         IJB2=IJ2-1 
         FACTOR=(YC(IJB1)-YC(IJ1))/(YC(IJ2)-YC(IJB2))
         ANUM=T(IJ1)+FACTOR*T(IJ2)
         DENOM=1.0D0+FACTOR
         T(IJB1)=ANUM/(DENOM+1.0D-20)
         T(IJB2)=T(IJB1)
      ENDDO 
C
       DO IB=1,MIF 
        DO I=1,NIM(IB)
          DO J=1,NJM(IB)
            IJ=LB(IB)+LI(IB,I)+J
            PSI(IJ)=F1(IJ)
          END DO
        END DO
       ENDDO
C
       PSI(1)=0.0D0
       DO IB=1,1 
        DO I=1,NIM(IB)
          II=LB(IB)+LI(IB,I)
          IF((I.NE.1))PSI(II+1)=PSI(II-NJ(IB)+1)-F2(II+1)
          DO J=2,NJM(IB)
            IJ=II+J
            PSI(IJ)=PSI(IJ-1)+PSI(IJ)
          END DO
        END DO
       ENDDO
C
      WRITE(7,*)'TITLE = "Example: Multi-Zone 2-D Plot"'
      WRITE(7,*)'VARIABLES = "X", "Y","XC","YC","U","V","P","T","S"'
      DO IB=1,MIF
      WRITE(FILOUT,'(a5,I1)') 'BLOCK',IB
      WRITE(7,*)'ZONE T=', FILOUT,' ','I=',NJ(IB),' ','J=',NI(IB),' ',
     *'DATAPACKING=POINT'
      DO I=1,NI(IB)
         DO J=1,NJ(IB)
            IJ=LB(IB)+LI(IB,I)+J
           WRITE(7,*)X(IJ),Y(IJ),XC(IJ),YC(IJ),U(IJ),V(IJ),P(IJ),T(IJ)
     *,PSI(IJ)
         ENDDO
      ENDDO 
      ENDDO
      CLOSE(7)
C
      OPEN(UNIT=8,FILE='tempi.plt')
      DO J=1,NJ(1)
         IJ=LB(1)+LI(1,NI(1))+J
         WRITE(8,*)YC(IJ),T(IJ)
      ENDDO
      CLOSE(8)  
C
      OPEN(UNIT=8,FILE='tempt.plt')
      DO I=1,NI(1)
         IJ=LB(1)+LI(1,I)+NJ(1)
         WRITE(8,*)XC(IJ),T(IJ)
      ENDDO
      CLOSE(8)  
C
      IB=1
      DO J=2,NJM(IB)
         IJ=LB(IB)+LI(IB,NIM(IB))+J
         IJB=IJ+NJ(IB)
         TEMPGRADS(IJ)=(T(IJB)-T(IJ))/(XC(IJB)-XC(IJ))
         TGRADS(IJB)=DTX(IJ)+(DTX(IJ)-DTX(IJ-NJ(IB)))
     **(1.0D0-FX(IJ-NJ(IB)))
      ENDDO
C
      DO I=2,NIM(IB)
         IJ=LB(IB)+LI(IB,I)+NJM(IB)
         IJB=IJ+1
         TEMPGRADT(IJ)=(T(IJB)-T(IJ))/(YC(IJB)-YC(IJ))
         TGRADT(IJB)=DTY(IJ)+(DTY(IJ)-DTY(IJ-1))
     **(1.0D0-FY(IJ-1))
      ENDDO
C
      OPEN(UNIT=8,FILE='heat_side.plt')
      IB=1
      HEAT_SIDE=0.0D0
      AREA=0.0D0
      AV_NU_SIDEWALL=0.0D0
      DO J=2,NJM(IB)
         IJ=LB(IB)+LI(IB,NIM(IB))+J
         IJB=IJ+NJ(IB)
         IJI=IJ-NJ(IB)
         SX=SEX(IJ)
         SY=SEY(IJ)
         ARR=SQRT(SX**2+SY**2+1.0D-20)
         HEAT_SIDE=HEAT_SIDE+TGRADS(IJB)*ARR
         AREA=AREA+ARR 
         SIDEWALL_NU=TGRADS(IJB)/T(IJB)
         AV_NU_SIDEWALL=AV_NU_SIDEWALL+SIDEWALL_NU*ARR
         WRITE(8,*)YC(IJ),TGRADS(IJB)*ARR,SIDEWALL_NU
      ENDDO
      CLOSE(8)
C
      AV_NU_SIDEWALL=AV_NU_SIDEWALL/AREA  
C
      OPEN(UNIT=8,FILE='heat_top.plt')
      IB=1
      HEAT_TOP=0.0D0
      AREA=0.0D0
      AV_NU_TOPWALL=0.0D0
      DO I=2,NIM(IB)
         IJ=LB(IB)+LI(IB,I)+NJM(IB)
         IJB=IJ+1
         IJI=IJ-1
         SX=SNX(IJ)
         SY=SNY(IJ)
         ARR=SQRT(SX**2+SY**2+1.0D-20)
         HEAT_TOP=HEAT_TOP+TGRADT(IJB)*ARR
         AREA=AREA+ARR 
         TOPWALL_NU=TGRADT(IJB)/T(IJB)
         AV_NU_TOPWALL=AV_NU_TOPWALL+TOPWALL_NU*ARR
         WRITE(8,*)XC(IJ),TGRADT(IJB)*ARR,TOPWALL_NU
      ENDDO
      CLOSE(8)
C
      AV_NU_TOPWALL=AV_NU_TOPWALL/AREA
C
      OPEN(UNIT=10,FILE='data.plt')
      WRITE(10,*)'HEAT_SIDE =',HEAT_SIDE,' HEAT_TOP =',HEAT_TOP
      WRITE(10,*)'AV_NU_SIDEWALL =',AV_NU_SIDEWALL
      WRITE(10,*)'AV_NU_TOPWALL =',AV_NU_TOPWALL
      CLOSE(10)
C 
      PRINT*,'HEAT_SIDE =',HEAT_SIDE,' HEAT_TOP =',HEAT_TOP
      PRINT*,'AV_NU_SIDE AV_NU_TOP =',AV_NU_SIDEWALL,AV_NU_TOPWALL    
C
      AV_TEMP=0.0D0
      VOLUME=0.0D0
      IB=1
      DO I=2,NIM(IB)
         DO J=2,NJM(IB)
            IJ=LB(IB)+LI(IB,I)+J
            VOLUME=VOLUME+VOL(IJ)
            AV_TEMP=AV_TEMP+T(IJ)*VOL(IJ)
         ENDDO
      ENDDO
      AV_TEMP=AV_TEMP/VOLUME
      PRINT*,'AV_TEMP =',AV_TEMP 
C 
      STOP
      END
