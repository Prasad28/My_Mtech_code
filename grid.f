      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NX=502,NY=502,NXY=NX*NY,MB=3)
      COMMON LI(MB,NX),LB(MB),NI(MB),NJ(MB),NIM(MB),NJM(MB),NIJ(MB)
      COMMON FD(NX),X(NXY),Y(NXY),XOO(NXY),YOO(NXY)

      DIMENSION NICV(MB),NJCV(MB)

      CHARACTER FILIN*12,FILOUT*12
C
      FILIN='grid.inp'
      OPEN (UNIT=5,FILE=FILIN)
      REWIND 5
C
      PRINT*,'ENTER NUMER OF BLOCK'
      READ(5,*) NB
      DO I=1,NB
         PRINT*,'DATA FOR BLOCK    ',I
         PRINT*,'ENTER THE NUMBER OF CONTROL VOLUME IN X-DIRECTION'
         READ(5,*)NICV(I)
         PRINT*,'ENTER THE NUMBER OF CONTROL VOLUME IN Y-DIRECTION'
         READ(5,*)NJCV(I)
      ENDDO
C
C.....EXPANSION FACTOR CALCULATION
      PI=4.0D0*ATAN(1.0D0)
C
      DO IB=1,NB
         NI(IB)=NICV(IB)+2
         NJ(IB)=NJCV(IB)+2 
         NIM(IB)=NI(IB)-1
         NJM(IB)=NJ(IB)-1   
         NIJ(IB)=NI(IB)*NJ(IB)
      ENDDO 
C
      DO IB=1,NB
      DO I=1,NI(IB)
         LI(IB,I)=(I-1)*NJ(IB) 
      ENDDO
      ENDDO
C
      LB(1)=0
      DO IB=2,NB
         LB(IB)=LB(IB-1)+NI(IB-1)*NJ(IB-1)
      ENDDO
C...
C
      DSIZERX=1.0D0
      DSIZERY=1.0D0
      WALL_THICKX=0.1D0
      WALL_THICKY=0.1D0 
C
      EXPW1=1.03D0
      EXPS1=1.03D0
      EXPSW1=1.1D0
      EXPW3=1.1D0
C
      DO IB=1,NB
C
      IF(IB.EQ.1)THEN  
C.....WEST CO-ORDINATES
      DLW=0.5D0*DSIZERY  
      XLS=0.0D0
      YLS=0.0D0
      XLE=0.0D0
      YLE=0.5D0*DSIZERY
      EXPW=EXPW1
      NJCV1=NJCV(IB)/2
      CALL DIVLINE(EXPW,NJCV1,DLW,XLS,YLS,XLE,YLE)
      DO I=1,NJCV1+1
         IJ=LB(IB)+LI(IB,1)+I
         XOO(IJ)=XLS+FD(I)*(XLE-XLS)
         YOO(IJ)=YLS+FD(I)*(YLE-YLS)
      END DO
C
      DLW=DSIZERY-0.5D0*DSIZERY  
      XLS=0.0D0
      YLS=0.5D0*DSIZERY
      XLE=0.0D0
      YLE=DSIZERY
      EXPW=2.0D0-EXPW1
      NJCV2=NJCV(IB)-NJCV1
      CALL DIVLINE(EXPW,NJCV2,DLW,XLS,YLS,XLE,YLE)
      DO I=NJCV1+2,NJCV1+NJCV2+1
         IJ=LB(IB)+LI(IB,1)+I
         XOO(IJ)=XLS+FD(I-NJCV1)*(XLE-XLS)
         YOO(IJ)=YLS+FD(I-NJCV1)*(YLE-YLS)
      END DO
C
C.....EAST CO-ORDINATES
      DO I=1,NJCV(IB)+1
         IJ=LB(IB)+LI(IB,NIM(IB))+I
         IJW=LB(IB)+LI(IB,1)+I
         XOO(IJ)=DSIZERX
         YOO(IJ)=YOO(IJW)
      END DO
C
C.....SOUTH CO-ORDINATE
      DLS=0.5D0*DSIZERX
      XLS=0.0D0
      YLS=0.0D0
      XLE=0.5D0*DSIZERX
      YLE=0.0D0
      EXPS=EXPS1
      NICV1=NICV(IB)/2
      CALL DIVLINE(EXPS,NICV1,DLS,XLS,YLS,XLE,YLE)
      DO I=1,NICV1+1
         IJ=LB(IB)+LI(IB,I)+1
         XOO(IJ)=XLS+FD(I)*(XLE-XLS)
         YOO(IJ)=YLS+FD(I)*(YLE-YLS)
      END DO
C
      DLS=DSIZERX-0.5D0*DSIZERX
      XLS=0.5D0*DSIZERX
      YLS=0.0D0
      XLE=DSIZERX
      YLE=0.0D0
      EXPS=2.0D0-EXPS1
      NICV2=NICV(IB)-NICV1
      CALL DIVLINE(EXPS,NICV2,DLS,XLS,YLS,XLE,YLE)
      DO I=NICV1+2,NICV1+NICV2+1
         IJ=LB(IB)+LI(IB,I)+1
         XOO(IJ)=XLS+FD(I-NICV1)*(XLE-XLS)
         YOO(IJ)=YLS+FD(I-NICV1)*(YLE-YLS)
      END DO
C
C.....NORTH CO-ORDINATE
      DO I=1,NICV(IB)+1
         IJ=LB(IB)+LI(IB,I)+NJM(IB)
         IJS=LB(IB)+LI(IB,I)+1
         XOO(IJ)=XOO(IJS)
         YOO(IJ)=DSIZERY
      END DO
C
CCCCCC
      DO I=1,NI(IB)
         IJ=LB(IB)+LI(IB,I)+NJ(IB)
         XOO(IJ)=XOO(IJ-1)
         YOO(IJ)=YOO(IJ-1)
      ENDDO
C
      DO J=1,NJ(IB)
         IJ=LB(IB)+LI(IB,NI(IB))+J
         XOO(IJ)=XOO(IJ-NJ(IB))
         YOO(IJ)=YOO(IJ-NJ(IB))
      ENDDO
C
      CALL CALXY(0,IB)
C
      OPEN(UNIT=4,FILE='coor1.dat')
      WRITE(4,*)'TITLE = "Example: Multi-Zone 2-D Plot"'
      WRITE(4,*)'VARIABLES = "X", "Y"'
      WRITE(FILOUT,'(a5,I1)') 'BLOCK',IB
      WRITE(4,*)'ZONE T=', FILOUT,' ','I=',NJ(IB),' ','J=',NI(IB),' ',
     *'DATAPACKING=POINT'
      DO I=1,NI(IB)
         DO J=1,NJ(IB)
            IJ=LB(IB)+LI(IB,I)+J
            WRITE(4,*)XOO(IJ),YOO(IJ)
         ENDDO
      ENDDO 
      CLOSE(4)
C
      ENDIF 
C
      IF(IB.EQ.2)THEN  
C.....WEST CO-ORDINATES
      DO I=1,NJCV(IB)+1
         IJ=LB(IB)+LI(IB,1)+I
         IJ1=LB(IB-1)+LI(IB-1,NIM(IB-1))+I
         XOO(IJ)=XOO(IJ1)
         YOO(IJ)=YOO(IJ1)
      END DO
C
C.....EAST CO-ORDINATES
      DO I=1,NJCV(IB)+1
         IJ=LB(IB)+LI(IB,NIM(IB))+I
         IJW=LB(IB)+LI(IB,1)+I
         XOO(IJ)=DSIZERX+WALL_THICKX
         YOO(IJ)=YOO(IJW)
      END DO
C
C.....SOUTH CO-ORDINATE
      DLS=0.5D0*WALL_THICKX
      XLS=DSIZERX
      YLS=0.0D0
      XLE=DSIZERX+0.5D0*WALL_THICKX
      YLE=0.0D0
      EXPS=EXPSW1
      NICV1=NICV(IB)/2 
      CALL DIVLINE(EXPS,NICV1,DLS,XLS,YLS,XLE,YLE)
      DO I=1,NICV1+1
         IJ=LB(IB)+LI(IB,I)+1
         XOO(IJ)=XLS+FD(I)*(XLE-XLS)
         YOO(IJ)=YLS+FD(I)*(YLE-YLS)
      END DO
C
      DLS=0.5D0*WALL_THICKX
      XLS=DSIZERX+0.5D0*WALL_THICKX
      YLS=0.0D0
      XLE=DSIZERX+WALL_THICKX
      YLE=0.0D0
      EXPS=2.0D0-EXPSW1
      NICV2=NICV(IB)-NICV1 
      CALL DIVLINE(EXPS,NICV2,DLS,XLS,YLS,XLE,YLE)
      DO I=NICV1+1,NICV1+NICV2+1
         IJ=LB(IB)+LI(IB,I)+1
         XOO(IJ)=XLS+FD(I-NICV1)*(XLE-XLS)
         YOO(IJ)=YLS+FD(I-NICV1)*(YLE-YLS)
      END DO
C
C.....NORTH CO-ORDINATE
      DO I=1,NICV(IB)+1
         IJ=LB(IB)+LI(IB,I)+NJM(IB)
         IJS=LB(IB)+LI(IB,I)+1
         XOO(IJ)=XOO(IJS)
         YOO(IJ)=DSIZERY
      END DO
C
CCCCCC
      DO I=1,NI(IB)
         IJ=LB(IB)+LI(IB,I)+NJ(IB)
         XOO(IJ)=XOO(IJ-1)
         YOO(IJ)=YOO(IJ-1)
      ENDDO
C
      DO J=1,NJ(IB)
         IJ=LB(IB)+LI(IB,NI(IB))+J
         XOO(IJ)=XOO(IJ-NJ(IB))
         YOO(IJ)=YOO(IJ-NJ(IB))
      ENDDO
C
      CALL CALXY(0,IB)
C
      OPEN(UNIT=4,FILE='coor2.dat')
      WRITE(4,*)'TITLE = "Example: Multi-Zone 2-D Plot"'
      WRITE(4,*)'VARIABLES = "X", "Y"'
      WRITE(FILOUT,'(a5,I1)') 'BLOCK',IB
      WRITE(4,*)'ZONE T=', FILOUT,' ','I=',NJ(IB),' ','J=',NI(IB),' ',
     *'DATAPACKING=POINT'
      DO I=1,NI(IB)
         DO J=1,NJ(IB)
            IJ=LB(IB)+LI(IB,I)+J
            WRITE(4,*)XOO(IJ),YOO(IJ)
         ENDDO
      ENDDO 
      CLOSE(4)
C
      ENDIF 
C
      IF(IB.EQ.3)THEN  
C.....WEST CO-ORDINATES
      DLW=0.5D0*WALL_THICKY  
      XLS=0.0D0
      YLS=DSIZERY
      XLE=0.0D0
      YLE=DSIZERY+0.5D0*WALL_THICKY
      EXPW=EXPW3
      NJCV1=NJCV(IB)/2
      CALL DIVLINE(EXPW,NJCV1,DLW,XLS,YLS,XLE,YLE)
      DO I=1,NJCV1+1
         IJ=LB(IB)+LI(IB,1)+I
         XOO(IJ)=XLS+FD(I)*(XLE-XLS)
         YOO(IJ)=YLS+FD(I)*(YLE-YLS)
      END DO
C
      DLW=WALL_THICKY-0.5D0*WALL_THICKY  
      XLS=0.0D0
      YLS=DSIZERY+0.5D0*WALL_THICKY
      XLE=0.0D0
      YLE=DSIZERY+WALL_THICKY
      EXPW=2.0D0-EXPW1
      NJCV2=NJCV(IB)-NJCV1
      CALL DIVLINE(EXPW,NJCV2,DLW,XLS,YLS,XLE,YLE)
      DO I=NJCV1+2,NJCV1+NJCV2+1
         IJ=LB(IB)+LI(IB,1)+I
         XOO(IJ)=XLS+FD(I-NJCV1)*(XLE-XLS)
         YOO(IJ)=YLS+FD(I-NJCV1)*(YLE-YLS)
      END DO
C
C.....EAST CO-ORDINATES
      DO I=1,NJCV(IB)+1
         IJ=LB(IB)+LI(IB,NIM(IB))+I
         IJW=LB(IB)+LI(IB,1)+I
         XOO(IJ)=DSIZERX+WALL_THICKX
         YOO(IJ)=YOO(IJW)
      END DO
C
C.....SOUTH CO-ORDINATE
      DO I=1,NICV(IB-2)+1
         IJ=LB(IB)+LI(IB,I)+1
         IJ1=LB(IB-2)+LI(IB-2,I)+NJM(IB-2)
         XOO(IJ)=XOO(IJ1)
         YOO(IJ)=DSIZERY
      END DO
C
      DO I=NICV(IB-2)+2,NICV(IB)+1
         IJ=LB(IB)+LI(IB,I)+1
         IJ2=LB(IB-1)+LI(IB-1,I-NICV(IB-2))+NJM(IB-1)
         XOO(IJ)=XOO(IJ2)
         YOO(IJ)=DSIZERY
      END DO
C
C.....NORTH CO-ORDINATE
      DO I=1,NICV(IB)+1
         IJ=LB(IB)+LI(IB,I)+NJM(IB)
         IJS=LB(IB)+LI(IB,I)+1
         XOO(IJ)=XOO(IJS)
         YOO(IJ)=DSIZERY+WALL_THICKY
      END DO
C
CCCCCC
      DO I=1,NI(IB)
         IJ=LB(IB)+LI(IB,I)+NJ(IB)
         XOO(IJ)=XOO(IJ-1)
         YOO(IJ)=YOO(IJ-1)
      ENDDO
C
      DO J=1,NJ(IB)
         IJ=LB(IB)+LI(IB,NI(IB))+J
         XOO(IJ)=XOO(IJ-NJ(IB))
         YOO(IJ)=YOO(IJ-NJ(IB))
      ENDDO
C
      CALL CALXY(0,IB)
C
      OPEN(UNIT=4,FILE='coor3.dat')
      WRITE(4,*)'TITLE = "Example: Multi-Zone 2-D Plot"'
      WRITE(4,*)'VARIABLES = "X", "Y"'
      WRITE(FILOUT,'(a5,I1)') 'BLOCK',IB
      WRITE(4,*)'ZONE T=', FILOUT,' ','I=',NJ(IB),' ','J=',NI(IB),' ',
     *'DATAPACKING=POINT'
      DO I=1,NI(IB)
         DO J=1,NJ(IB)
            IJ=LB(IB)+LI(IB,I)+J
            WRITE(4,*)XOO(IJ),YOO(IJ)
         ENDDO
      ENDDO 
      CLOSE(4)
C
      ENDIF 
C
      ENDDO
C
      NIJS=0
      DO I=1,NB
         NIJS=NIJS+NIJ(I)
      ENDDO
      PRINT*,'NO. OF GRID POINTS = ',NIJS
C
      OPEN(UNIT=6,FILE='2dgrid',FORM='UNFORMATTED')  
      WRITE(6)NB,NIJS
      DO IB=1,NB
         WRITE(6)NI(IB),NJ(IB),NIM(IB),NJM(IB),NIJ(IB)
      ENDDO
      WRITE(6)(XOO(IJ),IJ=1,NIJS),(YOO(IJ),IJ=1,NIJS)
      CLOSE(6)
C.....
C
      STOP
      END       
C
C###########################################################
      SUBROUTINE CALXY(IDD,IB) 
C###########################################################
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NX=502,NY=502,NXY=NX*NY,MB=3)
      COMMON LI(MB,NX),LB(MB),NI(MB),NJ(MB),NIM(MB),NJM(MB),NIJ(MB)
      COMMON FD(NX),X(NXY),Y(NXY),XOO(NXY),YOO(NXY)

      DIMENSION FL(NXY),FR(NXY)
C
C.....OFFSETS AT WEST AND EAST BOUNDARY
C
      IIE=LB(IB)+LI(IB,NIM(IB))
      IIW=LB(IB)+LI(IB,1)
C
C==================================================================
C.....CALCULATE COORDINATES OF INTERIOR POINTS (STRAIGHT LINES N-S)
C==================================================================
C
      IF(IDD.EQ.0) THEN
C
C.....DISTANCE BETWEEN CORNER POINTS (NORTH - SOUTH)
C
        DLR=1.0D0/SQRT((XOO(IIE+NJM(IB))-XOO(IIE+1))**2+
     *(YOO(IIE+NJM(IB))-YOO(IIE+1))**2)
        DLL=1.0D0/SQRT((XOO(IIW+NJM(IB))-XOO(IIW+1))**2+
     *(YOO(IIW+NJM(IB))-YOO(IIW+1))**2)
C
C.....DISTRIBUTION FUNCTION ON WEST AND EAST SIDE
C
        DO J=1,NJM(IB)
          FL(J)=SQRT((XOO(IIW+J)-XOO(IIW+1))**2+
     *(YOO(IIW+J)-YOO(IIW+1))**2)*DLL
          FR(J)=SQRT((XOO(IIE+J)-XOO(IIE+1))**2+
     *(YOO(IIE+J)-YOO(IIE+1))**2)*DLR
        END DO
C
C.....DISTANCE BETWEEN OPOSITE POINTS ON SOUTH AND NORTH SIDE
C
        RNIM=1.0D0/REAL(NIM(IB)-1)
        DO I=2,NIM(IB)-1
          II=LB(IB)+LI(IB,I)
          DX=XOO(II+NJM(IB))-XOO(II+1)
          DY=YOO(II+NJM(IB))-YOO(II+1)
C
C.....DISTRIBUTE POINTS USING INTERPOLATED WEST AND EAST DISTRIBUTIONS
C
          DO J=2,NJM(IB)-1
            FAC=(REAL(I-1)*FR(J)+REAL(NIM(IB)-I)*FL(J))*RNIM
            XOO(II+J)=XOO(II+1)+FAC*DX
            YOO(II+J)=YOO(II+1)+FAC*DY
          END DO
        END DO
C
C==================================================================
C.....CALCULATE COORDINATES OF INTERIOR POINTS (STRAIGHT LINES E-W)
C==================================================================
C
      ELSEIF(IDD.EQ.1) THEN
C
C.....DISTANCE BETWEEN CORNERS (EAST - WEST)
C
        DLR=1.0D0/SQRT((XOO(IIE+NJM(IB))-XOO(IIW+NJM(IB)))**2+
     *              (YOO(IIE+NJM(IB))-YOO(IIW+NJM(IB)))**2)
        DLL=1.0D0/SQRT((XOO(IIE+1)-XOO(IIW+1))**2+
     *(YOO(IIE+1)-YOO(IIW+1))**2)
C
C.....DISTRIBUTION FUNCTION ALONG SOUTH AND NORTH BOUNDARY
C
        DO I=1,NIM(IB)
          II=LB(IB)+LI(IB,I)
          FL(I)=SQRT((XOO(II+1)-XOO(IIW+1))**2+
     *(YOO(II+1)-YOO(IIW+1))**2)*DLL
          FR(I)=SQRT((XOO(II+NJM(IB))-XOO(IIW+NJM(IB)))**2+
     *            (YOO(II+NJM(IB))-YOO(IIW+NJM(IB)))**2)*DLR
        END DO
C
C.....DISTANCE BETWEEN OPOSITE POINTS ON EAST AND WEST BOUNDARY
C
        RNJM=1.0D0/REAL(NJM(IB)-1)
        DO J=2,NJM(IB)-1
          DX=XOO(IIE+J)-XOO(IIW+J)
          DY=YOO(IIE+J)-YOO(IIW+J)
C
C.....DISTRIBUTE POINTS USING INTERPOLATED NORTH AND SOUTH DISTRIBUTIONS
C
          DO I=2,NIM(IB)-1
            IJ=LB(IB)+LI(IB,I)+J
            FAC=(REAL(J-1)*FR(I)+REAL(NJM(IB)-J)*FL(I))*RNJM
            XOO(IJ)=XOO(IIW+J)+FAC*DX
            YOO(IJ)=YOO(IIW+J)+FAC*DY
          END DO
        END DO
C
      ENDIF
C
      RETURN
      END
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DIVLINE(EXP,NSEG,XYL,XLS,YLS,XLE,YLE)
CCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NX=502,NY=502,NXY=NX*NY,MB=3)
      COMMON LI(MB,NX),LB(MB),NI(MB),NJ(MB),NIM(MB),NJM(MB),NIJ(MB)
      COMMON FD(NX),X(NXY),Y(NXY),XOO(NXY),YOO(NXY)

      DIMENSION DL(NX)

        IF(ABS(EXP-1.0D0).LT.0.001D0) THEN
          DX1=XYL/REAL(NSEG)
        ELSE
          DX1=XYL*(1.0D0-EXP)/(1.0D0-EXP**NSEG)
        ENDIF
C
C.....CALCULATE SIZE OF INTERVALS AND SCALE TO THE RANGE (0 ... 1)
C
      XYLR=1.0D0/(XYL+1.0E-20)
      DL(1)=DX1
      FD(1)=0.0D0
      SUML=0.0D0
C
      DO I=2,NSEG
        DL(I)=DL(I-1)*EXP
        SUML=SUML+DL(I-1)
        FD(I)=SUML*XYLR
      END DO
      FD(NSEG+1)=1.0D0 
C
C
      RETURN
      END
C

