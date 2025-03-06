      BLOCK DATA
C
C   MAKES DEFAULT PARAMETER ASSIGNMENTS FOR VEGAS
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
         COMMON/BVEG1/NCALL,ITMX,NPRN,NDEV,XL(10),XU(10),ACC
         COMMON/BVEG2/IT,NDO,SI,SWGT,SCHI,XI(50,10)
         COMMON/BVEG3/ALPH,NDMX,MDS
         DATA NCALL/4999/,ITMX/5/,NPRN/5/,ACC/-1./,
     1        XL/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./,
     2        XU/1.,1.,1.,1.,1.,1.,1.,1.,1.,1./,
     3        ALPH/1.5/,NDMX/50/,MDS/1/,NDEV/6/,
     4        NDO/1/,XI/500*1./,IT/0/,SI,SWGT,SCHI/3*0./
      END
C
      SUBROUTINE VEGAS(NDIM,FXN,AVGI,SD,CHI2A)
C
C   SUBROUTINE PERFORMS NDIM-DIMENSIONAL MONTE CARLO INTEG'N
C      - BY G.P. LEPAGE    SEPT 1976/(REV)AUG 1979
C      - ALGORITHM DESCRIBED IN J COMP PHYS 27,192(1978)
C
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
         COMMON/BVEG1/NCALL,ITMX,NPRN,NDEV,XL(10),XU(10),ACC
         COMMON/BVEG2/IT,NDO,SI,SWGT,SCHI,XI(50,10)
         COMMON/BVEG3/ALPH,NDMX,MDS
         COMMON/BVEG4/CALLS,TI,TSI
         DIMENSION D(50,10),DI(50,10),XIN(50),R(50),DX(10),IA(10),
     1             KG(10),DT(10),X(10)
         REAL RAND(10)
         DATA ONE/1./
         SQRT(A)=DSQRT(A)
         ALOG(A)=DLOG(A)
         ABS(A)=DABS(A)
C
         NDO=1
         DO 1 J=1,NDIM
    1    XI(1,J)=ONE
C
       ENTRY VEGAS1(NDIM,FXN,AVGI,SD,CHI2A)
C         - INITIALIZES CUMULATIVE VARIABLES, BUT NOT GRID
         IT=0
         SI=0.
         SWGT=SI
         SCHI=SI
C
       ENTRY VEGAS2(NDIM,FXN,AVGI,SD,CHI2A)
C         - NO INITIALIZATION
         ND=NDMX
         NG=1
         IF(MDS.EQ.0) GO TO 2
         NG=(NCALL/2.)**(1./NDIM)
         MDS=1
         IF((2*NG-NDMX).LT.0) GO TO 2
         MDS=-1
         NPG=NG/NDMX+1
         ND=NG/NPG
         NG=NPG*ND
    2    K=NG**NDIM
         NPG=NCALL/K
         IF(NPG.LT.2) NPG=2
         CALLS=NPG*K
         DXG=ONE/NG
         DV2G=(CALLS*DXG**NDIM)**2/NPG/NPG/(NPG-ONE)
         XND=ND
         NDM=ND-1
         DXG=DXG*XND
         XJAC=ONE/CALLS
         DO 3 J=1,NDIM
            DX(J)=XU(J)-XL(J)
    3    XJAC=XJAC*DX(J)
C
C   REBIN, PRESERVING BIN DENSITY
         IF(ND.EQ.NDO) GO TO 8
         RC=NDO/XND
         DO 7 J=1, NDIM
            K=0
            XN=0.
            DR=XN
            I=K
    4       K=K+1
            DR=DR+ONE
            XO=XN
            XN=XI(K,J)
    5       IF(RC.GT.DR) GO TO 4
            I=I+1
            DR=DR-RC
            XIN(I)=XN-(XN-XO)*DR
            IF(I.LT.NDM) GO TO 5
            DO 6 I=1,NDM
    6       XI(I,J)=XIN(I)
    7    XI(ND,J)=ONE
         NDO=ND
C
    8    IF(NPRN.GE.0) WRITE(NDEV,200) NDIM,CALLS,IT,ITMX,ACC,NPRN,
     1                       ALPH,MDS,ND,(XL(J),XU(J),J=1,NDIM)
C
       ENTRY VEGAS3(NDIM,FXN,AVGI,SD,CHI2A)
C         - MAIN INTEGRATION LOOP
    9    IT=IT+1
         TI=0.
         TSI=TI
         DO 10 J=1,NDIM
            KG(J)=1
            DO 10 I=1,ND
               D(I,J)=TI
   10    DI(I,J)=TI
C
   11    FB=0.
         F2B=FB
         K=0
   12    K=K+1
         CALL RANDA(NDIM,RAND)
         WGT=XJAC
         DO 15 J=1,NDIM
            XN=(KG(J)-RAND(J))*DXG+ONE
            IA(J)=XN
            IF(IA(J).GT.1) GO TO 13
            XO=XI(IA(J),J)
            RC=(XN-IA(J))*XO
            GO TO 14
   13       XO=XI(IA(J),J)-XI(IA(J)-1,J)
            RC=XI(IA(J)-1,J)+(XN-IA(J))*XO
   14       X(J)=XL(J)+RC*DX(J)
   15    WGT=WGT*XO*XND
C
         F=WGT
         F=F*FXN(X,WGT)
         F2=F*F
         FB=FB+F
         F2B=F2B+F2
         DO 16 J=1,NDIM
            DI(IA(J),J)=DI(IA(J),J)+F
   16    IF(MDS.GE.0) D(IA(J),J)=D(IA(J),J)+F2
         IF(K.LT.NPG) GO TO 12
C
         F2B=SQRT(F2B*NPG)
         F2B=(F2B-FB)*(F2B+FB)
         TI=TI+FB
         TSI=TSI+F2B
         IF(MDS.GE.0) GO TO 18
         DO 17 J=1,NDIM
   17    D(IA(J),J)=D(IA(J),J)+F2B
   18    K=NDIM
   19    KG(K)=MOD(KG(K),NG)+1
         IF(KG(K).NE.1) GO TO 11
         K=K-1
         IF(K.GT.0) GO TO 19
C
C   COMPUTE FINAL RESULTS FOR THIS ITERATION
         TSI=TSI*DV2G
         TI2=TI*TI
         WGT=ONE/TSI
         SI=SI+TI*WGT
         SWGT=SWGT+WGT
         SCHI=SCHI+TI2*WGT
         AVGI=SI/SWGT
         CHI2A=(SCHI-SI*AVGI)/(IT-.9999)
         SD=SQRT(ONE/SWGT)
C
         IF(NPRN.LT.0) GO TO 21
         TSI=SQRT(TSI)
         WRITE(NDEV,201) IT,TI,TSI,AVGI,SD,CHI2A
         IF(NPRN.EQ.0) GO TO 21
         DO 20 J=1,NDIM
   20    WRITE(NDEV,202) J,(XI(I,J),DI(I,J),I=1+NPRN/2,ND,NPRN)
C
C   REFINE GRID
   21    DO 23 J=1,NDIM
            XO=D(1,J)
            XN=D(2,J)
            D(1,J)=(XO+XN)/2.
            DT(J)=D(1,J)
            DO 22 I=2,NDM
               D(I,J)=XO+XN
               XO=XN
               XN=D(I+1,J)
               D(I,J)=(D(I,J)+XN)/3.
   22       DT(J)=DT(J)+D(I,J)
            D(ND,J)=(XN+XO)/2.
   23    DT(J)=DT(J)+D(ND,J)
C
         DO 28 J=1,NDIM
            RC=0.
            DO 24 I=1,ND
               R(I)=0.
               IF(D(I,J).LE.0.) GO TO 24
               XO=DT(J)/D(I,J)
               R(I)=((XO-ONE)/XO/ALOG(XO))**ALPH
   24       RC=RC+R(I)
            RC=RC/XND
            K=0
            XN=0.
            DR=XN
            I=K
   25       K=K+1
            DR=DR+R(K)
            XO=XN
            XN=XI(K,J)
   26       IF(RC.GT.DR) GO TO 25
            I=I+1
            DR=DR-RC
            XIN(I)=XN-(XN-XO)*DR/R(K)
            IF(I.LT.NDM) GO TO 26
            DO 27 I=1,NDM
   27       XI(I,J)=XIN(I)
   28    XI(ND,J)=ONE
C
         IF(IT.LT.ITMX.AND.ACC*ABS(AVGI).LT.SD) GO TO 9
  200    FORMAT(/35H INPUT PARAMETERS FOR VEGAS:  NDIM=I3,8H  NCALL=F8.0
     1     /28X,5H  IT=I5,7H  ITMX=I5/28X,6H  ACC=G9.3
     2     /28X,7H  NPRN=I3,7H  ALPH=F5.2/28X,6H  MDS=I3,6H   ND=I4
     3     /28X,10H  (XL,XU)=(T40,2H( G12.6,3H , G12.6,2H )))
  201    FORMAT(///21H INTEGRATION BY VEGAS//14H ITERATION NO.I3,
     1     14H:   INTEGRAL =G14.8/21X,10HSTD DEV  =G10.4/
     2     34H ACCUMULATED RESULTS:   INTEGRAL =G14.8/
     3     24X,10HSTD DEV  =G10.4/24X,17HCHI**2 PER it'n =G10.4)
C     the apostrophe(') in the line above is needed for proper output
  202    FORMAT(/15H DATA FOR AXIS I2/25H    X       DELTA I       ,
     1     24H   X       DELTA I       ,18H   X       DELTA I
     2     /(1H F7.6,1X,G11.4,5X,F7.6,1X,G11.4,5X,F7.6,1X,G11.4))
         RETURN
      END

      SUBROUTINE RANDA(N,RAND)
         DIMENSION RAND(N)
         CALL RANDOM_NUMBER(RAND)
         RETURN
      END