      SUBROUTINE POSTAAR(SMIN,SMAX,TOL,NH,NR,NHO,NSO,IX,IY,
     1      POSTS,PX1,PY1,PX2,PY2,PX3,PY3,PX4,PY4)
      implicit none

C***************************************************************************
C
C  POSTHOLE
C
C  PURPOSE:
C  FINDS RECTANGLES FORMED OF FOUR POSTHOLES FROM THE LIST
C
C  OPERATION:
C  TAKES INPUT AND OUTPUT FILE NAMES AND ASKS FOR THREE PARAMETERS:
C  UPPER AND LOWER BOUNDS FOR THE RECTANGLE SIDES, AND THE TOLERANCE
C  OF POSTHOLE POSITION
C  GRAPHIC OUTPUT ASSUMES VGA COLOR DISPLAY
C
C  FILES:
C  LUN 3      INPUT, LIST OF POSTHOLES
C  LUN 2      OUTPUT, LIST OF RECTANGLES
C  LUN 7      SCRATCH
C
C  SUBROUTINES:
C  ALL SUBROUTINES USED ARE INCLUDED IN THE SOURCE CODE
C
C  GENERATION:
C  fl POSTHOLE.FOR /link GRAPHICS
C
C****************************************************************************
      
      INTEGER, INTENT (IN) :: SMIN,SMAX,TOL,NHO,IX,IY,NSO
      INTEGER, INTENT (INOUT) :: NR, NH
      INTEGER, INTENT (OUT), DIMENSION(NHO) :: POSTS,PX1,PY1,PX2,PY2,
     1      PX3,PY3,PX4,PY4
      INTEGER :: I,IB,IE,IH1,IH2,IH3,IH4,III,IJB1,
     1      IJB2,IUP,IW1,J,JJJ,K,L,LX,LY,NHI,NHIJ,NHIJ1,
     1      NHIJ2,IHUGE,IP,IP1,IWI,IXMA,IXMI,IYMA,IYMI,
     1      JUP
      REAL :: XSIZE,YSIZE,PI,XSQ,YSQ
      PARAMETER(PI=3.14159265)
      REAL, DIMENSION (NHO) :: WX,WY
      DIMENSION IP(0:NSO,NSO),IP1(0:NSO)
      CHARACTER(LEN=6) NAME,CW
      DIMENSION IX(NHO),IY(NHO),IW1(NHO),IWI(NHO),
     1      NAME(NHO),CW(NHO)
          

      NH=NHO
      
C ---  FIND THE SCALE FOR COORDINATES
C
C --- HUGE(I) IS THE MAXIMUM INTEGER
      IHUGE=HUGE(I)-1
      IXMI=IHUGE
      IYMI=IHUGE
      IXMA=-IHUGE
      IYMA=-IHUGE
      DO I=1,NH
      IXMI=MIN0(IXMI,IX(I))
      IYMI=MIN0(IYMI,IY(I))
      IXMA=MAX0(IXMA,IX(I))
      IYMA=MAX0(IYMA,IY(I))
      ENDDO
      XSIZE=IXMA-IXMI
      YSIZE=IYMA-IYMI


C --- DIVIDE THE WHOLE AREA INTO LY BY LX (APPROXIMATE) SQUARES
C --- OF SIZE SMAX BY SMAX EACH
C
C      CALL GETTIM(IHR,IMIN,ISEC,I100TH)
      LX=MAX0(1,IFIX(XSIZE/SMAX))
      LY=MAX0(1,IFIX(YSIZE/SMAX))
      XSQ=XSIZE/LX
      YSQ=YSIZE/LY      
      GO TO 60



C --- SORT THE ARRAY IX, AND ORDER IY AND NAME, TOO
C
60      CALL SORT3(NH,IX,IY,NAME,IW1,CW,IWI,NHO)
C
C --- CONSTRUCT THE COLUMN POINTER ARRAY IP1
C --- ENTRY IP1(JX-1) WILL POINT TO THE PLACE WHERE THE HOLES
C --- OF JXTH COLUMN OF SQUARES BEGIN IN THE ARRAYS IX, IY, AND NAME
C
      IP1(0)=0
      I=1
      K=1
70      IF(IX(K).LE.IXMI+INT(I*XSQ))THEN
          K=K+1
          IF(K.GT.NH)THEN
            DO III=I,LX
            IP1(III)=NH
            ENDDO
          ELSE
            GO TO 70
          ENDIF
      ELSE
          IP1(I)=K-1
          I=I+1
          IF(I.LE.LX)GO TO 70
      ENDIF
      
C
C --- CONSTRUCT A LY BY LX POINTER ARRAY IP
C --- ENTRY IP(JY-1,JX-1) WILL POINT TO THE PLACE WHERE THE HOLES
C --- OF SQUARE (JY,JX) BEGIN IN IX, IY, AND NAME
C
C --- TRIVIAL INITIAL VALUES
      DO I=1,LX
      IP(0,I)=IP1(I-1)
      ENDDO
C --- LOOP FOR THE ITH COLUMN OF IP
      DO I=1,LX
      IB=IP1(I-1)+1
      IE=IP1(I)
      NHI=IE-IB+1
      IF(NHI.EQ.0)THEN
          DO J=1,LY
          IP(J,I)=IE
          ENDDO
      ELSE
          IF(NHI.GT.1)
     1            CALL SORT3(NHI,IY(IB),IX(IB),NAME(IB),IW1,CW,IWI,NHO)
          J=1
          K=IB
80          IF(IY(K).LE.IYMI+INT(J*YSQ))THEN
            K=K+1
            IF(K.GT.IE)THEN
                DO JJJ=J,LY
                IP(JJJ,I)=IE
                ENDDO
            ELSE
                GO TO 80
            ENDIF
          ELSE
            IP(J,I)=K-1
            J=J+1
            IF(J.LE.LY)GO TO 80
          ENDIF
      ENDIF
      ENDDO
      
C
C --- RECTANGLES ARE LOOKED FOR IN THE LOOP OVER ALL
C --- THE LY BY LX SQUARES AND THE FOUND ONES ARE
C --- WRITTEN TO DIRECT ACCESS SCRATCH UNIT 17
C
      OPEN(UNIT=17,STATUS='SCRATCH',
     1 ACCESS='DIRECT',RECL=16,form="unformatted", 
     1 action="readwrite")
      JUP=MAX0(1,LY-1)
      IUP=MAX0(1,LX-1)
      DO 2010 J=1,JUP
      DO 2020 I=1,IUP
      IJB1=IP(J-1,I)+1
      IF(LY.GT.1)THEN
          NHIJ1=IP(J+1,I)-IJB1+1
      ELSE
          NHIJ1=IP(J,I)-IJB1+1
      ENDIF
      NHIJ=NHIJ1
      K=IJB1
      DO L=1,NHIJ1
      WX(L)=IX(K)
      WY(L)=IY(K)
      K=K+1
      ENDDO
      IF(LX.GT.1)THEN
          IJB2=IP(J-1,I+1)+1
          IF(LY.GT.1)THEN
            NHIJ2=IP(J+1,I+1)-IJB2+1
          ELSE
            NHIJ2=IP(J,I+1)-IJB2+1
          ENDIF
          NHIJ=NHIJ1+NHIJ2
          K=IJB2
          DO L=NHIJ1+1,NHIJ
          WX(L)=IX(K)
          WY(L)=IY(K)
          K=K+1
          ENDDO
      ELSE
          IJB2=1
      ENDIF
      IF(NHIJ.GE.4)CALL RECTAN(NHIJ,NHIJ1,IJB1-1,IJB2-1,WX,WY,
     1      NR,SMIN,SMAX,TOL,NHO)

2020  CONTINUE
2010  CONTINUE


C
C --- WRITE THE RESULTING NR RECTANGLES ON UNIT 9
C --- (RECTANGLE OUTPUT UNIT)
C

      DO I=1,NR
      READ(17,REC=I)IH1,IH2,IH3,IH4
      POSTS(I)=I
      PX1(I)=IX(IH1)
      PY1(I)=IY(IH1)
      PX2(I)=IX(IH2)
      PY2(I)=IY(IH2)
      PX3(I)=IX(IH3)
      PY3(I)=IY(IH3)
      PX4(I)=IX(IH4)
      PY4(I)=IY(IH4)
      ENDDO

      CLOSE (17)
      END
C***********************************************************************

C***********************************************************************
C
C   SUBROUTINE SORT3 SORTS AN ARRAY IA OF LENGTH N INTO ASCENDING
C      NUMERICAL ORDER WHILE MAKING THE CORRESPONDING REARRANGEMENTS
C      OF THE ARRAYS IB AND C.
C      AN INDEX TABLE IS CONSTRUCTED VIA THE SUBROUTINE INDEXX.
C   PARAMETERS
C      N IS THE LENGTH OF ALL THE ARRAYS, IA AND IB ARE INPUT ARRAYS,
C      KSP AND IWKSP ARE WORKING ARRAYS.
C      INPUT ARRAY C IS CHARACTER*6, WORKSPACE CWKSP AS WELL.
C      SEE NUMERICAL RECIPES
C
C***********************************************************************
      SUBROUTINE SORT3(N,IA,IB,C,KSP,CWKSP,IWKSP,NHO)
      implicit none
      INTEGER N,IA,IB,KSP,IWKSP,J,NHO
      DIMENSION IA(NHO),IB(NHO),KSP(NHO),IWKSP(NHO)
      CHARACTER*6 C(NHO),CWKSP(NHO)     
C --- MAKE THE INDEX TABLE
      CALL INDEXX(N,IA,IWKSP,NHO)
C --- SAVE THE ARRAY IA
      DO J=1,N
      KSP(J)=IA(J)
      ENDDO
C --- COPY IT BACK IN THE REARRANGED ORDER
      DO J=1,N
      IA(J)=KSP(IWKSP(J))
      ENDDO
C --- SAVE THE ARRAY IB
      DO J=1,N
      KSP(J)=IB(J)
      ENDDO
C --- COPY IT BACK IN THE REARRANGED ORDER
      DO J=1,N
      IB(J)=KSP(IWKSP(J))
      ENDDO
C --- SAVE THE ARRAY C
      DO J=1,N
      CWKSP(J)=C(J)
      ENDDO
C --- COPY IT BACK IN THE REARRANGED ORDER
      DO J=1,N
      C(J)=CWKSP(IWKSP(J))
      ENDDO    
      RETURN
      END
C***********************************************************************
C
C   SUBROUTINE INDEXX INDEXES AN ARRAY IARRIN OF LENGTH N,
C      I.E. OUTPUTS THE ARRAY INDX SUCH THAT IARRIN(INDX(J))
C      IS IN ASCENDING ORDER FOR J=1,...,N.
C   PARAMETERS
C      THE INPUT QUANTITIES N AND IARRIN ARE NOT CHANGED.
C      INDX IS ASSIGNED THE VALUES OF THE CORRESPONDING INDICES.
C      SEE NUMERICAL RECIPES
C
C***********************************************************************
      SUBROUTINE INDEXX(N,IARRIN,INDX,NHO)
      implicit none
      INTEGER N,IARRIN,INDX,J,INDXT,IQ,IR,L,I,NHO
      DIMENSION IARRIN(NHO),INDX(NHO)
C --- INITIALIZATION OF THE INDEX ARRAY
      DO J=1,N
          INDX(J)=J
      ENDDO
C --- HEAPSORT WITH INDIRECT INDEXING THROUGH INDX
C --- IN ALL REFERENCES TO IARRIN
      L=N/2+1
      IR=N
10      CONTINUE
          IF(L.GT.1)THEN
            L=L-1
            INDXT=INDX(L)
            IQ=IARRIN(INDXT)
          ELSE
            INDXT=INDX(IR)
            IQ=IARRIN(INDXT)
            INDX(IR)=INDX(1)
            IR=IR-1
            IF(IR.EQ.1)THEN
                INDX(1)=INDXT
                RETURN
            ENDIF
          ENDIF
          I=L
          J=L+L
20          IF(J.LE.IR)THEN
            IF(J.LT.IR)THEN
                IF(IARRIN(INDX(J)).LT.IARRIN(INDX(J+1)))J=J+1
            ENDIF
            IF(IQ.LT.IARRIN(INDX(J)))THEN
                INDX(I)=INDX(J)
                I=J
                J=J+J
            ELSE
                J=IR+1
            ENDIF
          GO TO 20
          ENDIF
          INDX(I)=INDXT
      GO TO 10
      END
C***********************************************************************
C
C   SUBROUTINE RECTAN FINDS RECTANGLES FORMED OF FOUR POSTHOLES
C      FROM THOSE PRESORTED AND SUPPLIED AS PARAMETERS (I.E. THOSE FROM
C      FOUR BASIC SQUARES)
C      THIS VERSION FINDS THE OTHER TWO VERTICES TO AN EDGE
C      USING RIGHT ANGLES
C   PARAMETERS
C      NH IS THE NUMBER OF POSTHOLES, NH1 IS THEIR NUMBER IN THE
C      LEFT-HAND PAIR OF SQUARES.
C      NB1 AND NB2 ARE THE POSITIONS OF THE POSTHOLES OF THE LEFT-HAND
C      AND RIGHT-HAND PAIRS OF SQUARES IN THE GLOBAL POSTHOLE ARRAY,
C      RESPECTIVELY, AND X, Y ARE ARRAYS OF COORDINATES OF POSTHOLES
C      JUST BEING PROCESSED.
C
C***********************************************************************
      SUBROUTINE RECTAN(NH,NH1,NB1,NB2,X,Y,NR,SMIN,SMAX,TOL,NHO)
      implicit none
C
C --- NR IS THE CURRENT NUMBER OF RECTANGLES FOUND.
C --- SMIN, SMAX, AND TOL ARE THE MINIMUM AND MAXIMUM LENGTHS OF
C --- RECTANGLE EDGES, AND THE TOLERANCE
C
      REAL :: PIH,ALPHA,DH12,DH13,DH14,DH23,DH24,DH34,PHI,
     1      SI3A,SI3I,SI4,SIMAX,SIMIN,TA3,TA4,X,Y
      PARAMETER(PIH=1.5707963)
      DIMENSION X(NHO),Y(NHO)
      LOGICAL SECOND
      INTEGER SMIN,SMAX,TOL,NH,NH1,NB1,NB2,NR,IH1,IH2,IH3,IH4,
     1      IIH1,IIH2,IIH3,IIH4,IG1,IG2,IG3,IG4,IR,NHO

C
C --- MAKE EDGES OF ADMISSIBLE PAIRS OF POSTHOLES
C --- AND LOOK FOR THE OTHER EDGES
C
      DO IH1=1,NH-1
      DO IH2=IH1+1,NH
C --- EDGE INADMISSIBLE
      DH12=SQRT((X(IH1)-X(IH2))**2+(Y(IH1)-Y(IH2))**2)
      IF(SMIN.GT.DH12.OR.DH12.GT.SMAX)GO TO 90
C --- THE ANGLE OF THE EDGE (RAY) IH1,IH2 WITH THE POSITIVE X-AXIS
      PHI=ATAN2(Y(IH2)-Y(IH1),X(IH2)-X(IH1))
C --- ALPHA IS THE ANGLE OF THE SOUGHT LATERAL EDGE.
C --- TWO CASES ARE TO BE CONSIDERED
      ALPHA=PHI-PIH
      SECOND=.TRUE.
      GO TO 20
10      ALPHA=PHI+PIH
      SECOND=.FALSE.
20      SI3I=X(IH1)+SMIN*COS(ALPHA)
      SI3A=X(IH1)+SMAX*COS(ALPHA)
      SIMIN=AMIN1(SI3I,SI3A)-TOL
      SIMAX=AMAX1(SI3I,SI3A)+TOL
C
C --- LOOK FOR THE ENDPOINT IH3 OF THE LATERAL EDGE
C
      DO IH3=1,NH
      IF(IH3.EQ.IH1.OR.IH3.EQ.IH2)GO TO 80
C --- THE ENDPOINT X-COORDINATE IS NOT IN THE RANGE
      IF(SIMIN.GT.X(IH3).OR.X(IH3).GT.SIMAX)GO TO 80
C --- THE ENDPOINT Y-COORDINATE DOES NOT AGREE
      TA3=Y(IH1)+(X(IH3)-X(IH1))*TAN(ALPHA)
      IF(ABS(TA3-Y(IH3)).GT.TOL)GO TO 80
C --- LATERAL EDGE INADMISSIBLE
      DH13=SQRT((X(IH1)-X(IH3))**2+(Y(IH1)-Y(IH3))**2)
      IF(SMIN.GT.DH13.OR.DH13.GT.SMAX)GO TO 80
C
C --- NOW THREE VERTICES ARE FOUND. CALCULATE THE POSITION OF THE
C --- FOURTH ONE AND LOOK FOR IT
C
      SI4=X(IH3)+X(IH2)-X(IH1)
      DO IH4=1,NH
C --- THE FOURTH VERTEX COINCIDES WITH ANY OF THE THREE
      IF(IH4.EQ.IH1.OR.IH4.EQ.IH2.OR.IH4.EQ.IH3)GO TO 70
C --- THE FOURTH VERTEX X-COORDINATE DOES NOT AGREE
      IF(ABS(SI4-X(IH4)).GT.TOL)GO TO 70
C --- THE FOURTH VERTEX Y-COORDINATE DOES NOT AGREE
      TA4=Y(IH3)+Y(IH2)-Y(IH1)
      IF(ABS(TA4-Y(IH4)).GT.TOL)GO TO 70
C
C --- RECTANGLE FOUND. TEST IT
C
C --- DIFFERENT LENGTHS OF EDGES
      DH34=SQRT((X(IH4)-X(IH3))**2+(Y(IH4)-Y(IH3))**2)
      IF(ABS(DH12-DH34).GT.TOL)GO TO 70
C --- EDGE INADMISSIBLE
      IF(SMIN.GT.DH34.OR.DH34.GT.SMAX)GO TO 70
C --- DIFFERENT LENGTHS OF LATERAL EDGES
      DH24=SQRT((X(IH4)-X(IH2))**2+(Y(IH4)-Y(IH2))**2)
      IF(ABS(DH13-DH24).GT.TOL)GO TO 70
C --- DIFFERENT LENGTHS OF DIAGONALS
      DH14=SQRT((X(IH4)-X(IH1))**2+(Y(IH4)-Y(IH1))**2)
      DH23=SQRT((X(IH2)-X(IH3))**2+(Y(IH2)-Y(IH3))**2)
      IF(ABS(DH14-DH23).GT.TOL)GO TO 70
C
C --- FIND GLOBAL VERTEX NUMBERS OF THE RECTANGLE
C
      IIH1=IH1+NB1
      IF(IH1.GT.NH1)IIH1=IH1+NB2-NH1
      IIH2=IH2+NB1
      IF(IH2.GT.NH1)IIH2=IH2+NB2-NH1
      IIH3=IH3+NB1
      IF(IH3.GT.NH1)IIH3=IH3+NB2-NH1
      IIH4=IH4+NB1
      IF(IH4.GT.NH1)IIH4=IH4+NB2-NH1
C
C --- CHECK FOR DUPLICATES AMONG THE RECTANGLES RECORDED
C
      DO IR=1,NR
      READ(17,REC=IR)IG1,IG2,IG3,IG4
C --- LOOK FOR IIH1 AMONG THE FOUR IG'S
      IF(IIH1.NE.IG1)GO TO 30
C --- IIH1=IG1, TRY TO IDENTIFY THE TWO RECTANGLES
      IF(IIH2.EQ.IG2.AND.IIH4.EQ.IG4.AND.IIH3.EQ.IG3)GO TO 70
      IF(IIH2.EQ.IG3.AND.IIH4.EQ.IG4.AND.IIH3.EQ.IG2)GO TO 70
      GO TO 60
30      IF(IIH1.NE.IG2)GO TO 40
C --- IIH1=IG2, TRY TO IDENTIFY THE TWO RECTANGLES
      IF(IIH2.EQ.IG4.AND.IIH4.EQ.IG3.AND.IIH3.EQ.IG1)GO TO 70
      IF(IIH2.EQ.IG1.AND.IIH4.EQ.IG3.AND.IIH3.EQ.IG4)GO TO 70
      GO TO 60
40      IF(IIH1.NE.IG3)GO TO 50
C --- IIH1=IG3, TRY TO IDENTIFY THE TWO RECTANGLES
      IF(IIH2.EQ.IG1.AND.IIH4.EQ.IG2.AND.IIH3.EQ.IG4)GO TO 70
      IF(IIH2.EQ.IG4.AND.IIH4.EQ.IG2.AND.IIH3.EQ.IG1)GO TO 70
      GO TO 60
50      IF(IIH1.NE.IG4)GO TO 60
C --- IIH1=IG4, TRY TO IDENTIFY THE TWO RECTANGLES
      IF(IIH2.EQ.IG3.AND.IIH4.EQ.IG1.AND.IIH3.EQ.IG2)GO TO 70
      IF(IIH2.EQ.IG2.AND.IIH4.EQ.IG1.AND.IIH3.EQ.IG3)GO TO 70
60      ENDDO
C
C --- THE RECTANGLE DOES NOT COINCIDE WITH OTHERS. RECORD IT
C
      NR=NR+1
      WRITE(17,REC=NR)IIH1,IIH2,IIH3,IIH4

70      ENDDO
C --- THE FOURTH VERTEX NOT FOUND
80      ENDDO
C --- THE THIRD VERTEX NOT FOUND
      IF(SECOND)GO TO 10
90      ENDDO
      ENDDO
      RETURN
      END