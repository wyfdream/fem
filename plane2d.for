*===================================================================*
*     ANALYSIS OF PLATES UNDER IPLANE LOADS                         *
*===================================================================*  
      PROGRAM MAIN
      COMMON A(100000)	          !��������������       
	CALL OPFILE()		          !���ļ��ӳ���     
	CALL INITDAT(NN,NE,ND,NFIX,NP)   !����������Ϣ�ӳ��� 
	CALL ALLOC(NP,NN,NE,ND,NFIX,MLOC,MCX,MCY,MIFIX,MP,MGS,MSTRES,
     &           MSTRAIN,MLP,MRP,MAT,MEND)		         !��ַ(����������)�����ӳ���
      IF(MEND.GT.100000) THEN
        WRITE(*,*)'*************** MEND=',MEND
        STOP 
      ENDIF      
      CALL CST(NN,NE,ND,NFIX,NP,A(MLOC),A(MCX),A(MCY),A(MIFIX),     
     &         A(MP),A(MGS),A(MSTRES),A(MSTRAIN),A(MLP),
     &         A(MRP),A(MAT))!�ӳ���       
      CALL CLFILE()		          !�ر��ļ��ӳ���   
	STOP				          !ֹͣ����       
	END

*===================================================================*
*                    SUBROUTINE OPEN FILE                           *
*===================================================================*       
      SUBROUTINE OPFILE()     !����������ļ�  
      CHARACTER*1024 FILENAME 
     
      
C	OPEN(10,FILE='INP_1.DAT')        
      WRITE(*,*)'please input the file name'      
      READ(*,2040) FILENAME
2040  FORMAT(A)


	OPEN(10,FILE=FILENAME, STATUS='OLD', IOSTAT=IER)      
      IF(IER.NE.0) THEN
        WRITE(*,*)'input file is not correct'
        STOP 
      ENDIF   	
      
	OPEN(30,FILE='OUTP.DAT')       
	
	RETURN       
	END
	
*===================================================================*
*                    SUBROUTINE CLOSE FILE                          *
*===================================================================*   
      SUBROUTINE CLFILE()     !�ر��ļ�д��ͨ��  
	CLOSE(10)       
	CLOSE(30)       
	RETURN       
	END

*===================================================================*
*                  SUBROUTINE READ INITIAL DATA                      *
*===================================================================*       
      SUBROUTINE INITDAT(NN,NE,ND,NFIX,NP)  !��ȡ�����ļ�����
*      NN:�����   NE:��Ԫ��  ND:���ɶ���  NFIX:Լ����  NP �غ���
      READ(10,*) NN,NE,NFIX,NP
      ND=NN*2       
	RETURN       
	END
	
	
*===================================================================*
*                    SUBROUTINE ALLOCATE MEMORY                     *
*===================================================================*       
      SUBROUTINE ALLOC(NP,NN,NE,ND,NFIX,MLOC,MCX,MCY,MIFIX,MP,MGS,	
     &	             MSTRES,MSTRAIN,MLP,MRP,MAT,MEND)       
      MLOC   =1       
	MCX    =MLOC +NE*3       
	MCY    =MCX  +NN       
	MIFIX  =MCY  +NN       
	MP     =MIFIX+NFIX       
	MGS    =MP   +ND       
	MSTRES =MGS  +ND*ND       
	MSTRAIN=MSTRES+NE*3
      MLP    =MSTRAIN+NE*3
      MRP    =MLP+NP
      MAT    =MRP+NP*2       
      MEND   =MAT+ND       
	RETURN       
	END


*===================================================================*
*                    SUBROUTINE CST                                 *
*===================================================================*       
      SUBROUTINE CST(NN,NE,ND,NFIX,NP,LOC,CX,CY,IFIX,P,GS,STRES,STRAIN,
     &               LP,RP,AT)	
      DIMENSION LOC(NE,3),CX(NN),CY(NN),IFIX(NFIX),P(ND),GS(ND,ND),   
     &          STRES(NE,3),STRAIN(NE,3),AA(6,6),N(6),CXL(3),CYL(3),    
     &	      QG(6),QL(6),DD(3,3),BB(3,6),DB(3,6),SLOC(3),AT(ND),
     &          LP(NP),RP(NP,2)       
	CALL READDAT(LOC,CX,CY,IFIX,LP,RP,NP,NE,NN,NFIX,ND,E,ANU,T)  !���������ӳ���       
      CALL INITIALDAT(GS,P,LP,RP,NP,ND)                    !��ʼ���ӳ���
      CALL STIF(NN,NE,ND,LOC,CX,CY,GS,E,T,ANU)             !�γɵ�Ԫ�նȾ���	
	CALL PROCEGS(GS,IFIX,P,NFIX,ND)                      !Լ�������ӳ���       
	CALL DECOMP(ND,GS,AT)			                     !��������ӳ���       
	CALL SOLVE(ND,GS,P)      
      CALL STRESS_STRAIN(NN,NE,ND,LOC,CX,CY,P,GS,STRES,STRAIN,E,T,ANU)   !�ش���ԪӦ���Ӧ��	
	CALL OUTPUT(P,STRES,STRAIN,NE,ND)	                 !�������ӳ���       
	RETURN								                 !���ӳ��򷵻�       
	END
C
      SUBROUTINE STIF(NN,NE,ND,LOC,CX,CY,GS,E,T,ANU)	
      DIMENSION LOC(NE,3),CX(NN),CY(NN),GS(ND,ND)
      DIMENSION CXL(3),CYL(3),AKL(6,6),AKG(6,6),ALAMDA(6,6),
     &          ALAMDAT(6,6),AA(6,6),N(6) 
	DO 100 I=1,NE       
	J =LOC(I,1)       
	K =LOC(I,2)       
	L =LOC(I,3)
*      �γ�ת������ LAMDA �ӳ���       
      CALL LAMBDA(J,K,L,NN,CX,CY,C1,C2,CXL,CYL,E,T,ANU,ALAMDA,AREA)       
	CALL CALCKL(AKL,C1,C2,CXL,CYL,ANU) !�γɾֲ�����ϵ�ĵ�Ԫ�ն� KL(��)       
	CALL SYMTRA(AKL,ALAMDA,ALAMDAT)	   !�γ�ת�� LAMDAT ������ KL       
	CALL MATMUL_EL(AKL,ALAMDA,AA,6,6,6)       
	CALL MATMUL_EL(ALAMDAT,AA,AKG,6,6,6)       
	CALL ASSEMB(J,K,L,N,AKG,GS,ND)	   !������װ�ܸ��ӳ���
100   CONTINUE       
      RETURN
	END
C
      SUBROUTINE STRESS_STRAIN(NN,NE,ND,LOC,CX,CY,P,GS,STRES,
     &                                             STRAIN,E,T,ANU)	
      DIMENSION LOC(NE,3),CX(NN),CY(NN),P(ND),STRES(NE,3),STRAIN(NE,3)
      DIMENSION CXL(3),CYL(3),AKL(6,6),AKG(6,6),ALAMDA(6,6),
     &          ALAMDAT(6,6),AA(6,6),N(6),QG(6),QL(6),BB(3,6),
     &          DD(3,3),DB(3,6),SLOC(3)
	DO 200 I=1,NE       
	J=LOC(I,1)     
	K=LOC(I,2)       
	L=LOC(I,3)       
	CALL LAMBDA(J,K,L,NN,CX,CY,C1,C2,CXL,CYL,E,T,ANU,ALAMDA,AREA)	   
	CALL TRANSF(J,K,L,ND,N,QG,QL,P,ALAMDA) !�γɾֲ�����ϵ��λ������ QL       
	CALL MATRIXB(CXL,CYL,AREA,BB)		   !�γ� B ����(BB)        
	CALL MATRIXD(DD,ANU,E)				   !�γɵ��Ծ��� DD       
	CALL MATMUL_EL(DD,BB,DB,3,3,6)
*      ��ԪӦ����Ӧ���ӳ����ӳ���	   
	CALL BKBRING(I,NE,SLOC,DB,BB,QL,ALAMDA,STRES,STRAIN)
200  	CONTINUE       
      RETURN
	END
      
*===================================================================*
*                SUBROUTINE INITIALIZE DATA VARIABLE                *
*===================================================================*
      SUBROUTINE INITIALDAT(GS,P,LP,RP,NP,ND)                    !��ʼ���ӳ���
      DIMENSION GS(ND,ND),P(ND),LP(NP),RP(NP,2)
	DO I=1,ND
         P(I)=0.0
	   DO J=1,ND
	     GS(I,J)=0.0
	   ENDDO
	ENDDO
      DO I=1,NP
        IX=LP(I)*2-1
        IY=LP(I)*2
        P(IX)=P(IX)+RP(I,1)
        P(IY)=P(IY)+RP(I,2)
      ENDDO
      RETURN
	END

*===================================================================*
*                SUBROUTINE LAMBDA                                  *
*===================================================================*      
      SUBROUTINE LAMBDA(J,K,L,NN,CX,CY,C1,C2,CXL,CYL,E,T,ANU,    
     &                       LAMDA,AREA)      
      REAL LAMDA(6,6)      
	DIMENSION CX(NN),CY(NN),CXL(3),CYL(3)      
	D12=SQRT((CX(K)-CX(J))**2+(CY(K)-CY(J))**2)      
	DCL12=(CX(K)-CX(J))/D12      
	DCM12=(CY(K)-CY(J))/D12      
	D14=DCL12*(CX(L)-CX(J))+DCM12*(CY(L)-CY(J))      
	D43=SQRT((CX(L)-CX(J))**2+(CY(L)-CY(J))**2-D14**2)      
	CXL(1)=0.0      
	CYL(1)=0.0      
	CXL(2)=0.0      
	CYL(2)=D12      
	CXL(3)=D43      
	CYL(3)=D14      
	DCL43=(CX(L)-CX(J)-DCL12*D14)/D43      
	DCM43=(CY(L)-CY(J)-DCM12*D14)/D43      
	AREA=((CXL(3)-CXL(2))*(CYL(2)-CYL(1))-    
     &      (CXL(2)-CXL(1))*(CYL(3)-CYL(2)))/2.0      
      C1=E*T/(4.0*AREA*(1.0-ANU**2))      
	C2=E*T/(8.0*AREA*(1.0+ANU)) 	
	DO II=1,6
	  DO JJ=1,6
	   LAMDA(II,JJ)=0.0
	  ENDDO
	ENDDO   
	LAMDA(1,1)=DCL43      
	LAMDA(1,2)=DCM43      
	LAMDA(2,1)=DCL12      
	LAMDA(2,2)=DCM12      
	LAMDA(3,3)=DCL43      
	LAMDA(3,4)=DCM43      
	LAMDA(4,3)=DCL12
      LAMDA(4,4)=DCM12      
	LAMDA(5,5)=DCL43      
	LAMDA(5,6)=DCM43      
	LAMDA(6,5)=DCL12      
	LAMDA(6,6)=DCM12      
	RETURN      
	END
	
*===================================================================*
*              SUBROUTINE READ DATA                                 *
*===================================================================*      
      SUBROUTINE READDAT(LOC,CX,CY,IFIX,LP,RP,NP,NE,NN,NFIX,ND,E,ANU,T)      
	DIMENSION LOC(NE,3),CX(NN),CY(NN),IFIX(NFIX),LP(NP),RP(NP,2)      
*      E :����ģ��  ANU:���ɱ�  T:���
	READ(10,*) E,ANU,T       
	READ(10,*)(IDX,(LOC(I,J),J=1,3),I=1,NE)      
	READ(10,*)(IDX,CX(I),CY(I),I=1,NN)      
	READ(10,*)(IDX,IFIX(J),J=1,NFIX)      
	READ(10,*)(IDX,LP(I),RP(I,1),RP(I,2),I=1,NP)      
	RETURN      
	END
*===================================================================*
*              SUBROUTINE CALCULATE KL                              *
*===================================================================*
      SUBROUTINE CALCKL(KL,C1,C2,CXL,CYL,ANU)      
	DIMENSION CXL(3),CYL(3)      
	REAL KL(6,6)      
	Y32 =CYL(3) -CYL(2)      
	X32 =CXL(3) -CXL(2)      
	Y31 =CYL(3) -CYL(1)      
	X31 =CXL(3) -CXL(1)      
	Y21 =CYL(2) -CYL(1)      
	X21 =CXL(2) -CXL(1)      
	KL(1,1) =C1*(Y32**2) +C2*(X32**2)      
	KL(2,1) =-C1*ANU*Y32*X32 -C2*X32*Y32      
	KL(2,2) =C1 *(X32**2) +C2 *(Y32**2)      
	KL(3,1) =-C1*Y32*Y31 -C2*X32*X31      
	KL(3,2) =C1*ANU*X32*Y31 +C2*Y32*X31      
	KL(3,3) =C1*(Y31**2) +C2*(X31**2)      
	KL(4,1) =C1*ANU*Y32*X31 +C2*X32*Y31      
	KL(4,2) =-C1*X32*X31 -C2*Y32*Y31      
	KL(4,3) =-C1*ANU*Y31*X31 -C2*X31*Y31      
	KL(4,4) =C1*(X31**2) +C2*(Y31**2)      
	KL(5,1) =C1*Y32*Y21 +C2*X32*X21      
	KL(5,2)=-C1*ANU*X32*Y21-C2*X21*Y32      
	KL(5,3)=-C1*Y31*Y21-C2*X31*X21      
	KL(5,4)=C1*ANU*X31*Y21+C2*Y31*X21      
	KL(5,5)=C1*(Y21**2)+C2*(X21**2)      
	KL(6,1)=-C1*ANU*Y32*X21-C2*X32*Y21      
	KL(6,2)=C1*X32*X21+C2*Y32*Y21      
	KL(6,3)=C1*ANU*Y31*X21+C2*X31*Y21      
	KL(6,4)=-C1*X31*X21-C2*Y31*Y21      
	KL(6,5)=-C1*ANU*Y21*X21-C2*X21*Y21      
	KL(6,6)=C1*(X21**2)+C2*(Y21**2)      
	RETURN      
	END


*===================================================================*
*              SUBROUTINE SYMETRY TRANSFORM                         *
*===================================================================*      
      SUBROUTINE SYMTRA(KL,LAMDA,LAMDAT)      
	REAL KL(6,6),LAMDA(6,6),LAMDAT(6,6)
      DO II=1,6
	  DO JJ=1,6
	  KL(II,JJ)=KL(JJ,II)
	  ENDDO
	ENDDO 
      DO II=1,6
	  DO JJ=1,6
	    LAMDAT(II,JJ)=LAMDA(JJ,II)
	  ENDDO
	ENDDO    
      
	RETURN      
	END
	
*===================================================================*
*              SUBROUTINE ASSEMBLE GS                               *
*===================================================================*     
      SUBROUTINE ASSEMB(J,K,L,N,KG,GS,ND)      
	DIMENSION N(6),GS(ND,ND)      
	REAL KG(6,6)      
	N(1)=J*2-1      
	N(2)=J*2      
	N(3)=K*2-1      
	N(4)=K*2      
	N(5)=L*2-1      
	N(6)=L*2      
	DO 40 II=1,6     
	DO 40 JJ=1,6      
	IK=N(II)      
	JK=N(JJ)      
	GS(IK,JK)=GS(IK,JK)+KG(II,JJ)
40    CONTINUE      
      RETURN      
	END


*===================================================================*
*              SUBROUTINE PROCESS GS                                *
*===================================================================*
      SUBROUTINE PROCEGS(GS,IFIX,P,NFIX,ND)      
	DIMENSION GS(ND,ND),IFIX(NFIX),P(ND)      
	DO 50 I=1,NFIX      
	IX=IFIX(I)      
	P(IX) =0.0 
	DO J=1,ND
	  GS(J,IX)=0.0
	  GS(IX,J)=0.0
	ENDDO       
	GS(IX,IX)=1.0
50    CONTINUE      
	RETURN      
	END

*===================================================================*
*              SUBROUTINE TRANSFORM                                 *
*===================================================================*
      SUBROUTINE TRANSF(J,K,L,ND,N,QG,QL,P,LAMDA)       
	DIMENSION N(6),QG(6),QL(6),P(ND)       
	REAL LAMDA(6,6)       
	N(1)=J*2-1       
	N(2)=J*2       
	N(3)=K*2-1       
	N(4)=K*2       
	N(5)=L*2-1       
	N(6)=L*2       
      DO II=1,6
	 NN=N(II)
	 QG(II)=P(NN)
	ENDDO     
	DO 70 II=1,6       
	QL(II)=0.0	      
	DO JJ=1,6
	  QL(II)=QL(II)+LAMDA(II,JJ)*QG(JJ) 
	ENDDO
70    CONTINUE       
	RETURN       
	END

*===================================================================*
*              SUBROUTINE MATRIXB                                   *
*===================================================================*
      SUBROUTINE MATRIXB(CXL,CYL,AREA,BB)       
	DIMENSION CXL(3),CYL(3),BB(3,6)       
	Y32=(CYL(3)-CYL(2))/(2.0*AREA)       
	X32=(CXL(3)-CXL(2))/(2.0*AREA)       
	Y31=(CYL(3)-CYL(1))/(2.0*AREA)       
	X31=(CXL(3)-CXL(1))/(2.0*AREA)       
	Y21=(CYL(2)-CYL(1))/(2.0*AREA)       
	X21=(CXL(2)-CXL(1))/(2.0*AREA)       
	DO 80 II=1,3       
	DO 80 JJ=1,6
80    BB(II,JJ)=0.0       
	BB(1,1)=Y32       
	BB(1,3)=-Y31       
	BB(1,5)=Y21       
	BB(2,2)=-X32       
	BB(2,4)=X31      
	BB(2,6)=-X21       
	BB(3,1)=-X32       
	BB(3,2)=Y32       
	BB(3,3)=X31       
	BB(3,4)=-Y31       
	BB(3,5)=-X21       
	BB(3,6)=Y21       
	RETURN       
	END

*===================================================================*
*              SUBROUTINE MATRIXD                                   *
*===================================================================*
      SUBROUTINE MATRIXD(DD,ANU,E)       
	DIMENSION DD(3,3)       
	ENU=E/(1.0-ANU**2)       
	DD(1,1)=ENU       
	DD(1,2)=ENU*ANU       
	DD(1,3)=0.0       
	DD(2,1)=ENU*ANU       
	DD(2,2)=ENU       
	DD(2,3)=0.0       
	DD(3,1)=0.0       
	DD(3,2)=0.0       
	DD(3,3)=ENU*(1.0-ANU)/2.0       
	RETURN       
	END
	
*===================================================================*
*              SUBROUTINE BACK BRING                                *
*===================================================================* 
      SUBROUTINE BKBRING(I,NE,SLOC,DB,BB,QL,LAMDA,STRES,STRAIN)       
	DIMENSION SLOC(3),DB(3,6),BB(3,6),QL(6),STRES(NE,3),STRAIN(NE,3)       
	REAL LAMDA(6,6)       
	DO 90 II=1,3       
	SLOC(II)=0.0       
	DO 90 JJ=1,6       
	SLOC(II)=SLOC(II)+DB(II,JJ)*QL(JJ)
90    CONTINUE       
      AL1=LAMDA(1,1)       
	AM1=LAMDA(1,2)       
	AL2=LAMDA(2,1)       
	AM2=LAMDA(2,2)
      STRES(I,1)=SLOC(1)*AL1**2+SLOC(2)*AL2**2+2.0*SLOC(3)*AL1*AL2       
	STRES(I,2)=SLOC(1)*AM1**2+SLOC(2)*AM2**2+2.0*SLOC(3)*AM1*AM2       
	STRES(I,3)=SLOC(1)*AL1*AM1+SLOC(2)*AL2*AM2+SLOC(3)*     
     &            (AL1*AM2+AL2*AM1)       
      DO 91 II=1,3       
	SLOC(II)=0.0       
	DO 91 JJ=1,6       
	SLOC(II)=SLOC(II)+BB(II,JJ)*QL(JJ)
91    CONTINUE       
	AL1=LAMDA(1,1)       
	AM1=LAMDA(1,2)       
	AL2=LAMDA(2,1)       
	AM2=LAMDA(2,2)       
	STRAIN(I,1)=SLOC(1)*AL1**2+SLOC(2)*AL2**2+2.0*SLOC(3)*AL1*AL2       
	STRAIN(I,2)=SLOC(1)*AM1**2+SLOC(2)*AM2**2+2.0*SLOC(3)*AM1*AM2       
	STRAIN(I,3)=SLOC(1)*AL1*AM1+SLOC(2)*AL2*AM2+SLOC(3)*     
     &            (AL1*AM2+AL2*AM1)       
      RETURN       
	END

*===================================================================*
*              SUBROUTINE OUTPUT                                    *
*===================================================================* 
      SUBROUTINE OUTPUT(P,STRES,STRAIN,NE,ND)      
	DIMENSION P(2,ND/2),STRES(NE,3),STRAIN(NE,3)      
	WRITE(30,*)'================DISPLACEMENT OF NODES=============='      
	WRITE(30,'(/,I,4X,2F17.5)')(I,(P(J,I),J=1,2),I=1,ND/2)
	WRITE(30,*)' '
	WRITE(30,*)'=================STRAINS IN ELEMENTS==============='      
	WRITE(30,'(/,I,4X,3F17.5)')(I,(STRAIN(I,J),J=1,3),I=1,NE)   
	WRITE(30,*)' '	   
	WRITE(30,*)'=================STRESSES IN ELEMENTS=============='      
	WRITE(30,'(/,I,4X,3F17.5)')(I,(STRES(I,J),J=1,3),I=1,NE)      
	RETURN      
	END
	
      SUBROUTINE OUTPUT_2(P,STRES,STRAIN,NE,ND)      
	DIMENSION P(ND),STRES(NE,3),STRAIN(NE,3)      
	WRITE(30,*)'================DISPLACEMENT OF NODES=============='      
	WRITE(30,'(/,9X,2E17.5)')P      
	WRITE(30,*)'=================STRESSES IN ELEMENTS=============='      
	WRITE(30,'(/,2X,3E17.5)')((STRES(I,J),J=1,3),I=1,NE)      
	WRITE(30,*)'=================STRAINS IN ELEMENTS==============='      
	WRITE(30,'(/,2X,3E17.5)')((STRAIN(I,J),J=1,3),I=1,NE)      
	RETURN      
	END	

*===================================================================*
*              SUBROUTINE DECOMP                                    *
*===================================================================* 
      SUBROUTINE DECOMP(ND,GS,AT)       
	DIMENSION GS(ND,ND),AT(ND)       
	DO 10 I=1,ND       
	DO 10 J=1,I       
	S=GS(I,J)       
	I1=J-1 
	DO K=1,I1
	  S=S-GS(J,K)*AT(K)   
	ENDDO     
	AT(J)=S       
	IF(I.EQ.J) GS(I,J)=S       
	IF(I.NE.J) GS(I,J)=S/GS(J,J)
10    CONTINUE       
      RETURN       
	END

*===================================================================*
*              SUBROUTINE SOLVE                                     *
*===================================================================* 
      SUBROUTINE SOLVE(ND,GS,P)      
	DIMENSION GS(ND,ND),P(ND)
	DO I=1,ND
	  DO J=1,I-1
       	P(I)=P(I)-GS(I,J)*P(J)	 
	  ENDDO
	ENDDO      
      
	DO I=1,ND
	 P(I)=P(I)/GS(I,I)
	ENDDO  
	
	DO I=1,ND-1
	  DO J=ND-I+1,ND
          P(ND-I)=P(ND-I)-GS(J,ND-I)*P(J)
	  ENDDO
	ENDDO 
	      
      RETURN      
	END

*===================================================================*
*               SUBROUTINE MATMUL_EL                                   *
*===================================================================*      
      SUBROUTINE MATMUL_EL(A,B,C,N1,N2,N3)	  !����˷��ӳ���      
	DIMENSION A(N1,N2), B(N2,N3),C(N1,N3)      
	DO I=1,N1
	  DO J=1,N3     
	   S=0.0
	   DO K=1,N2
	     S=S+A(I,K)*B(K,J)
	   ENDDO 
	   C(I,J)=S
	  ENDDO
	ENDDO     
      RETURN      
	END
