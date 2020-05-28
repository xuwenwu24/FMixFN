IMPLICIT NONE
REAL,ALLOCATABLE::MARKER(:,:),& ! marker matrix of training set
					Y(:,:),&    ! vector of phenotype, the second column is sex in training set
					PHENO(:),&	! used to calculated the variance of phenotypes
					YY(:,:),&	! vector of phenotype, the second column is sex in validation set
					UUU(:),&     ! vector of fixed effects
					G(:),&      ! vector of SNP effect
					P(:),&      ! used to calculate the allele frequency to standardize
					T_MATRIX(:,:),T_MM(:,:)  ! transformed (standardized) marker matrix
REAL,ALLOCATABLE::TBV(:),&    ! true breeding value
				    EBV(:),&    ! estimated breeding value
				    MM(:,:)       ! marker matrix of validation population
CHARACTER*20,ALLOCATABLE::INDI_TRAIN(:),INDI_VAL(:),&   ! the name of individuals in training and validation population
						  TRAIN(:),&    ! read the individuals in training_phenotype.txt, maybe not the same with genotype.txt
						  VAL(:),&		! read the individuals in validation_phenotype.txt, maybe not the same with genotype.txt
						  SNP(:)		! the name of SNP
INTEGER::NUMBER,RANDOM,N_FIX,&
         I,J,K,II,JJ,KK,&
		 NVAL,INDEX,INITIAL,CIRCLE
REAL::A,HERITABILITY,MEAN,SD2,XX,T,&
	  VarE,&	! variance of environmental
	  FIX_EFFECT,EBV_EFFECT,&
	  PI,&          ! the ratio of larger variance in genome
	  V1,V2,V3,V4,VG,VE,&    ! variances
	  COR,REG	! used to calculate accuracy
CHARACTER*10::ID,TT
character*100::date,time,zone
integer::begin(8),over(8)

OPEN(15,FILE='parameter.txt')
READ(15,*)NUMBER    ! number of training set
READ(15,*)RANDOM    ! number of SNP
READ(15,*)N_FIX     ! number of fixed effects
READ(15,*)PI        ! the ratio 
READ(15,*)HERITABILITY     ! heritability
READ(15,*)NVAL      ! number of validation population
READ(15,*)INDEX     ! if index is 1, predictive accuracy, else (0), doesn't predictive accuracy
CLOSE(15)

ALLOCATE(MARKER(NUMBER,RANDOM),Y(NUMBER,N_FIX+1),PHENO(NUMBER),G(RANDOM),P(RANDOM),T_MATRIX(NUMBER,RANDOM))
ALLOCATE(TBV(NVAL),EBV(NVAL),MM(NVAL,RANDOM),INDI_TRAIN(NUMBER),TRAIN(NUMBER),YY(NVAL,N_FIX+1),&
        T_MM(NVAL,RANDOM),INDI_VAL(NVAL),VAL(NVAL),SNP(RANDOM),UUU(N_FIX))

OPEN(11,FILE='training_genotype.txt')
OPEN(12,FILE='training_phenotype.txt')
OPEN(13,FILE='validation_genotype.txt')
READ(11,*)ID,SNP(:)    ! read snp name
DO I=1,NUMBER
	READ(11,*)INDI_TRAIN(I),MARKER(I,:)	! read name, marker genotype of training set
ENDDO

DO I=1,NUMBER
	READ(12,*)TRAIN(I),Y(I,:)	! read name, phenotype and sex of training set
	PHENO(I)=Y(I,1)
ENDDO
CALL VARIANCE(PHENO,NUMBER,MEAN,SD2)
VE=SD2*(1-HERITABILITY)
VG=SD2*HERITABILITY
DO I=1,NUMBER				! to make the ID consistent in INDI_TRAIN and TRAIN
	DO J=1,NUMBER
		IF(TRAIN(J)==INDI_TRAIN(I))THEN
			TT=TRAIN(I);TRAIN(I)=TRAIN(J);TRAIN(J)=TT
			DO K=1,N_FIX+1
				T=Y(I,K);Y(I,K)=Y(J,K);Y(J,K)=T
			ENDDO
		ENDIF
	ENDDO
ENDDO

READ(13,*)ID,SNP(:)					! read SNP name
DO I=1,NVAL
	READ(13,*)INDI_VAL(I),MM(I,:)				! read name, marker genotype of validation set
ENDDO
CLOSE(11);CLOSE(12);CLOSE(13)
!...................   .....Now, standardize the genotypic records
p=SUM(MARKER,1)/(2.0*NUMBER)                         	! calculate frequency of reference allele (=1) for each SNP in training set
T_MATRIX=0.0
DO J=1,RANDOM                                    ! standardize B(i,j) [ = (b-mean)/SD ]
	IF(P(J)==0.0 .OR. P(J)==1.0)THEN
		DO I=1,NUMBER
			T_MATRIX(I,J)=0.0
		ENDDO
	ELSE
		DO I=1,NUMBER
			T_MATRIX(I,J)=(MARKER(I,J)-2.0*p(J))/SQRT(2.0*P(J)*(1-P(J)))
		ENDDO
	ENDIF
ENDDO

p=SUM(MM,1)/(2.0*NVAL)                         	! calculate frequency of reference allele (=1) for each SNP in validation set
T_MM=0.0
DO J=1,RANDOM
	IF(P(J)==0.0 .OR. P(J)==1.0)THEN
		DO I=1,NVAL
			T_MM(I,J)=0.0
		ENDDO
	ELSE                                    ! standardize B(i,j) [ = (b-mean)/SD ]
		DO I=1,NVAL
			T_MM(I,J)=(MM(I,J)-2.0*p(J))/SQRT(2.0*P(J)*(1-P(J)))
		ENDDO
	ENDIF
ENDDO
!............................................
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
call date_and_time(date,time,zone,begin)
CALL MixP(NUMBER,RANDOM,N_FIX,Y,T_MATRIX,PI,HERITABILITY,VG,VE,G,UUU)	 ! Calculate the Iteration Conditional Expectation algorithm
call date_and_time(date,time,zone,over)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!............................................now, save the effects of SNP
OPEN(14,FILE='Effect of Marker.txt')
WRITE(14,*)'No.marker      Effect'
DO I=1,RANDOM
	WRITE(14,'(I12,f12.6,f12.6)')I,G(I)
ENDDO
CLOSE(14)
!............................................now, doing correlation analysis (predictive accuracy)
DO I=1,NVAL
	EBV_EFFECT=0.0
	DO J=1,RANDOM
		EBV_EFFECT=EBV_EFFECT+T_MM(I,J)*G(J)
	ENDDO
	EBV(I)=EBV_EFFECT
ENDDO

OPEN(17,FILE='ebv.txt')
WRITE(17,*)'Individual               EBV'
DO I=1,NVAL
	WRITE(17,'(A20,f10.4)')INDI_VAL(I),EBV(I)
ENDDO
CLOSE(17)

IF(INDEX==1)THEN     ! if need to calculate accuracy
    OPEN(16,FILE='validation_phenotype.txt')
    DO I=1,NVAL
        READ(16,*)VAL(I),YY(I,:)            ! reading TBV (or phenotype) from files
    ENDDO
    CLOSE(16)

    DO I=1,NVAL							! to make the ID consistent in INDI_VAL and VAL
        DO J=1,NVAL
            IF(VAL(J)==INDI_VAL(I))THEN
                TT=VAL(I);VAL(I)=VAL(J);VAL(J)=TT
                DO K=1,N_FIX+1
                    T=YY(I,K);YY(I,K)=YY(J,K);YY(J,K)=T
                ENDDO
            ENDIF
        ENDDO
    ENDDO

    DO I=1,NVAL
        FIX_EFFECT=0.0
        DO J=1,N_FIX
            FIX_EFFECT=FIX_EFFECT+YY(I,J+1)*UUU(J)
        END DO
        TBV(I)=YY(I,1)-FIX_EFFECT
    ENDDO
	CALL CORRELATE(NVAL,EBV,TBV,COR,REG)
	PRINT*,'The accuracy of MixP in this population is:',COR,REG
	OPEN(18,FILE='Predictive accuracy.txt')
	WRITE(18,*)'Training size (number of individuals in training set) is:',NUMBER
	WRITE(18,*)'Number of marker loci is:',RANDOM
	WRITE(18,*)'The accuracy of MixP in this population is:',COR,REG
	CLOSE(18)
ENDIF

write(*,101)begin(1:3),begin(5:7)
write(*,102)over(1:3),over(5:7)
101 format('the start time is: date = ',i4,'/',i2,'/',i2,', time = ',i2,':',i2,':',i2)
102 format('the end time is: date = ',i4,'/',i2,'/',i2,', time = ',i2,':',i2,':',i2)

END ! END all program

subroutine variance(a,n1,mean,cov)     ! calculate variance
integer::n1
real::b,c,cov,mean
real::a(n1)
mean=sum(a(:))/n1
b=sum(a(:)**2)
c=(sum(a(:)))**2/real(n1)
cov=(b-c)/real(n1-1)
end subroutine variance

SUBROUTINE Correlate(n,VecX,VecY,rXY,bYonX)		! correlate X and Y + get regression coefficient of Y on X
implicit none
integer::i
integer::n
real:: meanX,meanY,SDx,SDy,sumXY
real:: rXY, bYonX
real:: VecX(n), VecY(n)
sumXY=0
do i=1,n
	sumXY=sumXY+VecX(i)*VecY(i)
end do
call MeanSD(n,VecX,meanX,SDx)		        	! calculate mean & SD of X
call MeanSD(n,VecY,meanY,SDy)		        	! calculate mean & SD of Y
rXY  =(sumXY-n*meanX*meanY)/((n-1)*SDx*SDy)
bYonX=(sumXY-n*meanX*meanY)/((n-1)*SDx**2)
END SUBROUTINE Correlate

SUBROUTINE MeanSD(n,VecX,meanX,SDx)    		! calculate mean & SD of vector VecX
implicit none
integer::i,j
integer::n
real::meanX,SDx
real::sum1,sum2
real::VecX(n)       ! vector used to calculate mean & SD
sum1=0
sum2=0
do j=1,n
	sum1 = sum1 + VecX(j)
	sum2 = sum2 + VecX(j)**2
end do
meanX=sum1/n                               		! meanX = mean of the n elements of VecX
SDx=SQRT((sum2-sum1**2/n)/(n-1))           		!   SDx = SD   of the n elements of VecX
END SUBROUTINE MeanSD

SUBROUTINE MixP(N,SNP,N_FIX,PHENO,DESIGN,P,HERITABILITY,VG,VE,G,U)
IMPLICIT NONE
INTEGER::N,&	    ! number of training set
		 SNP,&	    ! number of SNP
		 N_FIX,&    ! number of fixed effects
		 K			! count the cycle
INTEGER::I,J,II
REAL::G(SNP),&			    ! the vector of SNP effects
		PHENO(N,N_FIX+1),&	! the first column is phenotype, the second column is sex
		DESIGN(N,SNP),&     ! the matrix of marker genotype
		FIX(N),&            ! sum of fixed effects
		XG(N),&             ! sum of random effects
		VG,&                ! total additive genetic variance
        HERITABILITY,&      !
		V1,V2,V3,V4,&       !  variances
		VE,&			    ! the variance of errors
		U(N_FIX),&			! overall mean
		Yi(N),&			    ! y-XU-bg(exclude the SNP needed to be calculated)
		GLAST(SNP),&		! used to judge the convergence
		ULAST(N_FIX),&				! used to judge the convergence
		COE_CONVER,COE_LAST,&		! convergent_coefficient
		COE_NUMERATOR,COE_DENOMINATOR
REAL::L,LL,BB,Y,SIGMA2,	&	! VE/(b'b)
		GAMMA,&				! the probability of marker having effects
		LAMDA				! parameter of exponential distribution
REAL::NUMERATOR,&		    ! numerator of E(g|Y)
		DENOMINATOR,&	    ! denominator of E(g|Y)
		P,&                 ! the ratio of SNP with larger variance
		PI,&
		SOR
if (HERITABILITY>=0.2 .and. HERITABILITY<=0.4) then
 V1=0.8367*VG/(P*SNP)
 V2=0.1246*VG/(P*SNP)
 V3=0.0342*VG/(P*SNP)
 V4=0.0043*VG/(P*SNP)
else if (HERITABILITY>0.4) then
 V1=0.8752*VG/(P*SNP)
 V2=0.0958*VG/(P*SNP)
 V3=0.0256*VG/(P*SNP)
 V4=0.0032*VG/(P*SNP)
else if (HERITABILITY<0.2) then
 V1=0.8225*VG/(P*SNP)
 V2=0.1413*VG/(P*SNP)
 V3=0.0324*VG/(P*SNP)
 V4=0.0036*VG/(P*SNP)
end if
PI=3.14159265358
SOR=1.0
G=0.0         ! the initial value of marker
U=0.0        ! the initial value of over_all
K=0           ! count the cycle
OPEN(19,FILE='Process of iteration.txt')
WRITE(19,*)'Iteration	convergent_coefficient'
COE_CONVER=10.0
!**************************************************************************************
XG=0.0
DO I=1,N
    DO J=1,SNP
        XG(I)=XG(I)+DESIGN(I,J)*G(J)
    ENDDO
ENDDO
FIX=0.0
DO I=1,N
    DO J=1,N_FIX
        FIX(I)=FIX(I)+PHENO(I,J+1)*U(J)
    END DO
END DO
!**************************************************************************************
DO WHILE(COE_CONVER > 1E-8)
	GLAST=G
	ULAST=U
    DO I=1,N_FIX
        DO J=1,N
            FIX(J)=FIX(J)-PHENO(J,I+1)*U(I)
        END DO
        L=0.0
        DO J=1,N
            L=L+PHENO(J,I+1)*(PHENO(J,1)-XG(J)-FIX(J))     ! calculate x'(y-xg)
        ENDDO
        LL=0.0
        DO J=1,N
            LL=LL+PHENO(J,I+1)*PHENO(J,I+1)
        ENDDO
        U(I)=L/LL						! calculate overall mean
        U(I)=(1.0-SOR)*ULAST(I)+SOR*U(I)
        DO J=1,N
            FIX(J)=FIX(J)+PHENO(J,I+1)*U(I)
        END DO
    ENDDO
!**************************************************************************************
	DO I=1,SNP                       ! sample the effect of G(I)
		DO J=1,N
			XG(J)=XG(J)-DESIGN(J,I)*G(I)    ! require G(I)=0.0
		ENDDO
		DO J=1,N
			Yi(J)=PHENO(J,1)-XG(J)-FIX(J)
		ENDDO
		BB=0.0
		DO J=1,N
			BB=BB+DESIGN(J,I)*DESIGN(J,I)
		ENDDO
		Y=0.0
		DO J=1,N
			Y=Y+DESIGN(J,I)*Yi(J)
		ENDDO
		Y=Y/BB
		SIGMA2=VE/BB
!*************************************************************Now, calculate E(g|Y)
        NUMERATOR=Y*V1/(V1+SIGMA2)+EXP(0.5*(Y**2/(V1+SIGMA2)-&
        Y**2/(V2+SIGMA2)))*SQRT(V1+SIGMA2)/SQRT(V2+SIGMA2)*Y*V2/(V2+SIGMA2)+&
        EXP(0.5*(Y**2/(V1+SIGMA2)-Y**2/(V3+SIGMA2)))*SQRT(V1+SIGMA2)/SQRT(V3+SIGMA2)*Y*V3/(V3+SIGMA2)+&
        EXP(0.5*(Y**2/(V1+SIGMA2)-Y**2/(V3+SIGMA2)))*SQRT(V1+SIGMA2)/SQRT(V4+SIGMA2)*Y*V4/(V4+SIGMA2)
        DENOMINATOR=1+EXP(0.5*(Y**2/(V1+SIGMA2)-Y**2/(V2+SIGMA2)))*SQRT(V1+SIGMA2)/SQRT(V2+SIGMA2)+&
        EXP(0.5*(Y**2/(V1+SIGMA2)-Y**2/(V3+SIGMA2)))*SQRT(V1+SIGMA2)/SQRT(V3+SIGMA2)+&
        EXP(0.5*(Y**2/(V1+SIGMA2)-Y**2/(V4+SIGMA2)))*SQRT(V1+SIGMA2)/SQRT(V4+SIGMA2)
		G(I)=NUMERATOR/DENOMINATOR
!*************************************************************
        G(I)=(1.0-SOR)*GLAST(I)+SOR*G(I)
		DO J=1,N
			XG(J)=XG(J)+DESIGN(J,I)*G(I)    ! require G(I)=0.0
		ENDDO
	ENDDO
!****************************************************************************************************
	K=K+1
	COE_NUMERATOR=0.0
	DO I=1,N_FIX
        COE_NUMERATOR=COE_NUMERATOR+(U(I)-ULAST(I))**2
    END DO
	DO I=1,SNP
		COE_NUMERATOR=COE_NUMERATOR+(G(I)-GLAST(I))**2
	ENDDO

	COE_DENOMINATOR=0.0
	DO I=1,N_FIX
        COE_DENOMINATOR=COE_DENOMINATOR+U(I)**2
    END DO
	DO I=1,SNP
		COE_DENOMINATOR=COE_DENOMINATOR+G(I)**2
	ENDDO
	COE_CONVER=COE_NUMERATOR/COE_DENOMINATOR
	IF(MOD(K,5)==0 .AND. COE_CONVER>1E-8 )PRINT*,'The cycle is running:',K,'  Convergent coefficient=',COE_CONVER
	IF(COE_CONVER<=1E-8)PRINT*,'The cycle is running:',K,'  Convergent coefficient=',COE_CONVER
	WRITE(19,'(I6,2f15.8)')K,COE_CONVER
	IF(K==1)COE_LAST=COE_CONVER
	IF(COE_CONVER > COE_LAST .AND. K > 1)SOR=SOR*0.5
	COE_LAST=COE_CONVER
	IF(K >= 500)EXIT
ENDDO
CLOSE(19)

END SUBROUTINE MixP
