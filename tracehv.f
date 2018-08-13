      REAL*8 FUNCTION TRACEHV(HESSINV,STPAR,IDIMHESS,NOF,NOL,
     ;IVPAR,NZPAR)
      IMPLICIT NONE
      
************************************************************
* PURPOSE:
* to calculate the trace of the matrix product (inverse of hessian)*covariance
* matrix, where the covariance matrixis the matrix of covariances associated 
* to a single parameter type. 
*                                                     external variables: scalars
      INTEGER*4 IDIMHESS, NOF,NOL,NZPAR
*                                                      external variables: arrays
      REAL*8  HESSINV(IDIMHESS),STPAR(NZPAR)
      INTEGER*4 IVPAR(NZPAR)
*                                                     internal variables: scalars
      INTEGER*4 I,POS
      REAL*8 LCOVINV
*                                                      internal variables: arrays
      
      TRACEHV = 0


      DO I=NOF,NOL
        IF (IVPAR(I) .GT. 0) THEN
          POS=IVPAR(I)*(IVPAR(I)+1)/2
          LCOVINV = 1D0/STPAR(I)/STPAR(I)
          TRACEHV = TRACEHV + HESSINV(POS)*LCOVINV
        ENDIF
      ENDDO
          
      RETURN
      END
