      SUBROUTINE WRIT_DEVICE
     ;(NDEVS,ND,IODEVICE,BUDAT,NUMTNOD,WTOBSU,WTOBSBU,MAINF,NUMTU,NUF,
     ;NUMTBU,NOBUF,IOBUTYP,NO,NUMTOBSC,VOBS,TOBS,COVINV,IDIMCOV,
     ;NUMTOBS,IOTINT,WTOBSN,INDEXNOD,NUMTNODC,IOUTYP,X,Y,Z,NUMNP)

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION IODEVICE(NDEVS+1,9),BUDAT(4,NUMTNOD),IOBUTYP(NUMTNOD),
     ;WTOBSU(NUMTNOD),WTOBSBU(NUMTNOD),VOBS(NUMTOBS),TOBS(NUMTOBS,2),
     ;COVINV(IDIMCOV),IOTINT(NUMTOBS),NOBUF(NUMTNOD+1),WTOBSN(NUMTNOD),
     ;INDEXNOD(NUMTNOD),IOUTYP(NUMTNOD),X(NUMNP),Y(NUMNP),Z(NUMNP)

*_______________________Writes group E3

      IF (IODEVICE(ND,3).EQ.1) THEN
        WRITE(MAINF,3000)BUDAT(1,NUMTBU),BUDAT(2,NUMTBU),BUDAT(3,NUMTBU)
 3000   FORMAT(5X,'POINT LOCATION',/,5X,3G17.8,//,5X
     ;        ,'NODES USED TO DESCRIBE POINT',/
     ; ,5X,'--------------------------',/,5X,
     ;'OBS_NODE GRID_NODE         X         Y         ',
     ;'     Z          BASIS FUNC.')
        SUM=0
        DO I=IODEVICE(ND,7),NUMTNODC
          J=INDEXNOD(I)
          WRITE(MAINF,3020)I,J,X(J),Y(J),Z(J),WTOBSN(I)
 3020     FORMAT(8X,I5,5X,I5,3G14.7,3X,G15.8)
          SUM=SUM+WTOBSN(I)
        ENDDO
        WRITE(MAINF,3025) SUM
 3025   FORMAT(63X,'SUM',G15.8)
      ENDIF

      IF (IODEVICE(ND,3).GE.2) THEN

        WRITE(MAINF,3030)
 3030   FORMAT(5X,'NODES USED TO DESCRIBE DEVICE',/,5X,'--------------
     ;---------------',/,5X,'OBS_NODE GRID_NODE         X         Y'
     ;,'         Z     N_WEIGHT')
        SUM=0
        DO I=IODEVICE(ND,7),NUMTNODC
          J=INDEXNOD(I)
          WRITE(MAINF,3040)I,J,X(J),Y(J),Z(J),WTOBSN(I)
 3040     FORMAT(8X,I5,5X,I5,F10.3,F10.3,F10.3,3X,F10.3)
          SUM=SUM+WTOBSN(I)
        ENDDO
        WRITE(MAINF,3050)SUM
 3050   FORMAT(53X,'SUM',F10.3)

        WRITE(MAINF,3100)
 3100   FORMAT(5X,/,5X,'   BU BU_TYPE BU_WEIGHT    U   U_TYPE  U_WEIGHT'
     ;,/,5X,'-----------------------------------------------')
        DO I=NUF,NUMTU
          DO J=NOBUF(I),NOBUF(I+1)-1
            WRITE(MAINF,3300) J,IOBUTYP(I),WTOBSBU(J),I,IOUTYP(I),
     ;WTOBSU(I)
 3300       FORMAT(5X,I5,3X,I5,F10.3,I5,4X,I5,5F10.3)
          ENDDO
        ENDDO
      ENDIF

      WRITE(MAINF,3390) ND,NUMTNODC-IODEVICE(ND,7)+1
 3390 FORMAT(/,5X,'TOTAL NUMBER OF NODES DEFINING DEVICE',I5,6X,':',I5)
      WRITE(MAINF,3400) ND,NUMTBU-NOBUF(NUF)+1
 3400 FORMAT(5X,'TOTAL NUMBER OF BASIC UNITS DEFINING DEVICE',I5,':',I5)
      WRITE(MAINF,3500) ND,NUMTU-NUF+1
 3500 FORMAT(5X,'TOTAL NUMBER OF UNITS DEFINING DEVICE',I5,6X,':',I5)


*_______________________Writes card E4.1

      WRITE(MAINF,3505)
 3505 FORMAT(/,5X,'NO VOBS      STDEVOBS  ',
     ;'TOBS      TOBSEND   IOTINT')
      DO NOB=NUMTOBSC+1,NO
         WRITE(MAINF,3510) NOB,VOBS(NOB),COVINV(NOB),TOBS(NOB,1),
     ;TOBS(NOB,2),IOTINT(NOB)
 3510   FORMAT(2X,I5,4E10.3,I5)
      ENDDO


      RETURN
      END
