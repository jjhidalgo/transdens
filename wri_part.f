      SUBROUTINE WRI_PART
     ;(IDPARM     ,IDSTPAR  ,IDCFPAREL
     ;,IDCFPARNP ,IDPARNP    ,IDPAREL  ,IDCOVPAR 
     ;,IDSOURCE  ,IDGRAV     ,IDGRAVEL
     ;,IDAFLU    ,IDBFLU     ,IDHCALAN
     ;,IDDFLU    ,IDALFA     ,IDHAUX1  ,IDHAUX2
     ;,IDCFLU    ,IDAFLUDSC  ,IDAFLUDSCF
     ;,IDDERH    ,IDBM_ND_FL ,IDBM_ZN_FL                     
     ;!------------------------------------------ Inverse problem general arrays
     ;,IDGRAD    ,IDHESSAUX  ,IDHESS   ,IDPAR
     ;,IDPARAUX   ,IDCOVINV ,IDVJAC   ,IDPARC
     ;,IDWGT_PAR  ,IDDERIV  ,IDWGT_UNK ,IDPARGOOD
     ;!---------------------------Coordinates and finite element integral arrays
     ;
     ;,IDCOORD   ,IDBIBI     ,IDVOLNOD  ,IDAREA
     ;,IDGRDFF
     ;!-------------------------------------Auxiliar arrays and Real time arrays
     ;,IDWORK     ,IDFNT     ,IDTIME     
     ;,IDCAUDAL   ,IDCONCFLOW
     ;!------------------------------------ Darcy velocity and related variables
     ;,IDVD      ,IDQXYZ     ,IDXNORVD
     ;,IDDQDFLU  ,IDDQDTRA   ,IDDVDH    ,IDDVDC
     ;!-----------------------------------------------Consistent velocity arrays
     ;,IDPOINTWEIGHT,IDGRADLOC
     ;,IDGP_COORD
     ;!---------------------------------------------------------transport arrays
     ;,IDWATVOL  ,IDDWDH     ,IDACTH            
     ;,IDDERC    ,IDCAUX1    ,IDCAUX2    ,IDDVDP        
     ;,IDATRA    ,IDATRADSC  ,IDATRADSCF ,IDBTRA              
     ;,IDDTRA    ,IDCCALAN   ,IDBM_ZN_TT ,IDBM_ND_TT
     ;
     ;!------------------------------------------------------ Observation arrays
     ;,IDVOBSC   ,IDVOBS     ,IDDVOBS  ,IDTOBS	
     ;,IDTIT     ,IDBUDAT    ,IDEXTNBU ,IDWTOBSBU                       
     ;,IDWTOBSU  ,IDWTOBSN   ,IDWTOBST                          
     ;!---------------------------------------------------------- FLOW nonlinear
     ;,IDHBASE   ,IDHPREV1   ,IDDFLUDFLU,IDDBFLUDFLU
     ;,IDHPREV2  ,IDDNODALRH ,IDDPARELDH
     ;,IDHINI    ,IDDAT_VD
     ;,IDCBASE   ,IDCPREV1    ,IDCPREV2     
     ;,IDCCALIT  ,IDDTRADTRA
     ;!-----------------------------------------------Coupled flow and transport
     ;,IDA_COUPL_DSC,IDA_COUPL_DSCF,IDDFLUDTRA,IDDTRADFLU,IDDPARELDC
     ;,IDBCOUPLED      ,IDDBFLUDTRA
     ;,IDDERVISC
     ;
     ;!------------------------- Common non linear variables (flow or transport
     ;,IDHCALIT  ,IDDTMXDS     ,IDPARACD
     ;
     ;!------------------------------------------------density dependency arrays
     ;,IDDENSITY ,IDVISCOSITY  ,IDDELTAITER
     ;
     ;!--------------------------------------------- Statistical analysis arrays
     ;,IDDEVICESTAT,IDEIGENVEC,IDRESID    ,IDRESIDPAR
     ;
     ;!----------------------------------------------------------- String arrays
     ;,IDDEVNAME
     ;
     ;!================================================================
     ;!--------------- Integer arrays
     ;,IDIBCOD  ,IDISOZ    ,IDNFTPAR    ,IDNFNLPAR          
     ;,IDIVPAR  ,IDIDMBLCVP
     ;,IDLXPAREL   ,IDIXPARNP
     ;
     ;!--------------------------------------------------- Discretization arrays
     ;,IDKXX     ,IDLTYPE   ,IDLNNDEL    ,IDLDIM
     ;,IDIBTCO   ,IDINDPAR
     ;!----------------------------------------------------Sparse storage arrays
     ;,IDIAD_S   ,IDIADD_S  ,IDIADN_S    ,IDIAD_D
     ;,IDIADD_D  ,IDIADN_D  ,IDIAFD_S    ,IDIAFDD_S
     ;,IDIAFDN_S ,IDIAFD_D  ,IDIAFDD_D   ,IDIAFDN_D
     ;
     ;!---------------- Integer arrays related to measurements and to dimensions
     ;!---------------- NUMTNOD, NUMTOBS and NUMTIT
     ;,IDIODEVICE,IDINDEXNOD ,IDIOBUTYP  ,IDIOCALBU
     ;,IDIOUTYP  ,IDNBUW     ,IDNOBUF    ,IDNOOBSIT
     ;,IDIOTINT
     ;
     ;!--------------- Time arrays & Auxiliary array & Non linear integer arrays
     ;,IDKINT    ,IDISOLEQ  ,IDINTAUX, IDNFNLTIP
     ;,IDNFNLPRG ,IDLCOORD,IDIPARTNER
     ;
     ;!-------------------------------------------- Statistical analysis arrays
     ;
     ;,IDOBSCLASS,IDITYPEPAR,IDIZPAR
     ;
     ;!------------------------------------------------------------------others
     ;,LASTII      ,LASTIR       ,MAINF)

*****************************************************************************
* PURPOSE
*     Writes location of all partitioned variables in arrays RV and IV
*
* DESCRIPTION
*     Writes location of all partitioned variables in arrays RV and IV
*
* EXTERNAL VARIABLES: SCALARS
*
*  ID*******              Position of variable *******

* HISTORY
*
*     AMS      4-1997     First coding
*     AMS      1-1998     Revision and addition of header
*     AAR     10-2001     Revision and inclusion of obs. related variables
*     AAR      4-2002     Revision and inclusion of geoestatistics related var.
*     AAR      2-2003     Revision
*
*****************************************************************************

      IMPLICIT NONE

      INTEGER*4
     ; IDPARC    ,IDPARM     ,IDSTPAR  ,IDCFPAREL
     ;,IDCFPARNP ,IDPARNP    ,IDPAREL  ,IDCOVPAR 
     ;,IDSOURCE  ,IDGRAV   ,IDGRAVEL
     ;,IDAFLU    ,IDBFLU     ,IDHCALAN
     ;,IDDFLU    ,IDALFA     ,IDHAUX1  ,IDHAUX2
     ;,IDCFLU    ,IDAFLUDSC,IDAFLUDSCF
     ;,IDDERH    ,IDBM_ND_FL ,IDBM_ZN_FL                    
     ;,IDGRAD    ,IDHESSAUX  ,IDHESS   ,IDPAR
     ;,IDPARAUX   ,IDCOVINV ,IDVJAC
     ;,IDCOORD   ,IDBIBI     ,IDVOLNOD  ,IDAREA
     ;,IDGRDFF
     ;,IDWORK     ,IDFNT     ,IDTIME     
     ;,IDCAUDAL   ,IDCONCFLOW
     ;,IDVD      ,IDQXYZ     ,IDXNORVD
     ;,IDDQDFLU,IDDQDTRA,IDDVDH,IDDVDC
     ;,IDPOINTWEIGHT,IDGRADLOC
     ;,IDGP_COORD
     ;,IDWATVOL  ,IDDWDH     ,IDACTH             
     ;,IDDERC    ,IDCAUX1    ,IDCAUX2  ,IDDVDP        
     ;,IDATRA    ,IDATRADSC  ,IDATRADSCF,IDBTRA               
     ;,IDDTRA    ,IDCCALAN   ,IDBM_ZN_TT ,IDBM_ND_TT
     ;,IDVOBSC   ,IDVOBS     ,IDDVOBS  ,IDTOBS	
     ;,IDTIT     ,IDBUDAT    ,IDEXTNBU ,IDWTOBSBU                       
     ;,IDWTOBSU  ,IDWTOBSN   ,IDWTOBST                          
     ;,IDHBASE   ,IDHPREV1   ,IDDFLUDFLU,IDDBFLUDFLU
     ;,IDHPREV2  ,IDDNODALRH ,IDDPARELDH
     ;,IDHINI    ,IDDAT_VD
     ;,IDCBASE   ,IDCPREV1    ,IDCPREV2     
     ;,IDCCALIT  ,IDDTRADTRA
     ;,IDA_COUPL_DSC,IDA_COUPL_DSCF,IDDFLUDTRA,IDDTRADFLU,IDDPARELDC
     ;,IDBCOUPLED
     ;,IDDBFLUDTRA,IDDERVISC
     ;,IDHCALIT  ,IDDTMXDS     ,IDPARACD
     ;,IDDENSITY ,IDVISCOSITY,IDDELTAITER
     ;,IDVARIO    ,IDROTMAT   ,IDPOS_MEAS ,IDPOS_MEAS2
     ;,IDPOS_PIPO ,IDVAR      ,IDVAR2     ,IDSECAT
     ;,IDPOS_BL   ,IDPOS_DIS  ,IDPOS_AUX  ,IDVSTATS
     ;,IDTMAX     ,IDSEARCH   ,IDA        ,IDR
     ;,IDRR       ,IDS        ,IDCLOSE    ,IDTMP
     ;,IDAUXKRIG  ,IDBL_WGT   ,IDCROSS_COV
     ;,IDSIMCOND  ,IDDIV_BL   ,IDPOSDAT_SC
     ;,IDDATA_SC  ,IDORDER
     ;,IDX,IDY,IDZ,IDTPREVINV
     ;,IDDEVICESTAT,IDEIGENVEC,IDRESID
     ;,IDDEVNAME
     ;,IDIBCOD  ,IDISOZ    ,IDNFTPAR    ,IDNFNLPAR          
     ;,IDIVPAR  ,IDIDMBLCVP,IDLXPAREL   ,IDIXPARNP
     ;,IDKXX     ,IDLTYPE   ,IDLNNDEL    ,IDLDIM
     ;,IDIBTCO   ,IDINDPAR
     ;,IDIAD_S   ,IDIADD_S  ,IDIADN_S    ,IDIAD_D
     ;,IDIADD_D  ,IDIADN_D  ,IDIAFD_S    ,IDIAFDD_S
     ;,IDIAFDN_S ,IDIAFD_D  ,IDIAFDD_D   ,IDIAFDN_D
     ;,IDIODEVICE,IDINDEXNOD ,IDIOBUTYP  ,IDIOCALBU
     ;,IDIOUTYP  ,IDNBUW     ,IDNOBUF    ,IDNOOBSIT
     ;,IDIOTINT
     ;,IDKINT    ,IDISOLEQ  ,IDINTAUX, IDNFNLTIP
     ;,IDNFNLPRG ,IDLCOORD  ,IDIPARTNER
     ;,IDIVARIO     ,IDIDRIF     ,IDIVA    ,IDIXSBTOSR
     ;,IDIYSBTOSR   ,IDIZSBTOSR  ,IDNISB   ,IDIBL_PP
     ;,IDIBL_N_WGT  ,IDNSAMPLE   ,IDINDVAR
     ;,IDOBSCLASS,IDITYPEPAR,IDIZPAR,IDRESIDPAR 
     ;,LASTII      ,LASTIR       ,MAINF
     ;,IDWGT_PAR  ,IDDERIV  ,IDWGT_UNK ,IDPARGOOD
C--------------- Initialize to -1 some variables
C--------------- not passed in the subroutine call.

      IDVARIO = -1
      IDROTMAT = -1
      IDVAR = -1
      IDSECAT = -1
      IDPOS_BL = -1
      IDPOS_DIS = -1
      IDPOS_AUX = -1
      IDTMAX = -1
      IDSEARCH = -1
      IDA = -1
      IDR = -1
      IDRR = -1
      IDCLOSE = -1
      IDTMP =-1
      IDAUXKRIG = -1
      IDSIMCOND  = -1
      IDDIV_BL = -1
      IDPOSDAT_SC = -1
      IDORDER = -1
      IDX = -1
      IDY = -1
      IDZ = -1
      IDTPREVINV = -1
      IDIVARIO = -1
      IDIDRIF = -1
      IDIVA = -1
      IDIXSBTOSR = -1
      IDIYSBTOSR = -1
      IDIZSBTOSR = -1
      IDNISB = -1
      IDIBL_PP = -1
      IDIBL_N_WGT = -1
      IDNSAMPLE = -1
      IDINDVAR = -1
      IDS = -1
      IDPOS_MEAS = -1
      IDPOS_MEAS2 = -1
      IDVAR2 = -1
      IDPOS_PIPO = -1
      IDVSTATS = -1
      IDBL_WGT = -1
      IDCROSS_COV = -1
      IDDATA_SC = -1
C--------------- Writes main header

       WRITE(MAINF,1500) 

 1500  FORMAT(////,10X,'PARTITION LIST',/,
     ;             10X,'========= ====',//)

C--------------- Writes all locations in array RV

       WRITE(MAINF,2000)
 2000  FORMAT(//,1X,'ARRAY RV',//)

       WRITE(MAINF,'(A16,I10)')
     ; '  STPAR         ',IDSTPAR     ,'  CFPAREL       ',IDCFPAREL
     ;,'  CFPARNP       ',IDCFPARNP   ,'  PARNP         ',IDPARNP
     ;,'  PAREL         ',IDPAREL     ,'  COVPAR        ',IDCOVPAR
     ;,'  SOURCE        ',IDSOURCE
     ;,'  GRAV          ',IDGRAV      ,'  GRAVEL        ',IDGRAVEL
     ;,'  AFLU          ',IDAFLU      ,'  BFLU          ',IDBFLU
     ;,'  HCALIT        ',IDHCALIT    ,'  HCALAN        ',IDHCALAN
     ;,'  DFLU          ',IDDFLU      ,'  ALFA          ',IDALFA
     ;,'  HAUX1         ',IDHAUX1     ,'  HAUX2         ',IDHAUX2
     ;,'  CFLU          ',IDCFLU
     ;,'  AFLUDSC       ',IDAFLUDSC   ,'  AFLUDSCF      ',IDAFLUDSCF  
     ;,'  DERH          ',IDDERH      
     ;,'  BM_ND_FL      ',IDBM_ND_FL  ,'  BM_ZN_FL      ',IDBM_ZN_FL  
     ;,'  GRAD          ',IDGRAD
     ;,'  HESSAUX       ',IDHESSAUX   ,'  HESS          ',IDHESS
     ;,'  PAR           ',IDPAR
     ;,'  PARAUX        ',IDPARAUX    ,'  COVINV        ',IDCOVINV
     ;,'  VJAC          ',IDVJAC
     ;,'  PARC          ',IDPARC      ,'  PARM          ',IDPARM
     ;,'  WGT_PAR       ',IDWGT_PAR   ,'  DERIV         ',IDDERIV
     ;,'  WGT_UNK       ',IDWGT_UNK   ,'  PARGOOD       ',IDPARGOOD
     ;,'  COORD         ',IDCOORD
     ;,'  BIBI          ',IDBIBI      ,'  VOLNOD        ',IDVOLNOD
     ;,'  AREA          ',IDAREA      ,'  GRDFF         ',IDGRDFF
     ;,'  WORK          ',IDWORK
     ;,'  FNT           ',IDFNT       ,'  TIME          ',IDTIME
     ;,'  CAUDAL        ',IDCAUDAL    ,'  CONCFLOW      ',IDCONCFLOW
     &,'  VD            ',IDVD
     ;,'  QXYZ          ',IDQXYZ      ,'  XNORVD        ',IDXNORVD
     ;,'  DQDFLU        ',IDDQDFLU    ,'  DQDTRA        ',IDDQDTRA
     ;,'  DVDH          ',IDDVDH      ,'  DVDC          ',IDDVDC
     ;,'  POINTWEIGHT   ',IDPOINTWEIGHT,' GRADLOC       ',IDGRADLOC
     ;,'  GP_COORD      ',IDGP_COORD 
     ;,'  WATVOL        ',IDWATVOL    ,'  DWDH          ',IDDWDH
     ;,'  ACTH          ',IDACTH      ,'  ATRA          ',IDATRA      
     ;,'  BTRA          ',IDBTRA      ,'  CCALIT        ',IDCCALIT    
     ;,'  CCALAN        ',IDCCALAN    ,'  DTRA          ',IDDTRA
     ;,'  CAUX1         ',IDCAUX1     ,'  CAUX2         ',IDCAUX2     
     ;,'  ATRADSC       ',IDATRADSC   ,'  ATRADSCF      ',IDATRADSCF
     ;,'  DERC          ',IDDERC      ,'  DVDP          ',IDDVDP      
     ;,'  BM_ND_TT      ',IDBM_ND_TT  ,'  BM_ZN_TT      ',IDBM_ZN_TT  
     ;,'  VOBS          ',IDVOBS      ,'  VOBSC         ',IDVOBSC
     ;,'  DVOBS         ',IDDVOBS     ,'  TOBS          ',IDTOBS
     ;,'  TIT           ',IDTIT       ,'  BUDAT         ',IDBUDAT
     ;,'  EXTNBU        ',IDEXTNBU    ,'  WTOBSBU       ',IDWTOBSBU
     ;,'  WTOBSU        ',IDWTOBSU    ,'  WTOBSN        ',IDWTOBSN
     ;,'  WTOBST        ',IDWTOBST    ,'  HBASE         ',IDHBASE
     ;,'  HPREV1        ',IDHPREV1    ,'  HPREV2        ',IDHPREV2    
     ;,'  DFLUDFLU      ',IDDFLUDFLU  ,'  DBFLUDFLU     ',IDDBFLUDFLU 
     ;,'  DBFLUDTRA     ',IDDBFLUDTRA ,'  DERVISC       ',IDDERVISC
     ;,'  DPARELDH      ',IDDPARELDH  ,'  DNODALRH      ',IDDNODALRH  
     ;,'  HINI          ',IDHINI      ,'  DAT_VD        ',IDDAT_VD
     ;,'  CBASE         ',IDCBASE     ,'  CPREV1        ',IDCPREV1
     ;,'  CPREV2        ',IDCPREV2
     ;,'  DTRADTRA      ',IDDTRADTRA  ,'  DPARELDC      ',IDDPARELDC
     ;,'  A_COUPL_DSC   ',IDA_COUPL_DSC
     ;,'  A_COUPL_DSCF  ',IDA_COUPL_DSCF,'  DFLUDTRA      ',IDDFLUDTRA
     ;,'  DTRADFLU      ',IDDTRADFLU 
     ;,'  BCOUPLED      ',IDBCOUPLED   ,'  DTMXDS        ',IDDTMXDS
     ;,'  PARACD        ',IDPARACD     ,'  DENSITY       ',IDDENSITY
     ;,'  VISCOSITY     ',IDVISCOSITY  ,'  DELTAITER     ',IDDELTAITER
     ;,'  VARIO         ',IDVARIO     ,'  ROTMAT        ',IDROTMAT    
     ;,'  POS_MEAS      ',IDPOS_MEAS  ,'  POS_MEAS2     ',IDPOS_MEAS2 
     ;,'  POS_PIPO      ',IDPOS_PIPO  ,'  VAR           ',IDVAR       
     ;,'  VAR2          ',IDVAR2      ,'  SECAT         ',IDSECAT     
     ;,'  POS_BL        ',IDPOS_BL    ,'  POS_DIS       ',IDPOS_DIS   
     ;,'  POS_AUX       ',IDPOS_AUX   ,'  VSTATS        ',IDVSTATS    
     ;,'  TMAX          ',IDTMAX      ,'  SEARCH        ',IDSEARCH    
     ;,'  A             ',IDA         ,'  R             ',IDR         
     ;,'  RR            ',IDRR        ,'  IDS           ',IDS
     ;,'  CLOSE         ',IDCLOSE     ,'  TMP           ',IDTMP       
     ;,'  AUXKRIG       ',IDAUXKRIG  
     ;,'  BL_WGT        ',IDBL_WGT    ,'  CROSS_COV     ',IDCROSS_COV 
     ;,'  SIMCOND       ',IDSIMCOND   
     ;,'  DIV_BL        ',IDDIV_BL    ,'  POSDAT_SC     ',IDPOSDAT_SC 
     ;,'  DATA_SC       ',IDDATA_SC   ,'  ORDER         ',IDORDER
     ,,'  X             ',IDX         ,'  Y             ',IDY
     ;,'  Z             ',IDZ         ,'  TPREVINV      ',IDTPREVINV     
     ;,'  DEVICESTAT    ',IDDEVICESTAT,'  EIGENVEC      ',IDEIGENVEC  
     ;,'  RESID         ',IDRESID     ,'  RESIDPAR      ',IDRESIDPAR    
     ;,'  LASTIR        ',LASTIR
               
   
    
C--------------- Writes all locations in array IV

       WRITE(MAINF,2100)
 2100  FORMAT(//,1X,'ARRAY IV',//)

       WRITE(MAINF,'(A16,I10)')
     ; '  IBCOD         ',IDIBCOD       ,'  ISOZ          ',IDISOZ
     ;,'  NFTPAR        ',IDNFTPAR      ,'  NFNLPAR       ',IDNFNLPAR
     ;,'  IVPAR         ',IDIVPAR       ,'  IDMBLCVP      ',IDIDMBLCVP
     ;,'  LXPAREL       ',IDLXPAREL     ,'  IXPARNP       ',IDIXPARNP
     ;,'  KXX           ',IDKXX         ,'  LTYPE         ',IDLTYPE
     ;,'  LNNDEL        ',IDLNNDEL      ,'  LDIM          ',IDLDIM
     ;,'  IAD_S         ',IDIAD_S       ,'  IADD_S        ',IDIADD_S
     ;,'  IADN_S        ',IDIADN_S      ,'  IAFD_S        ',IDIAFD_S      
     ;,'  IAFDD_S       ',IDIAFDD_S     ,'  IAFDN_S       ',IDIAFDN_S     
     ;,'  IAD_D         ',IDIAD_D
     ;,'  IADD_D        ',IDIADD_D      ,'  IADN_D        ',IDIADN_D
     ;,'  IAFD_D        ',IDIAFD_D
     ;,'  IAFDD_D       ',IDIAFDD_D     ,'  IAFDN_D       ',IDIAFDN_D
     ;,'  IBTCO         ',IDIBTCO       ,'  INDPAR        ',IDINDPAR
     ;,'  IODEVICE      ',IDIODEVICE    ,'  INDEXNOD      ',IDINDEXNOD
     ;,'  IOBUTYP       ',IDIOBUTYP     ,'  IOCALBU       ',IDIOCALBU
     ;,'  IOUTYP        ',IDIOUTYP      ,'  NBUW          ',IDNBUW
     ;,'  NOBUF         ',IDNOBUF       ,'  NOOBSIT       ',IDNOOBSIT
     ;,'  IOTINT        ',IDIOTINT      ,'  KINT          ',IDKINT
     ;,'  ISOLEQ        ',IDISOLEQ      ,'  INTAUX        ',IDINTAUX
     ;,'  NFNLTIP       ',IDNFNLTIP     ,'  NFNLPRG       ',IDNFNLPRG
     ;,'  LCOORD        ',IDLCOORD      ,'  IPARTNER      ',IDIPARTNER
     ;,'  IVARIO        ',IDIVARIO      ,'  IDRIF         ',IDIDRIF
     ;,'  IVA           ',IDIVA         ,'  IXSBTOSR      ',IDIXSBTOSR
     ;,'  IYSBTOSR      ',IDIYSBTOSR    ,'  IZSBTOSR      ',IDIZSBTOSR
     ;,'  NISB          ',IDNISB        ,'  IBL_PP        ',IDIBL_PP
     ;,'  IBL_N_WGT     ',IDIBL_N_WGT   ,'  NSAMPLE       ',IDNSAMPLE
     ;,'  INDVAR        ',IDINDVAR      ,'  OBSCLASS      ',IDOBSCLASS
     ;,'  ITYPEPAR      ',IDITYPEPAR    ,'  IZPAR         ',IDIZPAR
     ;,'  LASTII        ',LASTII
     ;


C--------------- Writes all locations in array KV

       WRITE(MAINF,2200)
 2200  FORMAT(//,1X,'ARRAY KV',//)

       WRITE(MAINF,'(A16,I10)')
     ;,'  DEVNAME         ',IDDEVNAME

      RETURN
      END
