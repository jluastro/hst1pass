 
 
c----------------------------------------------
c
c
c hst2collate
c
c this is an early, quick-and-dirty version of hst2collate
c
c hst2collate REF=stem_drz.UVW
c             *_flt.xympquvmUVM
c             MAT+
c             LNK+
c             ZPU=EXPT100
c
c
c ---> outputs: mat files --> mapping info
c               lnk files --> time series
c               avg file  --> xbar, ybar, mbar, xrms, yrms, mrms, xerr, yerr, merr,
c                             msat, muns, np, nf, ng, wg, RA DEC
c
c (c) copyright 2022 Jay Anderson, STScI
c
c--------------------------------------------------
 
 
c--------------------------------------------------
c
c This is version 1a.. the first release on 2022.07.05
c
c--------------------------------------------------
 
 
#define _PI_    3.141592653589793238462643d0
#define _LINUX_ .true.
 
#define _MMAX_  500   /* max number of hst images to collate; can be upped */
#define _GMAX_    1   /* max number of GAIA stars in HST field             */
 
 
 
      program hst2collate
      implicit none
 
      real*8,    dimension(:,:),   allocatable :: xr_nm
      real*8,    dimension(:,:),   allocatable :: yr_nm
      real*8,    dimension(:,:),   allocatable :: mr_nm
      real*8,    dimension(:,:),   allocatable :: qr_nm
      real*8,    dimension(:,:),   allocatable :: xc_nm
      real*8,    dimension(:,:),   allocatable :: yc_nm
      real*8,    dimension(:,:),   allocatable :: mc_nm
      real*8,    dimension(:,:),   allocatable :: uc_nm
      real*8,    dimension(:,:),   allocatable :: vc_nm
      real*8,    dimension(:,:),   allocatable :: ra_nm
      real*8,    dimension(:,:),   allocatable :: de_nm
 
      real*8,    dimension(:,:),   allocatable :: ut_nm
      real*8,    dimension(:,:),   allocatable :: vt_nm
      real*8,    dimension(:,:),   allocatable :: mt_nm
 
      real*8,    dimension(:), allocatable :: xref_r
      real*8,    dimension(:), allocatable :: yref_r
      real*8,    dimension(:), allocatable :: mref_r
 
      integer*2,      dimension(:,:), allocatable :: cvg_map
      integer*2,      dimension(:,:), allocatable :: min_map
      integer*2,      dimension(:,:), allocatable :: pk1_map
      integer*2,      dimension(:,:), allocatable :: pk2_map
      integer*2,      dimension(:,:), allocatable :: pk3_map
 
      character*200, dimension(:,:), allocatable :: s200_nm
      character*200  s200h
      character*200  s200o
      integer        strlen
 
      integer NMAX
      integer MMAX
 
      integer M, Ms
      integer N, Ns, Nu, Nx, NLINEs
 
      integer NARG, NARGs
      character*200 STRING_ARG(_MMAX_)
 
      character*200 FILE_REF
      character*200 FILE_LST
      character*20  COLS_REF
      character*20  STEM_REF
      character*80  MATFILE
      character*80  LNKFILE
 
      character*200    FILE_M(_MMAX_)
      character*20     COLS_M(_MMAX_), COLSu
      character*20     STEM_M(_MMAX_), STEMu
      real             EXPT_M(_MMAX_)
      integer         NCOLS_M(_MMAX_)
      integer           NEM_M(_MMAX_)
      integer          Ns_M(_MMAX_)
      integer            Ks_M(_MMAX_), K
      real      CRPIX1_INST_M(_MMAX_)
      real      CRPIX2_INST_M(_MMAX_)
      character*200  STDGDC_M(_MMAX_)
 
      character*200 STR200
      character*40  cdata(99)
 
      integer NCOL, NCOLs, c
      integer R, Rs
 
      integer i, isla, idot, iend
      integer j
 
      integer hist(101,101), hmax
      real  dx,  dy,  dd
      real ddx, ddy, ddd
 
      real*8 x1_l(99999), y1_l(99999), m1_l(99999)
      real*8 x2_l(99999), y2_l(99999), m2_l(99999)
      real   dm_l(99999), dmbar, dmsig
      integer uu_l(99999)
 
      integer L, Ls, Lu, Lt, Lm, Lmu
      real*8 AG, BG, CG, DG
      real*8 GA, GB, GC, GD
      real*8 x1o, y1o, x2o, y2o
      real*8 x1e, y1e, x2e, y2e
      real   RMS2D
      real   SIG_CLIP
      real   RES_MIN
      real   RES_MAX
 
      real*8  AG_MI(_MMAX_,5), x1o_MI(_MMAX_,5)
      real*8  BG_MI(_MMAX_,5), y1o_MI(_MMAX_,5)
      real*8  CG_MI(_MMAX_,5), x2o_MI(_MMAX_,5)
      real*8  DG_MI(_MMAX_,5), y2o_MI(_MMAX_,5)
      real*8  GA_MI(_MMAX_,5),  GB_MI(_MMAX_,5)
      real*8  GC_MI(_MMAX_,5),  GD_MI(_MMAX_,5)
      integer  N_MI(_MMAX_,5)
      real*8  ZP_MI(_MMAX_,5)
      real*8 eZP_MI(_MMAX_,5)
 
      integer NIT, NITs
 
      real*8  xe, ye
      real*8  XMIN, XMAX
      real*8  YMIN, YMAX
 
      integer IMIN, IMAX, IDIM
      integer JMIN, JMAX, JDIM
 
      integer NAXIS1_REF, NAXIS2_REF
      real*8  CRPIX1_REF, CRPIX2_REF
      real*8  CRVAL1_REF, CRVAL2_REF
      real*8   CD1_1_REF,  CD1_2_REF
      real*8   CD2_1_REF,  CD2_2_REF
 
      real*8  CRPIX1_GGG, CRPIX2_GGG
      real*8   CD1_1_GGG,  CD1_2_GGG
      real*8   CD2_1_GGG,  CD2_2_GGG
 
      logical inside_poly
      integer P
      real*8  BDRYX_PKM(4,4,_MMAX_), BDRYU_P(4)
      real*8  BDRYY_PKM(4,4,_MMAX_), BDRYV_P(4)
      real*8  xcu, ycu
      real*8  xrefu, yrefu
 
      integer PKMIN
      integer pkx_l(99999)
      integer pky_l(99999)
      integer pkf_l(99999)
      integer pkp_l(99999)
 
      integer, dimension(:,:), allocatable :: n_lm
 
      integer utrans_l(999999)
      real*8    xbar_l(999999)
      real*8    ybar_l(999999)
      real*8    mbar_l(999999)
      real*8    qbar_l(999999)
      real*8    xrms_l(999999)
      real*8    yrms_l(999999)
      real*8    drms_l(999999)
      real*8    mrms_l(999999)
      integer   nuse_l(999999)
 
      real    dmin
      integer nmin
 
      character*8 ZP_USE
 
      integer O, Os, Om, Of, Op, Ou
      real*8  x_o(_MMAX_)
      real*8  y_o(_MMAX_)
      real*8  m_o(_MMAX_)
      real*8  q_o(_MMAX_)
      integer u_o(_MMAX_)
      real*8  xbar, ybar, mbar, qbar
      real*8  xsig, ysig, msig, qsig
      real*8  asig, bsig
      real*8  DKEEP
      real*8  mbar_min
      logical SHOW
 
      real*8 xbarx, ybarx
      real*8 xsigx, ysigx, dsigx
 
      integer G, Gs
      real*8  RA, dRA, xy2r
      real*8  DE, dDE, xy2d
      real*8  r1, r2, d1, d2
      character*140         str140_g(_GMAX_)
      real*8     x_g(_GMAX_),    y_g(_GMAX_),   m_g(_GMAX_)
      real*8    ex_g(_GMAX_),   ey_g(_GMAX_),  em_g(_GMAX_)
      real*8     r_g(_GMAX_),    d_g(_GMAX_),  er_g(_GMAX_),  ed_g(_GMAX_)
      real*8    mg_g(_GMAX_),  emg_g(_GMAX_)
      real*8    mb_g(_GMAX_),  emb_g(_GMAX_)
      real*8    mr_g(_GMAX_),  emr_g(_GMAX_)
      real*8   pmr_g(_GMAX_), epmr_g(_GMAX_)
      real*8   pmd_g(_GMAX_), epmd_g(_GMAX_)
      real*8   par_g(_GMAX_), epar_g(_GMAX_)
      integer*8 id_g(_GMAX_)
 
      integer  U, Us, Uu
      real*8  x1_u(9999), y1_u(9999), m1_u(9999)
      real*8  x2_u(9999), y2_u(9999), m2_u(9999)
      integer uu_u(9999)
      integer gg_u(9999)
      integer ll_u(9999)
 
      real*8 rDAT
      real*8 EXPTu
 
      real*8 xre, yre
      real*8 xce, yce
      real*8 ute, vte, mte
 
      logical INPUT_LIST
      logical FIND_MAT
      logical FIND_LST
 
      INPUT_LIST      = .false.
      FIND_MAT        = .true.
      FIND_LST        = .true.
 
      ZP_USE = 'EXPT1000'
      EXPTu  = 0.
 
      Ms = 0
 
      NITs = 1
 
      NIT = 0
 
      NARGs = iargc()
 
      if (NARGs.lt.2) then
          print*,' '
          print*,'Barebones version of hst2collate takes 2+ args   '
          print*,'                                                 '
          print*,'./hst2collate.e               \                  '
          print*,'    REF=stem_drz.UVW          \                  '
          print*,'    *_flt.xympquvmUVM         \                  '
          print*,'   [ZPU=EXPT1000]             \                  '
          print*,'   [LST=starlist.UVW]         \  (no finding; just  '
          print*,'                                  collate ref stars)'
          print*,'                                                 '
          print*,'It cross matches the REFs UV frame              '
          print*,'to the UV frame of the individual                '
          print*,'_flt exposures and collates the stars in         '
          print*,'the _flt images into a catalog, reporting        '
          print*,'average photometry and RMS residuals for         '
          print*,'all.                                             '
          print*,'                                                 '
          stop
          endif
 
      print*,' '
      print*,'------------------------------------------'
      print*,'NARGS---'
      print*,' '
      do NARG = 1, NARGs
         call getarg(NARG,STRING_ARG(NARG))
         write(*,'('' NARG'',i3.3,1x,a80)') NARG, STRING_ARG(NARG)(1:80)
         enddo
      print*,' '
 
      print*,'NARGs: ',NARGs
      do NARG = 1, NARGs
         if (STRING_ARG(NARG)(1:3).eq.'IN=') then
             open(09,file=STRING_ARG(NARG)(4:80),status='old')
   91        read(09,*,end=92) STRING_ARG(NARGs+1)
             NARGs = NARGs + 1
             goto 91
   92        continue
             endif
         enddo
 
      print*,'NARGs: ',NARGs
 
      do NARG = 1, NARGs
         write(*,'(i3.3,1x,a80)') NARG, STRING_ARG(NARG)
         enddo
 
 
      do NARG = 1, NARGs
         if (STRING_ARG(NARG)(1:3).eq.'IN=') goto 1
         if (STRING_ARG(NARG)(1:9).eq.'FIND_MAT-') then
             FIND_MAT = .false.
             NITs = 0
             goto 1
             endif
         if (STRING_ARG(NARG)(1:4).eq.'LST=') then
             FIND_LST = .false.
             FILE_LST = STRING_ARG(NARG)(5:200)
             goto 1
             endif
         if (STRING_ARG(NARG)(1:4).eq.'REF=') then
             FILE_REF = STRING_ARG(NARG)(5:200)
             goto 1
             endif
         if (STRING_ARG(NARG)(1:4).eq.'ZPU=') then
             ZP_USE = STRING_ARG(NARG)(5:12)
             goto 1
             endif
         Ms = Ms + 1
         if (Ms.gt._MMAX_) stop 'not designed for > _MMAX_ images...'
         FILE_M(Ms) = STRING_ARG(NARG)
    1    continue
         enddo
 
      print*,'ZP_USE: ',ZP_USE
 
      if (ZP_USE(1:4).eq.'EXPT') then
          read(ZP_USE(5:8),*) EXPTu
          print*,'---> EXPTu: ',EXPTu
          endif
 
      idot = 0
      iend = 0
      isla = 0
      do i = 1, 200
         if (FILE_REF(i:i).eq.'.') idot = i
         if (FILE_REF(i:i).ne.' ') iend = i
         if (FILE_REF(i:i).eq.'/') isla = i
         enddo
 
      COLS_REF = FILE_REF(idot+1:iend)
      STEM_REF = FILE_REF(isla+1:idot-1)
 
      ncols = iend-idot
 
c     print*,'FILE_REF: ',FILE_REF(1:40)
c     print*,'  isla: ',isla
c     print*,'  idot: ',idot
c     print*,'  iend: ',iend
c     print*,' ncols: ',ncols
c     print*,'STEM_REF: ',STEM_REF
c     print*,'COLS_REF: ',COLS_REF
 
      NAXIS1_REF = -1
      NAXIS2_REF = -1
      CRPIX1_REF = -1
      CRPIX2_REF = -1
      CRVAL1_REF = -1
      CRVAL2_REF = -1
      CD1_1_REF  =  0.000
      CD1_2_REF  =  0.000
      CD2_1_REF  =  0.000
      CD2_2_REF  =  0.000
 
      open(11,file=FILE_REF,status='old')
      Rs = 0
 111  read(11,'(a200)',end=211) STR200
      if (STR200(1:11).eq.'#   NAXIS1:') read(STR200(12:40),*)NAXIS1_REF
      if (STR200(1:11).eq.'#   NAXIS2:') read(STR200(12:40),*)NAXIS2_REF
      if (STR200(1:11).eq.'#   CRPIX1:') read(STR200(12:40),*)CRPIX1_REF
      if (STR200(1:11).eq.'#   CRPIX2:') read(STR200(12:40),*)CRPIX2_REF
      if (STR200(1:11).eq.'#   CRVAL1:') read(STR200(12:40),*)CRVAL1_REF
      if (STR200(1:11).eq.'#   CRVAL2:') read(STR200(12:40),*)CRVAL2_REF
      if (STR200(1:11).eq.'#    CD1_1:') read(STR200(12:40),*) CD1_1_REF
      if (STR200(1:11).eq.'#    CD1_2:') read(STR200(12:40),*) CD1_2_REF
      if (STR200(1:11).eq.'#    CD2_1:') read(STR200(12:40),*) CD2_1_REF
      if (STR200(1:11).eq.'#    CD2_2:') read(STR200(12:40),*) CD2_2_REF
      if (STR200(1:11).eq.'#    CD2_2:') read(STR200(12:40),*) CD2_2_REF
      if (STR200(1:11).eq.'#     rDAT:') read(STR200(12:40),*) rDAT
      if (STR200(1:1).eq.'#') goto 111
      Rs = Rs + 1
      goto 111
 211  continue
      close(11)
 
      write(*,*)
      write(*,*) '--------------------------------------------'
      write(*,'('' FILE_REF -- '',a80)') FILE_REF(1:80)
      write(*,'('' '')')
      write(*,'('' '')')
      write(*,'(''               Rs: '',i6)')  Rs
      write(*,'('' '')')
      write(*,'(''            NAXIS: '',i6,'' x '',i4)') NAXIS1_REF,
     .                                                   NAXIS2_REF
      write(*,'('' '')')
      write(*,'(''            CRPIX: '',2f15.8)') CRPIX1_REF, CRPIX2_REF
      write(*,'(''            CRVAL: '',2f15.8)') CRVAL1_REF, CRVAL2_REF
      write(*,'(''             CD1_: '',2f15.8)')  CD1_1_REF,  CD1_2_REF
      write(*,'(''             CD2_: '',2f15.8)')  CD2_1_REF,  CD2_2_REF
      write(*,'('' '')')
      write(*,'(''             rDAT: '',f15.8)') rDAT
      write(*,'('' '')')
 
 
 
c------------------------------
c
c read in reference stars
c
 
      allocate(xref_r(Rs))
      allocate(yref_r(Rs))
      allocate(mref_r(Rs))
 
      print*,' '
      print*,'------------------------------------------'
      print*,'READ IN Rs: '
      print*,' '
      open(11,file=FILE_REF,status='old')
      R = 0
 311  read(11,'(a200)',end=411) STR200
      if (STR200(1:1).eq.'#') goto 311
      R = R + 1
      read(STR200,*) (cdata(ncol),ncol=1,NCOLs)
      read(cdata(1),*) xref_r(R)                      ! default read in from first three...
      read(cdata(2),*) yref_r(R)
      read(cdata(3),*) mref_r(R)
      do ncol = 1, ncols, 1
         if (COLS_REF(ncol:ncol).eq.'U') read(cdata(ncol),*) xref_r(R)
         if (COLS_REF(ncol:ncol).eq.'V') read(cdata(ncol),*) yref_r(R)
         if (COLS_REF(ncol:ncol).eq.'W') read(cdata(ncol),*) mref_r(R)
         enddo
      if (R.eq.1.or.
     .    R.eq.Rs.or.
     .    R.eq.R/(Rs/10)*(Rs/10))
     .    write(*,999) R, Rs, xref_r(R), yref_r(R), mref_r(R)
 999  format(5x,i5,' /',i5,4x,f9.3,1x,f9.3,1x,f8.4)
      goto 311
 411  continue
      close(11)
 
      open(69,file='hst2collate.00.info',status='unknown')
      write(69,319)
 
 551  continue
      print*,' '
      print*,'------------------------------------------'
      print*,'COUNT STARS FOR EACH CONTRIBUTING IMAGE...'
      print*,' '
      write(*,319)
 319  format('NIM',1x,'   Ns',1x,' Nmax',1x,'NLINEs',1x,
     .                ' EXPT_M ',1x,'CRPIX1_I',1x,'CRPIX2_I')
      Nx = 0
      do M = 1, Ms
         STDGDC_M(M) = 'NONE'
         CRPIX1_INST_M(M) = 0.000
         CRPIX2_INST_M(M) = 0.000
         Ks_M(M) = 0
         open(12,file=FILE_M(M),status='old')
         Ns = 0
         NLINEs = 0
 112     read(12,'(a200)',end=212) STR200
         NLINEs = NLINEs + 1
         if (STR200(01:11).eq.'#   CRPIX1:') then
             read(STR200(12:200),*) CRPIX1_INST_M(M)
             goto 112
             endif
         if (STR200(01:11).eq.'#   CRPIX2:') then
             read(STR200(12:200),*) CRPIX2_INST_M(M)
             goto 112
             endif
         if (STR200(01:11).eq.'#  GDC_INP:') then
             STDGDC_M(M) = STR200(13:200)
             goto 112
             endif
         if (STR200(01:11).eq.'#    EXPTE:') then
             read(STR200(12:200),*) EXPT_M(M)
             goto 112
             endif
         if (STR200(01:11).eq.'# BDRY_uv_K') then
             read(STR200(12:012),*) K
             read(STR200(14:200),*) BDRYX_PKM(1,K,M), BDRYY_PKM(1,K,M),
     .                              BDRYX_PKM(2,K,M), BDRYY_PKM(2,K,M),
     .                              BDRYX_PKM(3,K,M), BDRYY_PKM(3,K,M),
     .                              BDRYX_PKM(4,K,M), BDRYY_PKM(4,K,M)
             if (K.gt.Ks_M(M)) Ks_M(M) = K
             goto 112
             endif
         if (STR200(1:1).eq.'#') goto 112
         Ns = Ns + 1
         goto 112
 212     continue
         if (Ns.gt.Nx) Nx = Ns
         write(69,312) M,Ns,Nx,NLINEs,EXPT_M(M),
     .                 CRPIX1_INST_M(M),
     .                 CRPIX2_INST_M(M),FILE_M(M)(1:40)
         write( *,312) M,Ns,Nx,NLINEs,EXPT_M(M),
     .                 CRPIX1_INST_M(M),
     .                 CRPIX2_INST_M(M),FILE_M(M)(1:40)
 312     format(1x,i3,1x,i5,1x,i5,1x,i5,1x,f8.1,1x,2f8.1,1x,a40)
         close(12)
         enddo
      close(69)
 
      print*,'         '
      print*,'ALLOCATE --- '
      print*,'   Nx: ',Nx
      print*,'     x   '
      print*,'   Ms: ',Ms
      print*,'         '
 
      allocate(xr_nm(Nx,Ms))
      allocate(yr_nm(Nx,Ms))
      allocate(mr_nm(Nx,Ms))
      allocate(qr_nm(Nx,Ms))
      allocate(mt_nm(Nx,Ms))
 
      allocate(xc_nm(Nx,Ms))
      allocate(yc_nm(Nx,Ms))
      allocate(mc_nm(Nx,Ms))
      allocate(uc_nm(Nx,Ms))
      allocate(vc_nm(Nx,Ms))
      allocate(ut_nm(Nx,Ms))
      allocate(vt_nm(Nx,Ms))
      allocate(ra_nm(Nx,Ms))
      allocate(de_nm(Nx,Ms))
      allocate(s200_nm(Nx,Ms))
 
 
      print*,' '
      print*,'---------------------------------------------'
      print*,'READ IN STARS FROM EACH CONTRIBUTING IMAGE...'
      print*,' '
      do M = 1, Ms
         idot = 0
         iend = 0
         isla = 0
         do i = 1, 200
            if (FILE_M(M)(i:i).eq.'.') idot = i
            if (FILE_M(M)(i:i).ne.' ') iend = i
            if (FILE_M(M)(i:i).eq.'/') isla = i
            enddo
         COLSu = FILE_M(M)(idot+1:iend)
         COLS_M(M) = COLSu
         ncols = iend-idot
         NCOLS_M(M) = ncols
         STEMu = FILE_M(M)(isla+1:idot-1)
         STEM_M(M)  = STEMu
         NEM_M(M) = idot-1-isla
         open(12,file=FILE_M(M),status='old')
         Ns = 0
         Nu = 0
 412     read(12,'(a200)',end=512) STR200
         if (STR200(1:1).eq.'#') goto 412
         Ns = Ns + 1
         read(STR200,*) (cdata(ncol),ncol=1,NCOLs)
         s200_nm(Ns,M) = STR200
         do c = 1, ncols
            if (COLSu(c:c).eq.'x') read(cdata(c),*) xr_nm(Ns,M)
            if (COLSu(c:c).eq.'y') read(cdata(c),*) yr_nm(Ns,M)
            if (COLSu(c:c).eq.'m') read(cdata(c),*) mr_nm(Ns,M)
            if (COLSu(c:c).eq.'q') read(cdata(c),*) qr_nm(Ns,M)
            if (COLSu(c:c).eq.'u') read(cdata(c),*) xc_nm(Ns,M)
            if (COLSu(c:c).eq.'v') read(cdata(c),*) yc_nm(Ns,M)
            if (COLSu(c:c).eq.'r') read(cdata(c),*) ra_nm(Ns,M)
            if (COLSu(c:c).eq.'d') read(cdata(c),*) de_nm(Ns,M)
            if (COLSu(c:c).eq.'U') read(cdata(c),*) uc_nm(Ns,M)
            if (COLSu(c:c).eq.'V') read(cdata(c),*) vc_nm(Ns,M)
            enddo
         if ((.false.).and.Ns.eq.Ns/100*100)
     .   write(*,912) M, Ns,
     .                xr_nm(Ns,M), yr_nm(Ns,M), mr_nm(Ns,M),
     .                xc_nm(Ns,M), yc_nm(Ns,M),
     .                ra_nm(Ns,M), de_nm(Ns,M),
     .                uc_nm(Ns,M), vc_nm(Ns,M)
 912     format(1x,i3.3,1x,i4,
     .          5x,2f10.3,1x,f7.3,
     .          5x,2f10.3,
     .          5x,2f13.7,
     .          5x,2f10.3)
         Ns_M(M) = Ns
         mc_nm(Ns,M) = mr_nm(Ns,M)
         mt_nm(Ns,M) = mr_nm(Ns,M)
         ut_nm(Ns,M) = uc_nm(Ns,M)
         vt_nm(Ns,M) = vc_nm(Ns,M)
         if (EXPTu.gt.0)
     .   mt_nm(Ns,M) = mr_nm(Ns,M)
     .               - 2.5*log10(EXPTu/EXPT_M(M))
         if (xc_nm(Ns,M).gt.0.00.and.
     .       yc_nm(Ns,M).gt.0.00) Nu = Nu + 1
         goto 412
 512     continue
         close(12)
         write(*,'(i3,1x,i5,1x,a40)') M,Ns_M(M),FILE_M(M)
         if (Nu.eq.0) stop 'no stars with (uv) photy read in'
         enddo
      print*,' '
      print*,' '
 
 
 
      open(70,file='hst2collate.01.trans1',status='unknown')
      write(70,370)
      write(70,270)
      write(70,370)
 
      if (.not.FIND_MAT) then
         print*,' '
         print*,'---------------------------------------------'
         print*,'READ IN THE TRANS FOR EACH CONTRIBUTING IMAGE'
         print*,' '
         write( *,370)
         write( *,270)
         write( *,370)
         do M = 1, Ms
            MATFILE = FILE_M(M)(1:9) // '_mat.UVuvWw'
            open(91,file=MATFILE,status='old')
            Ls = 0
  554       continue
            read(91,'(a200)',end=553) STR200
            if (STR200(1:1).eq.'#') goto 554
            Ls = Ls + 1
            read(STR200,*) x1_l(Ls), y1_l(Ls),
     .                     x2_l(Ls), y2_l(Ls),
     .                     m1_l(Ls), m2_l(Ls)
            dm_l(Ls) = m2_l(Ls) - m1_l(Ls)
            uu_l(Ls) = 1
            goto 554
  553       continue
            RMS2D = 0.
            Lu = Ls
            call glob_fit6nrDPU(x1_l,y1_l,x2_l,y2_l,uu_l,Ls,    ! added Nu
     .                          AG,BG,CG,DG,x1o,y1o,x2o,y2o)
            AG_MI(M,1) = AG
            BG_MI(M,1) = BG
            CG_MI(M,1) = CG
            DG_MI(M,1) = DG
           x1o_MI(M,1) = x1o
           y1o_MI(M,1) = y1o
           x2o_MI(M,1) = x2o
           y2o_MI(M,1) = y2o
             N_MI(M,1) = Lu
            GA_MI(M,1) =  DG/(AG*DG-BG*CG)
            GB_MI(M,1) = -BG/(AG*DG-BG*CG)
            GC_MI(M,1) = -CG/(AG*DG-BG*CG)
            GD_MI(M,1) =  AG/(AG*DG-BG*CG)
            call rbarsigs(dm_l,Ls,dmbar,dmsig,Lmu,2.5)
             ZP_MI(M,1) = dmbar
            eZP_MI(M,1) = dmsig
            do N = 1, Ns_M(M)
               ut_nm(N,M) = x1o + GA_MI(M,1)*(xc_nm(N,M)-x2o)
     .                          + GB_MI(M,1)*(yc_nm(N,M)-y2o)
               vt_nm(N,M) = y1o + GC_MI(M,1)*(xc_nm(N,M)-x2o)
     .                          + GD_MI(M,1)*(yc_nm(N,M)-y2o)
               mt_nm(N,M) = mc_nm(N,M) - dmbar
               enddo
            write( *,170) M, Ns_M(M), Ls, Lu,
     .                    AG, BG, CG, DG,
     .                    x1o, y1o, x2o, y2o, RMS2D,
     .                    atan2(BG-CG,AG+DG)*180/3.1415927,
     .                    dmbar, dmsig, Lmu,
     .                    FILE_M(M)(1:40)
            write(70,170) M, Ns_M(M), Ls, Lu,
     .                    AG, BG, CG, DG,
     .                    x1o, y1o, x2o, y2o, RMS2D,
     .                    atan2(BG-CG,AG+DG)*180/3.1415927,
     .                    dmbar, dmsig, Lmu,
     .                    FILE_M(M)(1:40)
            enddo
         close(70)
         goto 555
         endif
 
      print*,' '
      print*,'----------------------------------------------------'
      print*,'FIND THE TRANSFORMATIONS FOR EACH CONTRIBUTING IMAGE'
      print*,' '
      write( *,370)
      write( *,270)
      write( *,370)
 
c     do R = 1, Rs
c        write(41,'(3f10.3)') xref_r(R), yref_r(R), mref_r(R)
c        enddo
c     close(41)
 
 
      do M = 1, Ms
c        do N = 1, Ns_M(1)
c           write(42,'(3f10.3)') ut_nm(N,M), vt_nm(N,M), mt_nm(N,M)
c           enddo
c        close(42)
         do i = 001, 101
         do j = 001, 101
            hist(i,j) = 0
            enddo
            enddo
         do N = 1, Ns_M(M)
         do R = 1, Rs
            i = 51 + (ut_nm(N,M) - xref_r(R))/5
            j = 51 + (vt_nm(N,M) - yref_r(R))/5
            if (i.ge.1.and.i.le.100.and.
     .          j.ge.1.and.j.le.100) then
                hist(i  ,j)   = hist(i  ,j  ) + 1
                hist(i+1,j)   = hist(i+1,j  ) + 1
                hist(i  ,j+1) = hist(i  ,j+1) + 1
                hist(i+1,j+1) = hist(i+1,j+1) + 1
                endif
            enddo
            enddo
 
         imax = 0
         jmax = 0
         hmax = 0
         do i = 1, 101
         do j = 1, 101
            if (hist(i,j).gt.hmax) then
                hmax = hist(i,j)
                imax = i
                jmax = j
                endif
            enddo
            enddo
 
         dx = (imax-51+0.5)*5
         dy = (jmax-51+0.5)*5
 
c        print*,' imax: ',imax,jmax,hmax
c        print*,'   dx: ',dx,dy
 
 
         Ls = 0
         do N = 1, Ns_M(M)
         do R = 1, Rs
            ddx = ut_nm(N,M)-xref_r(R)-dx
            ddy = vt_nm(N,M)-yref_r(R)-dy
            ddd = sqrt(ddx**2+ddy**2)
            if (ddd.lt.10.5) then
                Ls = Ls + 1
                if (Ls.gt.99999) stop 'Ls.gt.99999'
                x1_l(Ls) = xref_r(R)
                y1_l(Ls) = yref_r(R)
                m1_l(Ls) = mref_r(R)
                x2_l(Ls) = xc_nm(N,M)
                y2_l(Ls) = yc_nm(N,M)
                m2_l(Ls) = mt_nm(N,M)
                dm_l(Ls) = mt_nm(N,M)-mref_r(R)
                uu_l(Ls) = 1
c               write(44,144) x1_l(Ls), y1_l(Ls),
c    .                        x2_l(Ls), y2_l(Ls),
c    .                        m1_l(Ls), m2_l(Ls),
c    .                        dm_l(Ls), N, M
c 144           format(7f10.3,1x,2i8)
                endif
            enddo
            enddo
         SIG_CLIP = 3.50
         RES_MIN  = 0.05
         RES_MAX  = 1.00
         call glob_fit6nrDPU_REJ(
     .             x1_l,y1_l,x2_l,y2_l,uu_l,Ls,Lu,           ! added Nu
     .             AG,BG,CG,DG,x1o,y1o,x2o,y2o,
     .             RMS2D,SIG_CLIP,RES_MIN,RES_MAX) ! added sev
         GA =  DG/(AG*DG-BG*CG)
         GB = -BG/(AG*DG-BG*CG)
         GC = -CG/(AG*DG-BG*CG)
         GD =  AG/(AG*DG-BG*CG)
         x2e = CRPIX1_INST_M(M)
         y2e = CRPIX2_INST_M(M)
         if (x2e**2+y2e**2.le.1) then
             x2e = x2o
             y2e = y2o
             endif
         x1e = x1o + GA*(x2e-x2o)
     .             + GB*(y2e-y2o)
         y1e = y1o + GC*(x2e-x2o)
     .             + GD*(y2e-y2o)
         Lm = 0
         do L = 1, Ls
            if (uu_l(L).eq.1) then
                Lm = Lm + 1
                dm_l(Lm) = dm_l(L)
                endif
            enddo
         call rbarsigs(dm_l,Lm,dmbar,dmsig,Lmu,2.5)
         write(70,170) M, Ns_M(M), Ls, Lu,
     .                 AG, BG, CG, DG,
     .                 x1e, y1e, x2e, y2e, RMS2D,
     .                 atan2(BG-CG,AG+DG)*180/3.1415927,
     .                 dmbar, dmsig, Lmu
         write( *,170) M, Ns_M(M), Ls, Lu,
     .                 AG, BG, CG, DG,
     .                 x1e, y1e, x2e, y2e, RMS2D,
     .                 atan2(BG-CG,AG+DG)*180/3.1415927,
     .                 dmbar, dmsig, Lmu
  170    format(1x,i3.3,1x,i5,1x,i5,1x,i5,
     .          2x,4f11.7,
     .          2x,2f10.3,
     .          2x,2f8.1,
     .          2x,f8.4,2x,f8.3,
     .          2x,f8.4,1x,f8.4,1x,i5,1x,a40)
  270    format('#','  M',1x,' Ns_M',1x,' Lmat',1x,' Luse',
     .          2x,'      AG         BG         CG         DG   ',
     .          2x,'    XREF      YREF  ',
     .          2x,'   XLOC    YLOC ',
     .          2x,'  RMS2D ',2x,'  ANGLE ',
     .          2x,'  dZP   '1x,'  sigZP ',1x,'  Lm ')
  370    format('#','---',1x,'-----',1x,'-----',1x,'-----',
     .          2x,' ---------- ---------- ---------- ----------',
     .          2x,' --------- ---------',
     .          2x,' ------- -------',
     .          2x,'--------',2x,'--------',
     .          2x,'--------'1x,'--------',1x,'-----')
  470    format('#',i3.3,1x,i5,1x,i5,1x,i5,
     .          2x,4f11.7,
     .          2x,2f10.3,
     .          2x,2f8.1,
     .          2x,f8.4,2x,f8.3,
     .          2x,f8.4,1x,f8.4,1x,i5)
          AG_MI(M,1) = AG
          BG_MI(M,1) = BG
          CG_MI(M,1) = CG
          DG_MI(M,1) = DG
         x1o_MI(M,1) = x1o
         y1o_MI(M,1) = y1o
         x2o_MI(M,1) = x2o
         y2o_MI(M,1) = y2o
           N_MI(M,1) = Lu
          GA_MI(M,1) = GA
          GB_MI(M,1) = GB
          GC_MI(M,1) = GC
          GD_MI(M,1) = GD
         if (ZP_USE(1:3).eq.'REF') then
             do N = 1, Ns_M(M)
                mt_nm(N,M) = mt_nm(N,M) - dmbar
                enddo
             endif
         do N = 1, Ns_M(M)
            ut_nm(N,M) = x1o + GA_MI(M,1)*(xc_nm(N,M)-x2o)
     .                       + GB_MI(M,1)*(yc_nm(N,M)-y2o)
            vt_nm(N,M) = y1o + GC_MI(M,1)*(xc_nm(N,M)-x2o)
     .                       + GD_MI(M,1)*(yc_nm(N,M)-y2o)
            enddo
         enddo
      write( *,370)
      write( *,270)
      write( *,370)
      write(70,370)
      write(70,270)
      write(70,370)
      close(70)
 
 555  continue
 
      IDIM = NAXIS1_REF
      JDIM = NAXIS2_REF
 
      if (.not.FIND_LST) then
         open(81,file=FILE_LST,status='old')
         Ls = 0
   56    continue
         read(81,'(a200)',end=57) STR200
         if (STR200(1:1).eq.'#') goto 56
         Ls = Ls + 1
         read(STR200,*) xe, ye
         pkx_l(Ls) = int(xe+0.5)
         pky_l(Ls) = int(ye+0.5)
         pkf_l(Ls) = -1
         pkp_l(Ls) = -1
         goto 56
   57    continue
         print*,'FIND_LST:  Ls = ',Ls
         close(81)
         goto 666
         endif
 
 
c---------------------------------------------------
c
c make a starlist that comes from the _flt images...
c
c  1) find bounds of the output frame
c
c---------------------------------------------------
 
 
      print*,' '
      print*,'----------------------------------------------------'
      print*,'ALLOCATE THE PEAK MAP IMAGES:  ',IDIM,' x ',JDIM
      print*,' '
 
 
      allocate(cvg_map(IDIM,JDIM))
      allocate(min_map(IDIM,JDIM))
      allocate(pk1_map(IDIM,JDIM))
      allocate(pk2_map(IDIM,JDIM))
      allocate(pk3_map(IDIM,JDIM))
 
      do i = 1, IDIM
      do j = 1, JDIM
         cvg_map(i,j) = 0
         min_map(i,j) = 0
         pk1_map(i,j) = 0
         pk2_map(i,j) = 0
         pk3_map(i,j) = 0
         enddo
         enddo
 
      open(71,file='hst2collate.02.raw_peaks',status='unknown')
      write(71,371)
      write(71,271)
      write(71,371)
 
      print*,' '
      print*,'----------------------------------------------------'
      print*,'FILL THE PEAK MAPS                                  '
      print*,' '
      do M = 1, Ms
      if (M.le.10.or.M.eq.M/10*10) print*,' ---> M: ',M,' / ',Ms
      do N = 1, Ns_M(M)
         xe = ut_nm(N,M)
         ye = vt_nm(N,M)
         i = int(xe+0.5)
         j = int(ye+0.5)
         if (i.ge.1.and.i.le.IDIM.and.
     .       j.ge.1.and.j.le.JDIM) then
             pk1_map(i,j) = pk1_map(i,j) + 1
             write(71,171) xe, ye,M,N,i,j,
     .                     pk1_map(i,j),
     .                     pk2_map(i,j)
  171        format(1x,f8.2,1x,f8.2,1x,'M',i3.3,1x,i5,1x,i4.4,1x,i4.4,
     .                   1x,i2,1x,i2)
  271        format('#','   umap ',1x,'   vmap ',
     .               1x,'  M',1x,'   N ',1x,'imap',1x,'jmap',
     .               1x,'p1',1x,'p2')
  371        format('#','--------',1x,'--------',
     .               1x,'--',1x,'-----',1x,'----',1x,'----',
     .               1x,'--',1x,'--')
             i = int(xe)
             j = int(ye)
             pk2_map(i  ,j  ) = pk2_map(i  ,j  ) + 1
             pk2_map(i+1,j  ) = pk2_map(i+1,j  ) + 1
             pk2_map(i  ,j+1) = pk2_map(i  ,j+1) + 1
             pk2_map(i+1,j+1) = pk2_map(i+1,j+1) + 1
             endif
         enddo
         enddo
      write(71,371)
      write(71,271)
      write(71,371)
      close(71)
 
      print*,'Make P3 map...'
      do i = 1+1, IDIM-1
      do j = 1+1, JDIM-1
         pk3_map(i,j) = pk1_map(i-1,j+1)+pk1_map(i,j+1)+pk1_map(i+1,j+1)
     .                + pk1_map(i-1,j  )+pk1_map(i,j  )+pk1_map(i+1,j  )
     .                + pk1_map(i-1,j-1)+pk1_map(i,j-1)+pk1_map(i+1,j-1)
         enddo
         enddo
 
 
c
c make coverage map...
c
      print*,'Make cvg map...'
      do M = 1, Ms
      do K = 1, Ks_m(M)
         do P = 1, 4
            xcu = BDRYX_PKM(P,K,M)
            ycu = BDRYY_PKM(P,K,M)
            xrefu = x1o_MI(M,1) + GA_MI(M,1)*(xcu-x2o_MI(M,1))
     .                          + GB_MI(M,1)*(ycu-y2o_MI(M,1))
            yrefu = y1o_MI(M,1) + GC_MI(M,1)*(xcu-x2o_MI(M,1))
     .                          + GD_MI(M,1)*(ycu-y2o_MI(M,1))
            BDRYU_P(P) = xrefu
            BDRYV_P(P) = yrefu
            enddo
         do i = 1, IDIM
         do j = 1, JDIM
            if (inside_poly(i*1.0d0,j*1.0d0,BDRYU_P,BDRYV_P,4))
     .                 cvg_map(i,j) = cvg_map(i,j)+1
            enddo
            enddo
         enddo
         enddo
 
      write(*,*)
      write(*,'('' output finding images '')')
      write(*,*)
 
      write(*,'('' * hst2collate.03.pk1_map.fits: direct peak map'')')
      call writfits_i2('hst2collate.03.pk1_map.fits',pk1_map,IDIM,JDIM)
 
      write(*,'('' * hst2collate.03.pk2_map.fits: 2x2 peaked up'')')
      call writfits_i2('hst2collate.03.pk2_map.fits',pk2_map,IDIM,JDIM)
 
      write(*,'('' * hst2collate.03.pk3_map.fits: 3x3 convolved'')')
      call writfits_i2('hst2collate.03.pk3_map.fits',pk3_map,IDIM,JDIM)
 
      write(*,'('' * hst2collate.03.cvg_map.fits: coverage map'')')
      call writfits_i2('hst2collate.03.cvg_map.fits',cvg_map,IDIM,JDIM)
      write(*,*)
 
 
 
 
c-------------------------------------------------------
c
c now, make a new starlist from the peak map.
c
 
      PKMIN = 2
 
      open(72,file='hst2collate.04.starry_peaks',status='unknown')
      write(72,372)
      write(72,272)
      write(72,372)
 
      print*,' '
      print*,'----------------------------------------------------'
      print*,'MAKE A NEW AVERAGE STARLIST FROM PEAK MAP           '
      print*,' '
      Ls = 0
      do j = 1+1, JDIM-1
      do i = 1+1, IDIM-1
         if (pk2_map(i,j).ge.PKMIN.and.
     .       pk2_map(i,j).ge.0.5*cvg_map(i,j).and.
     .       pk2_map(i,j).ge.pk2_map(i-1,j+1).and.
     .       pk2_map(i,j).ge.pk2_map(i  ,j+1).and.
     .       pk2_map(i,j).ge.pk2_map(i+1,j+1).and.
     .       pk2_map(i,j).ge.pk2_map(i-1,j  ).and.
     .       pk2_map(i,j).ge.pk2_map(i+1,j  ).and.
     .       pk2_map(i,j).ge.pk2_map(i-1,j-1).and.
     .       pk2_map(i,j).ge.pk2_map(i  ,j-1).and.
     .       pk2_map(i,j).ge.pk2_map(i+1,j-1)) then
             Ls = Ls + 1
             pkx_l(Ls) = i
             pky_l(Ls) = j
             pkf_l(Ls) = pk2_map(i,j)
             pkp_l(Ls) = cvg_map(i,j)
             pk2_map(i,j) = pk2_map(i,j)+1               ! this is so that I don't find another star in an adjoining pixel
             write(72,172) i, j, pkf_l(Ls), Ls
  172        format( 1x,i4.4,1x,i4.4,1x,i3,1x,i6.6)
  272        format('#','imap',1x,'jmap',1x,'p2',1x,'    L ')
  372        format('#','----',1x,'----',1x,'--',1x,'------')
             endif
          enddo
          enddo
      write(72,372)
      write(72,272)
      write(72,372)
      close(72)
 
  666 continue
      print*,'    '
      print*,'  Rs: ',Rs
      print*,'  Ls: ',Ls
      print*,'    '
 
      if (Ls.eq.0) stop 'no stars found'
 
 
      allocate(n_lm(Ls,Ms))
 
      do M = 1, Ms
      do L = 1, Ls
         n_lm(L,M) = 0
         enddo
         enddo
 
      print*,' '
      print*,'----------------------------------------------------'
      print*,'COLLATE STARS FROM IMAGES WITH NEW LIST             '
      print*,' '
      do M = 1, Ms
         Lu = 0
         do L = 1, Ls
            dmin = 1.75
            nmin = 0
            do N = 1, Ns_M(M)
               dx = pkx_l(L)-ut_nm(N,M)
               dy = pky_l(L)-vt_nm(N,M)
               dd = sqrt(dx**2+dy**2)
               if (dd.lt.dmin) then
                   dmin = dd
                   nmin = n
                   endif
               enddo
            n_lm(L,M) = nmin
            if (nmin.ne.0) Lu = Lu + 1
            enddo
         write(*,'(5x,''M'',i3.3,1x,i5,''/'',i5)') M, Lu, Ns_M(M)
         enddo
      print*,' '
 
 
 
      open(73,file='hst2collate.05.avgphot1',status='unknown')
      write(73,373)
      write(73,273)
      write(73,373)
 
 
      print*,' '
      print*,'----------------------------------------------------'
      print*,'AVERAGE PHOTOMETRY AND ASTROMETRY NIT = ',NIT
      print*,' '
      write(* ,373)
      write(* ,273)
      write(* ,373)
      mbar_min = 1.00
      do L = 1, Ls
         Os = 0
         do M = 1, Ms
            if (n_lm(L,M).ne.0) then
               N  = n_lm(L,M)
               Os = Os + 1
               x_o(Os) = ut_nm(N,M)
               y_o(Os) = vt_nm(N,M)
               m_o(Os) = mt_nm(N,M)
               q_o(Os) = qr_nm(N,M)
               u_o(Os) = 1
               endif
            enddo
         Om = 0
         Of = pkf_l(L)
         Op = pkp_l(L)
         call barsig2(x_o,y_o,m_o,q_o,Os,
     .                xbar,ybar,mbar,qbar,
     .                xsig,ysig,msig,qsig,
     .                Om,Ou,asig,bsig,u_o,0.035d0)
         xbar_l(L) = xbar
         ybar_l(L) = ybar
         if (INPUT_LIST) then
             xbar_l(L) = xref_r(L)
             ybar_l(L) = yref_r(L)
             endif
         mbar_l(L) = mbar
         xrms_l(L) = xsig
         yrms_l(L) = ysig
         drms_l(L) = sqrt(xsig**2+ysig**2)
         mrms_l(L) = msig
         qbar_l(L) = qbar
         nuse_l(L) = Ou
         SHOW = .false.
         if (mbar.lt.mbar_min) then
             SHOW = .true.
             mbar_min = mbar
             endif
         if (SHOW.or.
     .       L.le.010.or.
     .       L.lt.0100.and.(L.eq.L/10*10).or.
     .       L.eq.L/100*100)
     .       write(*,173) xbar, ybar, mbar,
     .                    xsig, ysig, msig,
     .                    qbar, Op, Of, Os, Ou, Om, L
         write(73,173) xbar, ybar, mbar,
     .                 xsig, ysig, msig,
     .                 qbar, Op, Of, Os, Ou, Om, L
  173    format(' ',f9.3,1x,f9.3,1x,f8.4,
     .          3x,f7.4,1x,f7.4,1x,f7.4,
     .          3x,f6.3,
     .          3x,i3,1x,i3,1x,i3,1x,i3,1x,i3,
     .          3x,'L',i6.6)
  273    format('#','   xbar  ',1x,'   ybar  ',1x,'  mbar  ',
     .          3x,'  xrms ',1x,'  yrms ',1x,'  mrms ',
     .          3x,' qbar ',
     .          3x,' Op',1x,' Of',1x,' Os',1x,' Ou',1x,' Om',
     .          3x,' ','STARNO')
  373    format('#','---------',1x,'---------',1x,'--------',
     .          3x,'-------',1x,'-------',1x,'-------',
     .          3x,'------',
     .          3x,'---',1x,'---',1x,'---',1x,'---',1x,'---',
     .          3x,'-','------')
         enddo
      write( *,373)
      write( *,273)
      write( *,373)
      write(73,373)
      write(73,273)
      write(73,373)
      close(73)
 
      if (FIND_MAT) then
          print*,' '
          print*,'----------------------------------------------------'
          print*,'RECOMPUTE THE TRANSFORMATIONS                       '
          print*,' '
          open(74,file='hst2collate.06.trans2',status='unknown')
          write(74,370)
          write(74,270)
          write(74,370)
          write( *,370)
          write( *,270)
          write( *,370)
          endif
 
 
      NIT = 0
    2 continue
      NIT = NIT + 1
      if (NITs.eq.0) NIT = 0
 
      do L = 1, Ls
         utrans_l(L) = 0
         if (drms_L(L).le.0.07) utrans_l(L) = 1
         enddo
      do M = 1, Ms
         if (FIND_MAT) then
             Lt = 0
             do L = 1, Ls
                uu_l(L) = 0
                if (n_lm(L,M).ne.0) then
                    uu_l(L) = utrans_l(L)
                    x2_l(L) = xc_nm(n_lm(L,M),M)
                    y2_l(L) = yc_nm(n_lm(L,M),M)
                    m2_l(L) = mt_nm(n_lm(L,M),M)
                    dm_l(L) = mt_nm(n_lm(L,M),M)-mbar_l(L)
                    Lt = Lt + 1
                    endif
                enddo
             SIG_CLIP = 3.50
             RES_MIN  = 0.05
             RES_MAX  = 1.00
             call glob_fit6nrDPU_REJ(
     .                 xbar_l,ybar_l,x2_l,y2_l,uu_l,Ls,Lu,           ! added Nu
     .                 AG,BG,CG,DG,x1o,y1o,x2o,y2o,
     .                 RMS2D,SIG_CLIP,RES_MIN,RES_MAX) ! added sev
             AG_MI(M,NIT+1) = AG
             BG_MI(M,NIT+1) = BG
             CG_MI(M,NIT+1) = CG
             DG_MI(M,NIT+1) = DG
            x1o_MI(M,NIT+1) = x1o
            y1o_MI(M,NIT+1) = y1o
            x2o_MI(M,NIT+1) = x2o
            y2o_MI(M,NIT+1) = y2o
              N_MI(M,NIT+1) = Lu
             GA_MI(M,NIT+1) =  DG/(AG*DG-BG*CG)
             GB_MI(M,NIT+1) = -BG/(AG*DG-BG*CG)
             GC_MI(M,NIT+1) = -CG/(AG*DG-BG*CG)
             GD_MI(M,NIT+1) =  AG/(AG*DG-BG*CG)
             x2e = CRPIX1_INST_M(M)
             y2e = CRPIX2_INST_M(M)
             if (x2e**2+y2e**2.le.1) then
                 x2e = x2o
                 y2e = y2o
                 endif
             x1e = x1o + GA_MI(M,NIT+1)*(x2e-x2o)
     .                 + GB_MI(M,NIT+1)*(y2e-y2o)
             y1e = y1o + GC_MI(M,NIT+1)*(x2e-x2o)
     .                 + GD_MI(M,NIT+1)*(y2e-y2o)
             Lm = 0
             do L = 1, Ls
                if (uu_l(L).eq.1) then
                    Lm = Lm + 1
                    dm_l(Lm) = dm_l(L)
                    endif
                enddo
             call rbarsigs(dm_l,Lm,dmbar,dmsig,Lmu,2.5)
             write(74,170) M, Ns_M(M), Lt, Lu,
     .                     AG, BG, CG, DG,
     .                     x1e, y1e, x2e, y2e, RMS2D,
     .                     atan2(BG-CG,AG+DG)*180/3.1415927,
     .                     dmbar, dmsig,Lmu
             write( *,170) M, Ns_M(M), Lt, Lu,
     .                     AG, BG, CG, DG,
     .                     x1e, y1e, x2e, y2e, RMS2D,
     .                     atan2(BG-CG,AG+DG)*180/3.1415927,
     .                     dmbar, dmsig,Lmu
             do N = 1, Ns_M(M)
                ut_nm(N,M) = x1o + GA_MI(M,NIT+1)*(xc_nm(N,M)-x2o)
     .                           + GB_MI(M,NIT+1)*(yc_nm(N,M)-y2o)
                vt_nm(N,M) = y1o + GC_MI(M,NIT+1)*(xc_nm(N,M)-x2o)
     .                           + GD_MI(M,NIT+1)*(yc_nm(N,M)-y2o)
                mt_nm(N,M) = mt_nm(N,M) - dmbar
                enddo
 
             MATFILE = FILE_M(M)(1:9) // '_mat.UVuvWw'
             open(91,file=MATFILE,status='unknown')
             write(91,'(''#'')')
             write(91,'(''# Generated by hst2collate ---'')')
             write(91,'(''#'')')
             write(91,'(''#'')')
             write(91,'(''#   NAXIS1: '',i6)')    NAXIS1_REF
             write(91,'(''#   NAXIS2: '',i6)')    NAXIS2_REF
             write(91,'(''#     rDAT: '',f15.8)')  rDAT
             write(91,'(''#'')')
             write(91,'(''#   CRPIX1: '',f15.8)') CRPIX1_REF
             write(91,'(''#   CRPIX2: '',f15.8)') CRPIX2_REF
             write(91,'(''#   CRVAL1: '',f15.8)') CRVAL1_REF
             write(91,'(''#   CRVAL2: '',f15.8)') CRVAL2_REF
             write(91,'(''#'')')
             write(91,'(''#    CD1_1: '',f15.12)') CD1_1_REF
             write(91,'(''#    CD1_2: '',f15.12)') CD1_2_REF
             write(91,'(''#    CD2_1: '',f15.12)') CD2_1_REF
             write(91,'(''#    CD2_2: '',f15.12)') CD2_2_REF
             write(91,'(''#'')')
             write(91,370)
             write(91,270)
             write(91,370)
             write(91,470) M, Ns_M(M), Lt, Lu,
     .                     AG, BG, CG, DG,
     .                     x1e, y1e, x2e, y2e, RMS2D,
     .                     atan2(BG-CG,AG+DG)*180/3.1415927,
     .                     dmbar, dmsig,Lmu
             write(91,'(''#'')')
             write(91,'(''#   STDGDC: '',a80)') STDGDC_M(M)(1:80)
             write(91,'(''#'')')
             write(91,391)
             write(91,291)
             write(91,391)
             do L = 1, Ls
                if (uu_l(L).eq.1) then
                    N = n_lm(L,M)
                    x1e = x1o + GA_MI(M,NIT+1)*(xc_nm(N,M)-x2o)
     .                        + GB_MI(M,NIT+1)*(yc_nm(N,M)-y2o)
                    y1e = y1o + GC_MI(M,NIT+1)*(xc_nm(N,M)-x2o)
     .                        + GD_MI(M,NIT+1)*(yc_nm(N,M)-y2o)
                    x2e = x2o + AG_MI(M,NIT+1)*(xbar_l(L) -x1o)
     .                        + BG_MI(M,NIT+1)*(ybar_l(L) -y1o)
                    y2e = y2o + CG_MI(M,NIT+1)*(xbar_l(L) -x1o)
     .                        + DG_MI(M,NIT+1)*(ybar_l(L) -y1o)
                    write(91,191) xbar_l(L), ybar_l(L),
     .                            xc_nm(N,M), yc_nm(N,M),
     .                            mbar_l(L),  mc_nm(N,M),
     .                            L, N,
     .                            xbar_l(L)-x1e,
     .                            ybar_l(L)-y1e,
     .                            xc_nm(N,M)-x2e,
     .                            yc_nm(N,M)-y2e
  191               format(1x,f9.3,1x,f9.3,
     .                     3x,f9.3,1x,f9.3,
     .                     3x,f8.4,1x,f8.4,
     .                     3x,'L',i5.5,1x,'N',i5.5,
     .                     3x,f8.4,1x,f8.4,
     .                     3x,f8.4,1x,f8.4)
  291               format('#','    UREF ',1x,'    VREF ',
     .                     3x,'    XC   ',1x,'    YC   ',
     .                     3x,'  MREF  ' ,1x,'  MC    ',
     .                     3x,'REF ID',1x,'  N_M ',
     .                     3x,' d(UREF)',1x,' d(VREF)',
     .                     3x,'  d(XC) ',1x,'  d(XC) ')
  391               format('#','---------',1x,'---------',
     .                     3x,'---------',1x,'---------',
     .                     3x,'--------',1x,'--------',
     .                     3x,'------',1x,'------',
     .                     3x,'--------',1x,'--------',
     .                     3x,'--------',1x,'--------')
                    endif
                enddo
             write(91,391)
             write(91,291)
             write(91,391)
             close(91)
             endif
 
c
c LNKFILE output
c
         LNKFILE = FILE_M(M)(1:9) // '_lnk.UVWuvXY_'
     .          // COLS_M(M)(1:NCOLS_M(M))
         open(92,file=LNKFILE,status='unknown')
         write(92,'(''#'')')
         write(92,'(''# Generated by hst2collate ---'')')
         write(92,'(''#'')')
         write(92,370)
         write(92,270)
         write(92,370)
            x1o_MI(M,NIT+1) = x1o
            y1o_MI(M,NIT+1) = y1o
            x2o_MI(M,NIT+1) = x2o
            y2o_MI(M,NIT+1) = y2o
         write(92,470) M, Ns_M(M), Lt, Lu,
     .                  AG_MI(M,NIT+1),  BG_MI(M,NIT+1),
     .                  CG_MI(M,NIT+1),  DG_MI(M,NIT+1),
     .                 x1o_MI(M,NIT+1), y1o_MI(M,NIT+1),
     .                 x2o_MI(M,NIT+1), y2o_MI(M,NIT+1),
     .                 RMS2D,
     .                 atan2(BG-CG,AG+DG)*180/3.1415927,
     .                 ZP_MI(M,NIT+1), eZP_MI(M,NIT+1)
         write(92,'(''#'')')
         write(92,'(''# FILE_M: '',200a)')
     .                  FILE_M(M)(1:strlen(FILE_M(M)))
         write(92,'(''# STDGDC: '',200a)')
     .                  STDGDC_M(M)(1:strlen(STDGDC_M(M)))
         write(92,'(''#'')')
         write(92,'(''#  COL01: UMAP - mapped ref frame pos, U'')')
         write(92,'(''#  COL02: VMAP - mapped ref frame pos, V'')')
         write(92,'(''#  COL03: MMAP - mapped ref frame mag, W'')')
         write(92,'(''#  COL04:  XCE - mapped dist cor pos, u'')')
         write(92,'(''#  COL05:  YCE - mapped dist cor pos, v'')')
         write(92,'(''#  COL06:  XRE - mapped raw chip pos, X'')')
         write(92,'(''#  COL07:  YRE - mapped raw chip pos, Y'')')
         write(92,'(''#  COL08:    L - star number in list   '')')
         write(92,'(''#  COL09+        line from hst1pass file'')')
         write(92,'(''#'')')
         write(92,392)
         write(92,292)
         write(92,392)
         write(s200h,'(200('' ''))')
         do i = 1, NCOLS_M(M)
            s200h(1+(i-1)*2:1+(i-1)*2) = '0'
            enddo
         do L = 1, Ls
            N = n_lm(L,M)
            xce = x2o + AG_MI(M,NIT+1)*(xbar_l(L)-x1o)
     .                + BG_MI(M,NIT+1)*(ybar_l(L)-y1o)
            yce = y2o + CG_MI(M,NIT+1)*(xbar_l(L)-x1o)
     .                + DG_MI(M,NIT+1)*(ybar_l(L)-y1o)
            if (xce.gt.9999.0) xce = 9999.0
            if (yce.gt.9999.0) yce = 9999.0
            if (xce.lt.-999.0) xce = -999.0
            if (yce.lt.-999.0) yce = -999.0
            call xcyc2xryr_stdgc(xce,yce,xre,yre,STDGDC_M(M))
            if (xre.gt.9999.0) xre = 9999.0
            if (yre.gt.9999.0) yre = 9999.0
            if (xre.lt.-999.0) xre = -999.0
            if (yre.lt.-999.0) yre = -999.0
            ute = 0.
            vte = 0.
            mte = 0.
            s200o = s200h
            if (N.ne.0) then
                ute = ut_nm(N,M)
                vte = vt_nm(N,M)
                mte = mt_nm(N,M)
                s200o = s200_nm(N,M)
                endif
            iend = 1
            do i = 1, 200
               if (s200o(i:i).ne.' ') iend = i
               enddo
            write(92,192) ute, vte, mte,
     .                    xce, yce,
     .                    xre, yre,
     .                    L,
     .                    s200o(1:iend)
  192       format(1x,f9.4,1x,f9.4,1x,f8.4,
     .             3x,f8.3,1x,f8.3,
     .             3x,f8.3,1x,f8.3,
     .             3x,'N',i6.6,
     .             3x,200a)
  292       format('#','  UMAP    ',1x,'  VMAP   ',1x,'  MMAP ',
     .              3x,'  XCE   ',1x,'  YCE   '
     .              3x,'  XRE   ',1x,'  YRE   '
     .              3x,' LIST# ',
     .              3x,'LINE FROM HST1PASS OUTPUT FILE')
  392       format('#','---------',1x,'---------',1x,'--------',
     .              3x,'--------',1x,'--------'
     .              3x,'--------',1x,'--------'
     .              3x,'-------',
     .              3x,'------------------------------------')
            enddo
         enddo!M
 
      if (FIND_MAT) then
          write( *,370)
          write( *,270)
          write( *,370)
          print*,' '
          write(74,370)
          write(74,270)
          write(74,370)
          close(74)
          endif
 
      if (NITs.eq.0) goto 3
 
 
 
      open(75,file='hst2collate.07.avgphot2',status='unknown')
      write(75,'(''# '')')
      write(75,'(''#   NAXIS1: '',i8   )') NAXIS1_REF
      write(75,'(''#   NAXIS2: '',i8   )') NAXIS2_REF
      write(75,'(''#   CRPIX1: '',f16.7)') CRPIX1_REF
      write(75,'(''#   CRPIX2: '',f16.7)') CRPIX2_REF
      write(75,'(''#   CRVAL1: '',f16.7)') CRVAL1_REF
      write(75,'(''#   CRVAL2: '',f16.7)') CRVAL2_REF
      write(75,'(''#    CD1_1: '',f16.7)')  CD1_1_REF
      write(75,'(''#    CD1_2: '',f16.7)')  CD1_2_REF
      write(75,'(''#    CD2_1: '',f16.7)')  CD2_1_REF
      write(75,'(''#    CD2_2: '',f16.7)')  CD2_2_REF
      write(75,'(''#    CD2_2: '',f16.7)')  CD2_2_REF
      write(75,'(''#     rDAT: '',f16.7)')  rDAT
      write(75,'(''# '')')
      write(75,373)
      write(75,273)
      write(75,373)
      write(*,373)
      write(*,273)
      write(*,373)
 
      print*,' '
      print*,'----------------------------------------------------'
      print*,'RE-AVERAGE PHOTOMETRY AND ASTROMETRY                '
      print*,' '
      mbar_min = 1.00
      do L = 1, Ls
         Os = 0
         do M = 1, Ms
            if (n_lm(L,M).ne.0) then
               N  = n_lm(L,M)
               Os = Os + 1
               x_o(Os) = ut_nm(N,M)
               y_o(Os) = vt_nm(N,M)
               m_o(Os) = mt_nm(N,M)
               q_o(Os) = qr_nm(N,M)
               u_o(Os) = 1
               endif
            enddo
         Om = 0
         Of = pkf_l(L)
         Op = pkp_l(L)
         call barsig2(x_o,y_o,m_o,q_o,Os,
     .                xbar,ybar,mbar,qbar,
     .                xsig,ysig,msig,qsig,
     .                Om,Ou,asig,bsig,u_o,0.035d0)
         xbar_l(L) = xbar
         ybar_l(L) = ybar
         mbar_l(L) = mbar
         xrms_l(L) = xsig
         yrms_l(L) = ysig
         drms_l(L) = sqrt(xsig**2+ysig**2)
         mrms_l(L) = msig
         qbar_l(L) = qbar
         nuse_l(L) = Ou
         SHOW = .false.
         if (mbar.lt.mbar_min) then
             SHOW = .true.
             mbar_min = mbar
             endif
         if (SHOW.or.
     .       L.le.010.or.
     .       L.lt.0100.and.(L.eq.L/10*10).or.
     .       L.eq.L/100*100)
     .       write(*,173) xbar, ybar, mbar,
     .                    xsig, ysig, msig,
     .                    qbar, Op, Of, Os, Ou, Om, L
         write(75,173) xbar, ybar, mbar,
     .                 xsig, ysig, msig,
     .                 qbar, Op, Of, Os, Ou, Om, L
         enddo
      write( *,373)
      write( *,273)
      write( *,373)
      write(75,373)
      write(75,273)
      write(75,373)
      close(75)
 
      print*,'          '
      print*,'    NIT : ',NIT
      print*,'    NITs: ',NITs
      print*,'          '
 
      if (NIT.lt.NITs) goto 2
 
    3 continue  ! done with iterations...
      print*,'done!'
 
      STOP
      end
 
 
 
                                              
c*********************************************
c****                                         
c**** #include "SUBROUTINES/barsig2.f"
c****                                         
c*********************************************
                                              
 
 
c#include "/user/jayander/FORTRAN/ROUTINES/SORT/cosrt2ri.f"
c#include "/user/jayander/FORTRAN/ROUTINES/SORT/cosrt2ii.f"
c#include "/user/jayander/FORTRAN/ROUTINES/SORT/cosort.f"
c#include "/user/jayander/FORTRAN/ROUTINES/SORT/ord_copy.f"
c#include "/user/jayander/FORTRAN/ROUTINES/SORT/rbubble.f"
 
 
c-----------------------------------------------------------
c
c this routine will take the x and y arrays
c of peaks (one in each image...) and will find
c the N (1st test: 5) closest and will analyze
c these 5 to produce an average/sigma.  There will
c of course be a bias in this result, but what I
c am doing here is essentially a MODE in 2 dimensions.
c
c I'll carry along the mag info so that I can average
c that for the same reported peaks
c
c
c NTOT : tot number of measts
c NUSE : the number included in the initial "mode"
c        finding (the closest NUSE...)
c NUSED: the actual number included (>=NUSE); I
c        will include those further and further
c        out provided they are consistent with
c        the inner points
c
 
      subroutine barsig2(xl,yl,ml,ql,NTOT,
     .                   xbar,ybar,mbar,qbar,
     .                   xsig,ysig,msig,qsig,
     .                   NMIN,NUSED,
     .                   asig,bsig,wl,DKEEP)
      implicit none
 
      integer NTOT
      real*8 xl(NTOT)
      real*8 yl(NTOT)
      real*8 ml(NTOT)
      real*8 ql(NTOT)
      real*8 xbar, ybar, mbar, qbar
      real*8 xsig, ysig, msig, qsig
      integer NMIN, NUSED
      real*8 asig, bsig
      integer wl(NTOT)
      integer n
      real*8 xsum, ysum, msum, qsum
      real*8 asum, bsum
      real*8 dsig
      real*8  dmax,dmaxx
      integer nmax
      real*8  DKEEP
 
      real*8 xbarx, ybarx
      real*8 xsigx, ysigx, dsigx
 
      real dx, dy, dd
      integer wsum, wsumx
 
 
      xbar = 0.
      ybar = 0.
      mbar = 0.
      qbar = 0.
      xsig = 1.
      ysig = 1.
      msig = 1.
      qsig = 1.
      NUSED = 0
      if (NTOT.lt.1) return
      if (NMIN.gt.NTOT) then
          print*,'NMIN.gt.NTOT'
          return
          endif
 
      if (NTOT.eq.1) then
         xbar = xl(1)
         ybar = yl(1)
         mbar = ml(1)
         qbar = ql(1)
         NUSED = 1
         return
         endif
 
      if (NTOT.eq.2) then
         xbar = (xl(1)+xl(2))/2
         ybar = (yl(1)+yl(2))/2
         mbar = (ml(1)+ml(2))/2
         qbar = (ql(1)+ql(2))/2
         xsig = abs(xl(1)-xl(2))
         ysig = abs(yl(1)-yl(2))
         msig = abs(ml(1)-ml(2))
         qsig = abs(ql(1)-ql(2))
         NUSED = 2
         return
         endif
 
      if (NMIN.gt.NTOT) stop 'BIG PROBLEM'
      if (NTOT.gt.1000) stop 'QQUIT'
 
      do n = 1, NTOT
         wl(n) = 1
         enddo
 
      xbar = 000
      ybar = 000
      dsig = 9e9
      dmax = 9e9
  1   continue
 
      xsum = 0.
      ysum = 0.
      msum = 0.
      qsum = 0.
      wsum = 0.
      do n = 1, NTOT
         xsum = xsum + wl(n)*xl(n)
         ysum = ysum + wl(n)*yl(n)
         msum = msum + wl(n)*ml(n)
         qsum = qsum + wl(n)*ql(n)
         wsum = wsum + wl(n)
         enddo
      if (wsum.lt.2) then
         xbar = 0.
         ybar = 0.
         mbar = 0.
         xsig = 1.
         ysig = 1.
         msig = 1.
         qsig = 1.
         return
         endif
      xbar = xsum/wsum
      ybar = ysum/wsum
      mbar = msum/wsum
      qbar = qsum/wsum
 
      xsum = 0.
      ysum = 0.
      msum = 0.
      qsum = 0.
      asum = 0.
      bsum = 0.
      dmax = 0.
      nmax = 0
      do n = 1, NTOT
         xsum = xsum + wl(n)*(xl(n)-xbar)**2
         ysum = ysum + wl(n)*(yl(n)-ybar)**2
         msum = msum + wl(n)*(ml(n)-mbar)**2
         qsum = qsum + wl(n)*(ql(n)-qbar)**2
         asum = asum + wl(n)*((xl(n)-xbar+yl(n)-ybar)/1.41421)**2
         bsum = bsum + wl(n)*((xl(n)-xbar-yl(n)+ybar)/1.41421)**2
         dx = xl(n)-xbar
         dy = yl(n)-ybar
         dd = sqrt(dx**2+dy**2)
         if (dd.gt.dmax.and.wl(n).eq.1) then
            dmax = dd
            nmax = n
            endif
         enddo
      xsig = sqrt(xsum/(wsum-1))
      ysig = sqrt(ysum/(wsum-1))
      asig = sqrt(asum/(wsum-1))
      bsig = sqrt(bsum/(wsum-1))
      msig = sqrt(msum/(wsum-1))
      qsig = sqrt(qsum/(wsum-1))
      dsig = sqrt(xsig**2+ysig**2)
 
      xsum = 0.
      ysum = 0.
      wsumx= 0.
      do n = 1, NTOT
         if (n.ne.nmax) then
            xsum = xsum + wl(n)*xl(n)
            ysum = ysum + wl(n)*yl(n)
            wsumx= wsumx+ wl(n)
            endif
         enddo
      xbarx = xsum/wsumx
      ybarx = ysum/wsumx
      dmaxx = sqrt((xl(NMAX)-xbarx)**2+(yl(NMAX)-ybarx)**2)
      xsum = 0.
      ysum = 0.
      do n = 1, NTOT
         if (n.ne.nmax) then
            xsum = xsum + wl(n)*(xl(n)-xbarx)**2
            ysum = ysum + wl(n)*(yl(n)-ybarx)**2
            endif
         enddo
      xsigx = sqrt(xsum/(wsumx-1))
      ysigx = sqrt(ysum/(wsumx-1))
      dsigx = sqrt(xsigx**2+ysigx**2)
 
 
c 7/20/2007: the 4*2dsig was waaaay too generous
c
c     if (dmaxx.gt.4.00*dsigx.and.
      if (dmaxx.gt.DKEEP.and.
     .    dmaxx.gt.2.00*dsigx.and.
     .    wsum.gt.NMIN.and.
     .    nmax.gt.0) then
c        print*,'REJECT! '
         wl(nmax) = 0
         goto 1
         endif
 
      NUSED = wsum
 
 
c     if (.not.(xbar.gt.  0.0)) xbar = 0.000
c     if (.not.(ybar.gt.  0.0)) ybar = 0.000
      if (.not.(mbar.gt.-99.0)) mbar = 0.000
      if (.not.(qbar.gt. -0.0)) qbar = 9.999
      if (.not.(xsig.lt.  9.0)) xsig = 9.999
      if (.not.(ysig.lt.  9.0)) ysig = 9.999
      if (.not.(msig.lt.  9.0)) msig = 9.999
      if (.not.(qsig.lt.  9.0)) qsig = 9.999
      if (.not.(asig.lt.  9.0)) asig = 9.999
      if (.not.(bsig.lt.  9.0)) bsig = 9.999
 
 
      return
      end
 
 
 
 
c      subroutine bail(STRING)
c      implicit none
c      character*80 STRING
c      integer i
c      print*,' '
c      print*,'BAIL: '
c      print*,' '
c      i = 0
c  1   continue
c      i = i + 1
c      if (string(i:i).ne.'\n'.and.i.le.79) goto 1
c      print*,' '
c      print*,STRING(1:i)
c      print*,' '
c      print*,' '
c      stop
c      end
c
 
 
c----------------------------------------------------
c
c This routine will take a list of NTOT values and
c will find the average and the std deviation of
c the points... right now, it takes the avg of the
c absolute value, not the std devn.  This is more
c robust.
c
      subroutine barsig(xlist,NTOT,bar,sig,NUSE)
      implicit none
 
      integer NTOT
      real*8 xlist(NTOT)
      real*8 bar
      real*8 sig
      integer NUSE
 
      integer n
      real*8    bsum, ssum
      integer nsum
      integer NIT
 
 
      bar = 0.e0
      sig = 9e9
      do NIT = 1, 8
         bsum = 0.
         ssum = 0.
         nsum = 0.
         do n = 1, NTOT
            if (abs(xlist(n)-bar).le.3*sig) then
               bsum = bsum + xlist(n)
               ssum = ssum + abs(xlist(n)-bar)
               nsum = nsum + 1
               endif
            enddo
         if (nsum.gt.0) bar = bsum / nsum
         if (nsum.gt.1) sig = ssum/(nsum-1)
         enddo
      NUSE = nsum
      if (nsum.le.1) sig = 9.999
 
      return
      end
 
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/PARLOR/islinux.f"
c****                                         
c*********************************************
                                              
 
c      program q
c      implicit none
c
c      logical islinux, islinux_var
c
c      islinux_var = islinux()
c
c      islinux_var = islinux()
c
c      stop
c      end
 
 
      logical function islinux()
      implicit none
 
      logical islinux_save
 
      integer i
      data i/0/
      data islinux_save/.true./
 
      common /islinux_/i,islinux_save
 
      byte    b(2)
      equivalence(i,b)
 
      if (i.eq.1) goto 1
 
      i = 1
      islinux_save = .false.
      if (b(1).eq.1) islinux_save = .true.
 
    1 islinux = islinux_save
      return
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/FITSIO/GEN/writfits_b1.f"
c****                                         
c*********************************************
                                              
 
 
      subroutine writfits_b1(FILE,pix,PXDIMX,PXDIMY)
      implicit none
 
      character*(*) FILE
      integer PXDIMX,PXDIMY
      byte pix(PXDIMX,PXDIMY)
 
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffc,buffb)
 
      integer ifirst, i1, i2
 
      integer np1, np2, npt
 
      character*70 HDR(25)
      common/HDR/HDR
 
      character*80 FILEU
 
      FILEU = FILE
      do i = 75,2,-1
         if (FILE(i:i+4).eq.'.fits') FILEU = FILE(1:i+4)
         enddo
 
 
 
c     write(*,'(''ENTER   writfits_b1: '',80a)') FILEU
 
      open(10,file=FILEU,
     .     status='unknown',
     .     err   =900,
     .     recl  =2880,
     .     form  ='UNFORMATTED',
     .     access='DIRECT')
 
 
      i = 1
      write(buffc( 0*80+1: 1*80),'(''SIMPLE  = T    '')')
      write(buffc( 1*80+1: 2*80),'(''BITPIX  =  8  '')')
      write(buffc( 2*80+1: 3*80),'(''NAXIS   = '',i12)') 2
      write(buffc( 3*80+1: 4*80),'(''NAXIS1  = '',i12)') PXDIMX
      write(buffc( 4*80+1: 5*80),'(''NAXIS2  = '',i12)') PXDIMY
      write(buffc( 5*80+1: 6*80),'(''DATATYPE= '',9a)')
     .                          ' ''BYTE    '''
      write(buffc( 6*80+1: 7*80),'(''DATE    = '',11a)')
     .                          ' ''28/01/00'''
      write(buffc( 7*80+1: 8*80),'(''BSCALE  = '',i12)') 00001
      write(buffc( 8*80+1: 9*80),'(''BZERO   = '',i12)') 00000
      write(buffc(09*80+1:10*80),'(''CRPIX1  = '',a20)') HDR(01)
      write(buffc(10*80+1:11*80),'(''CRPIX2  = '',a20)') HDR(02)
      write(buffc(11*80+1:12*80),'(''CRVAL1  = '',a20)') HDR(03)
      write(buffc(12*80+1:13*80),'(''CRVAL2  = '',a20)') HDR(04)
      write(buffc(13*80+1:14*80),'(''CTYPE1  = '',a20)') HDR(05)
      write(buffc(14*80+1:15*80),'(''CTYPE2  = '',a20)') HDR(06)
      write(buffc(15*80+1:16*80),'(''CD1_1   = '',a20)') HDR(07)
      write(buffc(16*80+1:17*80),'(''CD1_2   = '',a20)') HDR(08)
      write(buffc(17*80+1:18*80),'(''CD2_1   = '',a20)') HDR(09)
      write(buffc(18*80+1:19*80),'(''CD2_2   = '',a20)') HDR(10)
      write(buffc(19*80+1:20*80),'(''ORIENTAT= '',a20)') HDR(11)
      write(buffc(20*80+1:21*80),'(''PA_APER = '',a20)') HDR(12)
      write(buffc(21*80+1:22*80),'(''PA_V3   = '',a20)') HDR(13)
      write(buffc(22*80+1:23*80),'(''DATE-OBS= '',a20)') HDR(14)
      write(buffc(23*80+1:24*80),'(''TIME-OBS= '',a20)') HDR(15)
      write(buffc(24*80+1:25*80),'(''EXPTIME = '',a20)') HDR(16)
      write(buffc(25*80+1:26*80),'(''ROOTNAME= '',a20)') HDR(17)
      write(buffc(26*80+1:27*80),'(''TARGNAME= '',a20)') HDR(18)
      write(buffc(27*80+1:28*80),'(''RA_TARG = '',a20)') HDR(19)
      write(buffc(28*80+1:29*80),'(''DEC_TARG= '',a20)') HDR(20)
      write(buffc(29*80+1:30*80),'(''PROPOSID= '',a20)') HDR(21)
      write(buffc(30*80+1:31*80),'(''FILTER1 = '',a20)') HDR(22)
      write(buffc(31*80+1:32*80),'(''FILTER2 = '',a20)') HDR(23)
      write(buffc(32*80+1:33*80),'(''VAFACTOR= '',a20)') HDR(24)
      write(buffc(33*80+1:34*80),'(''COMMENT   '')')
      write(buffc(34*80+1:35*80),'(''COMMENT   '')')
      write(buffc(35*80+1:36*80),'(''END       '')')
 
c     print*,buffc
      write(10,rec=i,iostat=ios) buffc
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      nbper  =   PXDIMX*PXDIMY
      npt    =   PXDIMX*PXDIMY
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
 
c      print*,'nbyte1: ',nbyte1
c      print*,'nbyte2: ',nbyte1
c      print*,'    i1: ',i1
c      print*,'    i2: ',i2
      do i = i1, i2, 1
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1) + 1
         np2 = (nbyteE-nbyte1)+ 1
        !if (i.lt.0010) print*,'i: ',i,i1,i2,np1,np2,nbyte0
        !if (i.gt.2900) print*,'i: ',i,i1,i2,np1,np2,nbyte0
         call pix2buff_b1(pix,buffb,np1,npt)
         write(10,rec=i,iostat=ios) buffc
         enddo
 
      close(10)
 
      return
 
 900  continue
      print*,'WRITFITS.f ERROR'
      stop
 
      end
 
      subroutine pix2buff_b1(pix,buffb,n1,nt)
      implicit none
      byte pix(1)
      byte buffb(2880)
 
      integer   n1,nt
 
      integer i
      integer npu
 
      byte ipval
 
      do i = 1, 2880
         npu = (1+(n1-1)*2)+i-1
         npu = (1+(n1-1)*1)+i-1
         if (npu.ge.1.and.npu.le.nt) then
            ipval = pix(npu)
            buffb(i) = ipval
            endif
         enddo
 
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/FITSIO/GEN/writfits_i2.f"
c****                                         
c*********************************************
                                              
 
      subroutine writfits_i2(FILE,pix,PXDIMX,PXDIMY)
      implicit none
 
      integer PXDIMX,PXDIMY
      character*(*) FILE
      integer*2 pix(PXDIMX,PXDIMY)
 
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffb,buffc)
 
      integer ifirst, i1, i2
 
      integer np1, np2, npt
 
      character*70 HDR(25)
      common/HDR/HDR
 
      character*80 FILEU
 
      FILEU = FILE
      do i = 75,2,-1
         if (FILE(i:i+4).eq.'.fits') FILEU = FILE(1:i+4)
         enddo
 
      !write(*,'(''ENTER   writfits_i2: '',80a)') FILEU
 
      open(10,file=FILEU,
     .     status='unknown',
     .     err   =900,
     .     recl  =2880,
     .     form  ='UNFORMATTED',
     .     access='DIRECT')
 
 
      i = 1
      write(buffc( 0*80+1: 1*80),'(''SIMPLE  = T    '')')
      write(buffc( 1*80+1: 2*80),'(''BITPIX  =  16  '')')
      write(buffc( 2*80+1: 3*80),'(''NAXIS   = '',i12)') 2
      write(buffc( 3*80+1: 4*80),'(''NAXIS1  = '',i12)') PXDIMX
      write(buffc( 4*80+1: 5*80),'(''NAXIS2  = '',i12)') PXDIMY
      write(buffc( 5*80+1: 6*80),'(''DATATYPE= '',9a)')
     .                          ' ''INTEGER*4 '''
      write(buffc( 6*80+1: 7*80),'(''DATE    = '',11a)')
     .                          ' ''28/01/00'''
      write(buffc( 7*80+1: 8*80),'(''BSCALE  = '',i12)') 00001
      write(buffc( 8*80+1: 9*80),'(''BZERO   = '',i12)') 00000
      write(buffc(09*80+1:10*80),'(''CRPIX1  = '',a20)') HDR(01)
      write(buffc(10*80+1:11*80),'(''CRPIX2  = '',a20)') HDR(02)
      write(buffc(11*80+1:12*80),'(''CRVAL1  = '',a20)') HDR(03)
      write(buffc(12*80+1:13*80),'(''CRVAL2  = '',a20)') HDR(04)
      write(buffc(13*80+1:14*80),'(''CTYPE1  = '',a20)') HDR(05)
      write(buffc(14*80+1:15*80),'(''CTYPE2  = '',a20)') HDR(06)
      write(buffc(15*80+1:16*80),'(''CD1_1   = '',a20)') HDR(07)
      write(buffc(16*80+1:17*80),'(''CD1_2   = '',a20)') HDR(08)
      write(buffc(17*80+1:18*80),'(''CD2_1   = '',a20)') HDR(09)
      write(buffc(18*80+1:19*80),'(''CD2_2   = '',a20)') HDR(10)
      write(buffc(19*80+1:20*80),'(''ORIENTAT= '',a20)') HDR(11)
      write(buffc(20*80+1:21*80),'(''PA_APER = '',a20)') HDR(12)
      write(buffc(21*80+1:22*80),'(''PA_V3   = '',a20)') HDR(13)
      write(buffc(22*80+1:23*80),'(''DATE-OBS= '',a20)') HDR(14)
      write(buffc(23*80+1:24*80),'(''TIME-OBS= '',a20)') HDR(15)
      write(buffc(24*80+1:25*80),'(''EXPTIME = '',a20)') HDR(16)
      write(buffc(25*80+1:26*80),'(''ROOTNAME= '',a20)') HDR(17)
      write(buffc(26*80+1:27*80),'(''TARGNAME= '',a20)') HDR(18)
      write(buffc(27*80+1:28*80),'(''RA_TARG = '',a20)') HDR(19)
      write(buffc(28*80+1:29*80),'(''DEC_TARG= '',a20)') HDR(20)
      write(buffc(29*80+1:30*80),'(''PROPOSID= '',a20)') HDR(21)
      write(buffc(30*80+1:31*80),'(''FILTER1 = '',a20)') HDR(22)
      write(buffc(31*80+1:32*80),'(''FILTER2 = '',a20)') HDR(23)
      write(buffc(32*80+1:33*80),'(''VAFACTOR= '',a20)') HDR(24)
      write(buffc(33*80+1:34*80),'(''COMMENT   '')')
      write(buffc(34*80+1:35*80),'(''COMMENT   '')')
      write(buffc(35*80+1:36*80),'(''END       '')')
 
      write(10,rec=i,iostat=ios) buffc
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      nbper  = 2*PXDIMX*PXDIMY
      npt    =   PXDIMX*PXDIMY
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
 
      do i = i1, i2, 1
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/4 + 1
         np2 = (nbyteE-nbyte1)/4 + 1
         np1 = (nbyte0-nbyte1)/2 + 1
         np2 = (nbyteE-nbyte1)/2 + 1
        !if (i.lt.0010) print*,'i: ',i,i1,i2,np1,np2,nbyte0
        !if (i.gt.2900) print*,'i: ',i,i1,i2,np1,np2,nbyte0
         call pix2buff_i2(buffb,pix,np1,npt)
         write(10,rec=i,iostat=ios) buffc
         enddo
 
      close(10)
 
      return
 
 900  continue
      print*,'WRITFITS.f ERROR'
      stop
 
      end
 
 
      subroutine pix2buff_i2(buff,pix,n1,nt)
      implicit none
      byte buff(2880)
      integer*2 pix(*)
      integer n1,nt
 
      byte b(2)
      integer*2 ii
      equivalence(ii,b)
 
      integer i, npu, nbu
      logical islinux
 
      do i = 1, 1440
         npu = n1+i-1
         nbu = (i-1)*2
         if (npu.ge.1.and.npu.le.nt) ii = pix(npu)
         if (.not.(islinux())) then
            buff(nbu+1) = b(1)
            buff(nbu+2) = b(2)
           else
            buff(nbu+1) = b(2)
            buff(nbu+2) = b(1)
            endif
         enddo
 
      return
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/WCS/inside_poly.f"
c****                                         
c*********************************************
                                              
 
 
c-----------------------------------------------
c
c reports whether a point is inside a convex
c polygon
c
      logical function inside_poly(x,y,xl,yl,Ns)
      implicit none
 
      real*8 x,y
      integer Ns
      real*8 xl(Ns), yl(Ns)
 
      integer n
      real*8  cross1
      real*8  crossN
 
      inside_poly = .true.
      cross1 = (x-xl(Ns))*(yl(1)-yl(Ns))
     .       - (y-yl(Ns))*(xl(1)-xl(Ns))
      do n = 1, Ns-1
           crossN = (x-xl(N))*(yl(N+1)-yl(N))
     .            - (y-yl(N))*(xl(N+1)-xl(N))
           if (cross1*crossN.le.0) inside_poly = .false.
           enddo
      return
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/ROUTINES/TRANS/glob_fit6nrDPU_REJ.f"
c****                                         
c*********************************************
                                              
 
 
c---------------------------------------------------------------------
c
c    RMS2D:  report the fitting residual
c SIG_CLIP:  the sigma clipping threshold
c  RES_MIN:  never reject anything within this residual (where to stop)
c  RES_MAX:  never accept anything within this residual
c
c---------------------------------------------------------------------
 
      subroutine glob_fit6nrDPU_REJ(x1,y1,x2,y2,uu,Ns,Nu,           ! added Nu
     .                              A,B,C,D,x1o,y1o,x2o,y2o,
     .                              RMS2D,SIG_CLIP,RES_MIN,RES_MAX) ! added sev
      implicit none
 
      real*8 x1(*), y1(*)
      real*8 x2(*), y2(*)
      integer uu(*)
      integer Ns
      integer Nu
      real*8 A, B, C, D
      real*8 x1o, y1o, x2o, y2o
      real   RMS2D
      real   SIG_CLIP
      real   RES_MIN
      real   RES_MAX
 
      integer n, Us
      real dx, dy, dd
 
      real    dmax
      integer nmax
      real    drms, drmsX
 
      if (Ns.gt.99999) stop 'Ns.gt.99999'
 
      drms = 0.
 
    1 continue
      Us = 0
      do N = 1, Ns
         Us = Us + uu(N)
         enddo
      if (Us.le.3) goto 3
 
      call glob_fit6nrDPU(x1,y1,x2,y2,uu,Ns,
     .                    A,B,C,D,x1o,y1o,x2o,y2o)
 
      drms = 0.
      dmax = 0.
      nmax = 0
      do N = 1, Ns
         if (uu(N).eq.1) then
             dx = SNGL(x2(N)-x2o-A*(x1(N)-x1o)-B*(y1(N)-y1o))
             dy = SNGL(y2(N)-y2o-C*(x1(N)-x1o)-D*(y1(N)-y1o))
             dd = sqrt(dx**2+dy**2)
             if (dd.gt.dmax) then
                 dmax = dd
                 nmax = n
                 endif
             drms = drms + dd**2
             endif
         enddo
 
      drmsX  = drms - dmax**2
      drmsX  = sqrt(drmsX/(Us-2))
      drms   = sqrt(drms /(Us-1))
 
      if (dmax.le.RES_MIN) goto 3 ! everything is all within specs
 
      if (dmax.gt.RES_MAX) then ! one point is beyond specs
          uu(nmax) = 0
          goto 1
          endif
 
      if (dmax.gt.drmsX*SIG_CLIP) then ! sigma-clip this point
          uu(nmax) = 0
          goto 1
          endif
 
    3 continue
      call glob_fit6nrDPU(x1,y1,x2,y2,uu,Ns,
     .                    A,B,C,D,x1o,y1o,x2o,y2o)
 
      Nu = Us
      RMS2D = drms
 
 
      return
      end
 
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/ROUTINES/TRANS/glob_fit6nrDPU.f"
c****                                         
c*********************************************
                                              
 
 
      subroutine glob_fit6nrDPU(x1,y1,x2,y2,uu,Ns,
     .                          A,B,C,D,x1o,y1o,x2o,y2o)
      implicit none
 
      real*8 x1(*), y1(*)
      real*8 x2(*), y2(*)
      integer uu(*)
      integer Ns
      real*8 A, B, C, D
      real*8 x1o, y1o, x2o, y2o
 
      real*8 sxx, sx, swx, szx
      real*8 syy, sy, swy, szy
      real*8 sw,  sz, sxy
      real*8 dlta
      real*8 dsxx, dsyy, dsxy
      real*8 dswx, dswy, dszx, dszy
 
      integer n, nnn, NNu
 
      real*8 x1a, y1a, x1b, y1b
      real*8 x2a, y2a, x2b, y2b
      real*8 x2ae,y2ae
      real*8 dx1, dy1
      real*8 dx2, dy2
 
      if (Ns.lt.3) goto 999
 
      x1a = 0.
      x1b = 0.
      y1a = 0.
      y1b = 0.
 
      x2a = 0.
      x2b = 0.
      y2a = 0.
      y2b = 0.
 
      x1o = 0
      y1o = 0
      x2o = 0
      y2o = 0
      nnn = 0
      do n = 1, Ns
           x1o = x1o + uu(n)*x1(n)
           y1o = y1o + uu(n)*y1(n)
           x2o = x2o + uu(n)*x2(n)
           y2o = y2o + uu(n)*y2(n)
           nnn = nnn + uu(n)
           enddo
 
      x1o = x1o/NNN
      y1o = y1o/NNN
      x2o = x2o/NNN
      y2o = y2o/NNN
 
      if (NNN.eq.1) then
         A = 1.000
         B = 0.000
         C = 0.000
         D = 1.000
         return
         endif
 
      if (NNN.eq.2) then
          NNu = 0
          do N = 1, Ns
             if (uu(N).eq.1.and.NNu.eq.2) stop 'prob:  NNu.eq.2'
             if (uu(N).eq.1.and.NNu.eq.1) then
                 x1b = x1(n)
                 y1b = y1(n)
                 x2b = x2(n)
                 y2b = y2(n)
                 NNu = 2
                 endif
             if (uu(N).eq.1.and.NNu.eq.0) then
                 x1a = x1(n)
                 y1a = y1(n)
                 x2a = x2(n)
                 y2a = y2(n)
                 NNu = 1
                 endif
             enddo
          x1o = (x1a+x1b)/2.0
          y1o = (y1a+y1b)/2.0
          x2o = (x2a+x2b)/2.0
          y2o = (y2a+y2b)/2.0
          dx1 = (x1b-x1a)
          dy1 = (y1b-y1a)
          dx2 = (x2b-x2a)
          dy2 = (y2b-y2a)
          A =  (dx1*dx2+dy1*dy2)/(dx1**2+dy1**2)
          B =  (dy1*dx2-dx1*dy2)/(dx1**2+dy1**2)
          C = -B
          D =  A
          x2ae = x2o + A*(x1a-x1o) + B*(y1a-y1o)
          y2ae = y2o + C*(x1a-x1o) + D*(y1a-y1o)
          return
          endif
 
 
      sxx = 0.0
      sx  = 0.0
      syy = 0.0
      sy  = 0.0
      swx = 0.0
      swy = 0.0
      szx = 0.0
      szy = 0.0
      sw  = 0.0
      sz  = 0.0
      sxy = 0.0
      do n = 1, Ns
           sxy = sxy + uu(n)*(x1(n)-x1o)*(y1(n)-y1o)
           sxx = sxx + uu(n)*(x1(n)-x1o)*(x1(n)-x1o)
           sx  = sx  + uu(n)*(x1(n)-x1o)
           syy = syy + uu(n)*(y1(n)-y1o)*(y1(n)-y1o)
           sy  = sy  + uu(n)*(y1(n)-y1o)
           swx = swx + uu(n)*(x2(n)-x2o)*(x1(n)-x1o)
           swy = swy + uu(n)*(x2(n)-x2o)*(y1(n)-y1o)
           sw  = sw  + uu(n)*(x2(n)-x2o)
           szx = szx + uu(n)*(y2(n)-y2o)*(x1(n)-x1o)
           szy = szy + uu(n)*(y2(n)-y2o)*(y1(n)-y1o)
           sz  = sz  + uu(n)*(y2(n)-y2o)
           enddo
 
      dsxx = sx*sx - NNN*sxx
      dsyy = sy*sy - NNN*syy
      dsxy = sx*sy - NNN*sxy
 
      dlta = dsxx*dsyy - dsxy*dsxy
 
      if (dlta.eq.0) goto 999
 
      dswx = sw*sx - NNN*swx
      dswy = sw*sy - NNN*swy
      dszx = sz*sx - NNN*szx
      dszy = sz*sy - NNN*szy
 
      A =   (dswx*dsyy-dswy*dsxy)/dlta
      B =   (dswy*dsxx-dswx*dsxy)/dlta
      C =   (dszx*dsyy-dszy*dsxy)/dlta
      D =   (dszy*dsxx-dszx*dsxy)/dlta
 
 
      return
 
 999  continue
      A  = 1
      B  = 0
      C  = 0
      D  = 1
      x1o = 0
      y1o = 0
      x2o = 0
      y2o = 0
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/ROUTINES/BARSIG/rbarsigs.f"
c****                                         
c*********************************************
                                              
 
 
      subroutine rbarsigs(xlist,NTOT,bar,sig,NUSE,SIGCLIP)
      implicit none
 
      integer NTOT
      real xlist(NTOT)
      real bar
      real sig
      integer NUSE
      real SIGCLIP
 
      integer n
      real*8    bsum, ssum
      integer nsum, nsumo
      integer NIT
 
 
      nsum = 0
      bar = 0.e0
      sig = 9e9
      do NIT = 1, 20
         nsumo = nsum
         bsum = 0.
         ssum = 0.
         nsum = 0.
         do n = 1, NTOT
            if (abs(xlist(n)-bar).le.SIGCLIP*sig) then
               bsum = bsum + xlist(n)
               ssum = ssum + abs(xlist(n)-bar)
               nsum = nsum + 1
               endif
            enddo
         if (nsum.gt.0) bar = SNGL(bsum / nsum)
         if (nsum.gt.1) sig = SNGL(ssum/(nsum-1))
         if (nsum.lt.0.35*NTOT.and.NIT.ge.4) goto 1
         if (nsum.eq.nsumo.and.NIT.ge.4) goto 1
         enddo
  1   continue
      NUSE = nsum
      if (nsum.le.1) sig = 0.999
 
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/GDC/STDGDCs/extract_stdgc.f"
c****                                         
c*********************************************
                                              
 
c------------------------------------------------------
c
c
c
      subroutine extract_stdgc(xr,yr,xc,yc,mc,FILE_STDGC,SENSE)
      implicit none
 
      real*8 xr, yr
      real*8 xc, yc, mc
      character*200 FILE_STDGC
      integer SENSE
 
      save
 
      logical FIRST
      data    FIRST/.true./
 
      character*200 FILE_STDGCu
      character*200 FILE_STDGCc
 
      integer   NDIM_XGC, NDIM_YGC
      integer*4 XGC_0, NDIM_XCG
      integer*4 YGC_0, NDIM_YCG
 
      integer i, ios
      integer k
      character*08   FIELD
      character*20   STREAM
      character*2880 buffc
 
      real*8  fx,   fy
      real*8  xcu,  ycu
      integer ipix, jpix
 
      integer*4, dimension(:,:), allocatable :: xgc
      integer*4, dimension(:,:), allocatable :: ygc
      integer*2, dimension(:,:), allocatable :: mgc
      integer*4, dimension(:,:), allocatable :: xcg
      integer*4, dimension(:,:), allocatable :: ycg
 
      if (FILE_STDGC(1:4).eq.'NONE'.or.
     .    FILE_STDGC(1:4).eq.'NULL') then
         if (SENSE.eq.0) then
            mc = 0.
            return
            endif
         if (SENSE.eq.1) then ! forward...
            xc = xr
            yc = yr
            return
            endif
         if (SENSE.eq.-1) then ! forward...
            xr = xc
            yr = yc
            return
            endif
         endif
 
      if (FIRST) FILE_STDGCc = ' '
 
      FILE_STDGCu = FILE_STDGC
      do k = 196,2,-1
         if (FILE_STDGCu(k:k+4).eq.'.fits')
     .       FILE_STDGCu = FILE_STDGCu(1:k+4)
         enddo
 
      if (FILE_STDGCc.ne.FILE_STDGCu) then
          FILE_STDGCc = FILE_STDGCu
c         print*,'EXTRACT_STDGC --- OPEN(10): ',FILE_STDGCc
          open(10,file=FILE_STDGCc,
     .         status='old',
     .         err   =900,
     .         recl  =2880,
     .         form  ='UNFORMATTED',
     .         access='DIRECT')
          i = 1
          read(10,rec=i,iostat=ios) buffc
          do k = 0, 35
             FIELD  = buffc(k*80+01:k*80+08)
             STREAM = buffc(k*80+11:k*80+30)
             if (FIELD.eq.'NDIM_XGC') read(STREAM,*) NDIM_XGC
             if (FIELD.eq.'NDIM_YGC') read(STREAM,*) NDIM_YGC
             if (FIELD.eq.'NDIM_XCG') read(STREAM,*) NDIM_XCG
             if (FIELD.eq.'NDIM_YCG') read(STREAM,*) NDIM_YCG
             if (FIELD.eq.'XGC_0   ') read(STREAM,*) XGC_0
             if (FIELD.eq.'YGC_0   ') read(STREAM,*) YGC_0
             enddo
          if (.not.FIRST) deallocate(xgc,ygc,mgc,xcg,ycg)
          allocate(xgc(NDIM_XGC,NDIM_YGC))
          allocate(ygc(NDIM_XGC,NDIM_YGC))
          allocate(mgc(NDIM_XGC,NDIM_YGC))
          allocate(xcg(NDIM_XCG,NDIM_YCG))
          allocate(ycg(NDIM_XCG,NDIM_YCG))
          i = 2
          call getext_i4(xgc,NDIM_XGC,NDIM_YGC,i)
          call getext_i4(ygc,NDIM_XGC,NDIM_YGC,i)
          call getext_i2(mgc,NDIM_XGC,NDIM_YGC,i)
          call getext_i4(xcg,NDIM_XCG,NDIM_YCG,i)
          call getext_i4(ycg,NDIM_XCG,NDIM_YCG,i)
          close(10)
          FIRST = .false.
          endif
 
      if (SENSE.eq.1) then ! forward...
          ipix = int(xr)
          jpix = int(yr)
          xc = -9999
          yc = -9999
 
          if (ipix.lt.       1-1) return
          if (ipix.gt.NDIM_XGC+1) return
          if (jpix.lt.       1-1) return
          if (jpix.gt.NDIM_YGC+1) return
 
          if (ipix.lt.1) ipix = 1
          if (jpix.lt.1) jpix = 1
          if (ipix.gt.NDIM_XGC-1) ipix = NDIM_XGC-1
          if (jpix.gt.NDIM_YGC-1) jpix = NDIM_YGC-1
 
          fx = xr-ipix
          fy = yr-jpix
          xc = (1-fx)*(1-fy)*xgc(ipix  ,jpix  )/1d5
     .       + (1-fx)*( fy )*xgc(ipix  ,jpix+1)/1d5
     .       + ( fx )*(1-fy)*xgc(ipix+1,jpix  )/1d5
     .       + ( fx )*( fy )*xgc(ipix+1,jpix+1)/1d5
          yc = (1-fx)*(1-fy)*ygc(ipix  ,jpix  )/1d5
     .       + (1-fx)*( fy )*ygc(ipix  ,jpix+1)/1d5
     .       + ( fx )*(1-fy)*ygc(ipix+1,jpix  )/1d5
     .       + ( fx )*( fy )*ygc(ipix+1,jpix+1)/1d5
          return
          endif
 
      if (SENSE.eq.0) then ! find the zeropoint...
          ipix = int(xr)
          jpix = int(yr)
          mc = 0.000
          if (ipix.lt.         1) return
          if (ipix.gt.NDIM_XGC-1) return
          if (jpix.lt.         1) return
          if (jpix.gt.NDIM_YGC-1) return
          fx = xr-ipix
          fy = yr-jpix
          mc = (1-fx)*(1-fy)*mgc(ipix  ,jpix  )/1d4
     .       + (1-fx)*( fy )*mgc(ipix  ,jpix+1)/1d4
     .       + ( fx )*(1-fy)*mgc(ipix+1,jpix  )/1d4
     .       + ( fx )*( fy )*mgc(ipix+1,jpix+1)/1d4
          return
          endif
 
      if (SENSE.eq.-1) then ! reverse
          xcu =  1 + (xc-XGC_0)
          ycu =  1 + (yc-YGC_0)
          ipix = int(xcu)
          jpix = int(ycu)
          xr = -1
          yr = -1
          if (ipix.lt.         1) return
          if (ipix.gt.NDIM_XCG-1) return
          if (jpix.lt.         1) return
          if (jpix.gt.NDIM_YCG-1) return
          if (xcg(ipix  ,jpix  ).le.-2000000000) return
          if (xcg(ipix+1,jpix  ).le.-2000000000) return
          if (xcg(ipix  ,jpix+1).le.-2000000000) return
          if (xcg(ipix+1,jpix+1).le.-2000000000) return
           fx = xcu-ipix
           fy = ycu-jpix
           xr = (1-fx)*(1-fy)*xcg(ipix  ,jpix  )/1d5
     .        + (1-fx)*( fy )*xcg(ipix  ,jpix+1)/1d5
     .        + ( fx )*(1-fy)*xcg(ipix+1,jpix  )/1d5
     .        + ( fx )*( fy )*xcg(ipix+1,jpix+1)/1d5
           yr = (1-fx)*(1-fy)*ycg(ipix  ,jpix  )/1d5
     .        + (1-fx)*( fy )*ycg(ipix  ,jpix+1)/1d5
     .        + ( fx )*(1-fy)*ycg(ipix+1,jpix  )/1d5
     .        + ( fx )*( fy )*ycg(ipix+1,jpix+1)/1d5
           return
           endif
 
      print*,'ILLEGAL SENSE',SENSE
      stop
 
  900 continue
      print*,'file open error'
      print*,'  FILE_STDGCc: ',FILE_STDGCc
      print*,'  FILE_SGDGC : ',FILE_STDGC
      stop
      end
 
 
c--------------------------------------------------------------------
c
c
      subroutine getext_i4(pix,NDIMX,NDIMY,i)
      implicit none
 
      integer NDIMX, NDIMY
      integer*4 PIX(NDIMX,NDIMY)
      integer i
 
      integer o
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer ifirst, i1, i2, ios
      integer np1, np2, npt
      integer nbu, npu
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffc,buffb)
 
      byte      b(4)
      integer*4 ii
      equivalence(ii,b)
 
      integer ipix, jpix
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      nbper  = 4*NDIMX*NDIMY
      npt    =   NDIMX*NDIMY
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
      do i = i1, i2, 1
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/4 + 1
         np2 = (nbyteE-nbyte1)/4 + 1
         read(10,rec=i,iostat=ios) buffc
         do o = 1, 720
            npu = np1+o-1
            nbu = (o-1)*4
            ii  = 0
            ipix = npu - (npu-1)/NDIMX*NDIMX
            jpix =   1 + (npu-1)/NDIMX
            if (ipix.ge.1.and.ipix.le.NDIMX.and.
     .          jpix.ge.1.and.jpix.le.NDIMY) then
                if (.not.(_LINUX_)) then
                   b(1) = buffb(nbu+1)
                   b(2) = buffb(nbu+2)
                   b(3) = buffb(nbu+3)
                   b(4) = buffb(nbu+4)
                   endif
                if ((_LINUX_)) then
                   b(4) = buffb(nbu+1)
                   b(3) = buffb(nbu+2)
                   b(2) = buffb(nbu+3)
                   b(1) = buffb(nbu+4)
                   endif
                pix(ipix,jpix) = ii
                endif
            enddo!o
         enddo!i
 
      i = i2 + 1
 
      return
      end
 
 
 
c--------------------------------------------------------------------
c
c
      subroutine getext_i2(pix,NDIMX,NDIMY,i)
      implicit none
 
      integer NDIMX, NDIMY
      integer*2 PIX(NDIMX,NDIMY)
      integer i
 
      integer o
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer ifirst, i1, i2, ios
      integer np1, np2, npt
      integer nbu, npu
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffc,buffb)
 
      byte      b(2)
      integer*2 ii
      equivalence(ii,b)
 
      integer ipix, jpix
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      nbper  = 2*NDIMX*NDIMY
      npt    =   NDIMX*NDIMY
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
      do i = i1, i2, 1
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/2 + 1
         np2 = (nbyteE-nbyte1)/2 + 1
         read(10,rec=i,iostat=ios) buffc
         do o = 1, 1440
            npu = np1+o-1
            nbu = (o-1)*2
            ii  = 0
            ipix = npu - (npu-1)/NDIMX*NDIMX
            jpix =   1 + (npu-1)/NDIMX
            if (ipix.ge.1.and.ipix.le.NDIMX.and.
     .          jpix.ge.1.and.jpix.le.NDIMY) then
                if (.not.(_LINUX_)) then
                   b(1) = buffb(nbu+1)
                   b(2) = buffb(nbu+2)
                   endif
                if ((_LINUX_)) then
                   b(2) = buffb(nbu+1)
                   b(1) = buffb(nbu+2)
                   endif
                pix(ipix,jpix) = ii
                endif
            enddo!o
         enddo!i
 
      i = i2 + 1
 
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/GDC/STDGDCs/xryr2mc_stdgc.f"
c****                                         
c*********************************************
                                              
 
c------------------------------------------------------
c
c mc should be the pixel area correction
c
c
      subroutine xryr2mc_stdgc(xr,yr,mc,FILE_STDGC)
      implicit none
 
      real*8 xr, yr
      real*8 mc
      character*(*) FILE_STDGC
      real*8 xc, yc
 
      call extract_stdgc(xr,yr,xc,yc,mc,FILE_STDGC,0)  ! 0 = pix area correction
 
      return
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/GDC/STDGDCs/xryr2xcyc_stdgc.f"
c****                                         
c*********************************************
                                              
 
c------------------------------------------------------
c
      subroutine xryr2xcyc_stdgc(xr,yr,xc,yc,FILE_STDGC)
      implicit none
 
      real*8 xr, yr
      real*8 xc, yc
      character*(*) FILE_STDGC
      real*8 mc
 
      call extract_stdgc(xr,yr,xc,yc,mc,FILE_STDGC,1)  ! 1 = forward
 
      return
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/GDC/STDGDCs/xcyc2xryr_stdgc.f"
c****                                         
c*********************************************
                                              
 
c------------------------------------------------------
c
      subroutine xcyc2xryr_stdgc(xc,yc,xr,yr,FILE_STDGC)
      implicit none
 
      real*8 xc, yc
      real*8 xr, yr
      character*(*) FILE_STDGC
      real*8 mc
 
      call extract_stdgc(xr,yr,xc,yc,mc,FILE_STDGC,-1)  ! -1 = reverse
 
      return
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/STRING/strlen.f"
c****                                         
c*********************************************
                                              
 
      integer function strlen(STRING)
      implicit none
 
      character*(*) STRING
 
      integer i
 
      i = 0
    1 continue
      i = i + 1
      if (STRING(i:i).eq.' ') then
          strlen = i-1
          return
          endif
      goto 1
 
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/WCS/xy2r.f"
c****                                         
c*********************************************
                                              
 
c-----------------------------------------------
c
c This routine will take a point (x,y) in the
c tangent plane, and the coordinate of the
c tangent-point (ra0,dec0), and will compute
c the ra coordinate that corresponds to (x,y).
c
      real*8 function xy2r(x,y,r0,d0)
      implicit none
 
      real*8 x, y
      real*8 r0,d0
 
c     real*8 cosde, sinde
      real*8 cosd0, sind0
      real*8 tandr, dr
      real*8 xrad, yrad
 
      xrad = x/180.0d0*(_PI_)
      yrad = y/180.0d0*(_PI_)
 
      cosd0 = cos(d0*(_PI_)/180.0d0)
      sind0 = sin(d0*(_PI_)/180.0d0)
 
      tandr = xrad/(cosd0-yrad*sind0)
      dr    = atan(tandr)
      xy2r  = r0 + dr*180.0d0/(_PI_)
 
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/WCS/xy2d.f"
c****                                         
c*********************************************
                                              
 
c-----------------------------------------------
c
c This routine will take a point (x,y) in the
c tangent plane, and the coordinate of the
c tangent-point (ra0,dec0), and will compute
c the ra coordinate that corresponds to (x,y).
c
      real*8 function xy2d(x,y,r0,d0)
      implicit none
 
      real*8 x, y
      real*8 r0,d0
 
      real*8 cosd0, sind0
      real*8 tandr, dr
      real*8 cosdr
      real*8 tande
      real*8 xrad, yrad
 
      if (r0.gt.0) continue
 
      xrad = x/180.0d0*(_PI_)
      yrad = y/180.0d0*(_PI_)
 
      cosd0 = cos(d0*(_PI_)/180.0d0)
      sind0 = sin(d0*(_PI_)/180.0d0)
 
      tandr = xrad/(cosd0-yrad*sind0)
      dr    = atan(tandr)
 
      cosdr = cos(dr)
 
      tande = cosdr*(sind0+yrad*cosd0)/
     .              (cosd0-yrad*sind0)
 
      !print*,'---> tande: ',tande
      xy2d = atan(tande)*180.0d0/(_PI_)
 
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/WCS/rd2x.f"
c****                                         
c*********************************************
                                              
c-----------------------------------------------
c
c This routine will take an (ra,dec) and a
c tangent-point (ra0,dec0) and will compute
c the x coordinate in the tangent plane.  The
c tangent-plane x-axis is directed along -RA.
c The units for x will be degrees, just lilke
c those for RA and DEC.
c
      real*8 function rd2x(r,d,r0,d0)
      implicit none
 
      real*8 r, d
      real*8 r0,d0
 
      real*8 cosra, sinra
      real*8 cosde, sinde
      real*8 cosd0, sind0
      real*8 rrrr
      real*8 xrad
      real*8 x, y, z
      real*8 xx,yy,zz
 
      cosra = cos((r-r0)*(_PI_)/180.0d0)
      sinra = sin((r-r0)*(_PI_)/180.0d0)
 
      cosde = cos(d *(_PI_)/180.0d0)
      sinde = sin(d *(_PI_)/180.0d0)
      cosd0 = cos(d0*(_PI_)/180.0d0)
      sind0 = sin(d0*(_PI_)/180.0d0)
 
      rrrr = sind0*sinde + cosd0*cosde*cosra
 
 
      xrad = cosde*sinra/rrrr
      rd2x = xrad*180.0d0/_PI_
 
      x  = cosde*cos(r *(_PI_)/180.0d0)
      y  = cosde*sin(r *(_PI_)/180.0d0)
      z  = sinde
      xx = cosd0*cos(r0*(_PI_)/180.0d0)
      yy = cosd0*sin(r0*(_PI_)/180.0d0)
      zz = sind0
 
      if (x*xx + y*yy + z*zz.lt.0) rd2x = 90.0
 
 
      return
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/WCS/rd2y.f"
c****                                         
c*********************************************
                                              
 
c-----------------------------------------------
c
c This routine will take an (ra,dec) and a
c tangent-point (ra0,dec0) and will compute
c the y coordinate in the tangent plane.  The
c tangent-plane y-axis is directed along +DEC.
c The units for y will be degrees, just lilke
c those for RA and DEC.
c
      real*8 function rd2y(r,d,r0,d0)
      implicit none
 
      real*8 r, d
      real*8 r0,d0
 
      real*8 cosra, sinra
      real*8 cosde, sinde
      real*8 cosd0, sind0
      real*8 rrrr
      real*8 yrad
      real*8 x, y, z
      real*8 xx,yy,zz
 
      cosra = cos((r-r0)*(_PI_)/180.0d0)
      sinra = sin((r-r0)*(_PI_)/180.0d0)
      cosde = cos(d *(_PI_)/180.0d0)
      sinde = sin(d *(_PI_)/180.0d0)
      cosd0 = cos(d0*(_PI_)/180.0d0)
      sind0 = sin(d0*(_PI_)/180.0d0)
 
      rrrr = sind0*sinde + cosd0*cosde*cosra
 
      yrad = (cosd0*sinde-sind0*cosde*cosra)/rrrr
      rd2y = yrad*180.0d0/_PI_
 
      x  = cosde*cos(r *(_PI_)/180.0d0)
      y  = cosde*sin(r *(_PI_)/180.0d0)
      z  = sinde
      xx = cosd0*cos(r0*(_PI_)/180.0d0)
      yy = cosd0*sin(r0*(_PI_)/180.0d0)
      zz = sind0
 
      if (x*xx + y*yy + z*zz.lt.0) rd2y = 90
 
      return
      end
 
 
 
