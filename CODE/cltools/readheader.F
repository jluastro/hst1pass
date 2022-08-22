 
 
 
      program imginfo
      implicit none
 
      integer i
      character filename*200
      integer NCHAR
 
      integer iargc, NARG, NARGs
      integer NIM
 
      integer iEXPT
      real    rEXPT
      integer iPA
      real    rPA
      integer iPID
      character*6 FILTSTR
      character*12 DATESTRIN
      character*9 ROOTSTR
      real*8 rRA, hRA
      real*8 rDEC, aDEC
      real*8 rRA_TARG
      real*8 rDEC_TARG
      character*16 APERNAME
 
      character*8 DATESTR
      character*8 TIMESTR
 
      integer iRAH, iRAM, iRAS, iFRA
      integer iDEG, iMIN, iSEC, iFDC
      character decsgn
      real*8  rRA0, rDC0
      real    dRA , dDC
      real    rSUNA
      real    rTEMP
      real    vaFAC
 
      integer iDAY
      integer iMON
      integer iYIR
      integer iHRS
 
      real*8  rDAT
 
      character*8  field
      character*20 stream
 
      integer ios, k
 
      character*2880 buff
 
      integer nread
 
 
      character*200 FILELIST
      character*200 FILENAME_N(999999)
 
      character*11 OBJNAME
 
      integer DOY(12)
      data DOY / 000, 031, 059, 090, 120, 151,
     .           181, 212, 243, 273, 304, 334/
 
      logical DEBUG
      real rPT1, rPT2
      integer NEND
      integer NENDs
      integer NAXIS
 
      character*08 KEYWORD
      character*20 KEYSTREAM
 
      logical SHOW_REG
      logical SHOW_APER
      logical SHOW_POSTARG
      logical SHOW_KEYWORD
      logical DO_LIST
 
      DEBUG = .false.
      SHOW_REG     = .true.
      SHOW_APER    = .false.
      SHOW_POSTARG = .false.
      SHOW_KEYWORD = .false.
      DO_LIST      = .false.
      KEYWORD      = 'NONE'
 
      if (iargc().eq.0) then
         print*,' '
         print*,'TO USE THIS ROUTINE, YOU NEED TO'
         print*,'GIVE IT SOME FITS IMAGES ON THE'
         print*,'COMMAND LINE, USING ASTERISKS IF'
         print*,'DESIRED.'
         print*,' '
         print*,'FOR EXAMPLE, '
         print*,' '
         print*,'% readheader.e [POSTARG/APER/KEY=KEYWORD] *.fits     '
         print*,'% readheader.e [POSTARG/APER/KEY=KEYWORD] LIST=a.txt '
         print*,' '
         print*,'WILL TABULATE THE INFO FOR ALL THE '
         print*,'FITS FILES IN THE CURRENT DIRECTORY. '
         print*,' '
         stop
         endif
 
      NARGs = 0
      do NARG = 1, iargc()
         call getarg(NARG,FILENAME)
 
         if (FILENAME(1:5).eq.'LIST=') then
             DO_LIST = .true.
             FILELIST = FILENAME(6:200)
             print*,'---> FILELIST: ',FILELIST(1:40)
             open(19,file=FILENAME(6:200),status='old')
    1        read(19,'(a200)',end=77) FILENAME_N(NARGs+1)
             NARGs = NARGs + 1
             if (NARGs.ge.999999) stop 'TOO MANY FILES IN LIST >999999'
             goto 1
             endif
         if (FILENAME(1:5).eq.'DEBUG') then
            DEBUG = .true.
            goto 77
            endif
         if (FILENAME(1:4).eq.'APER') then
            SHOW_APER = .true.
            SHOW_REG  = .false.
            goto 77
            endif
         if (FILENAME(1:4).eq.'KEY=') then
            SHOW_KEYWORD = .true.
            SHOW_REG     = .false.
            KEYWORD      = FILENAME(5:12)
            goto 77
            endif
         if (FILENAME(1:7).eq.'POSTARG') then
            SHOW_POSTARG = .true.
            SHOW_REG  = .false.
            goto 77
            endif
         if (DO_LIST) stop 'cannot have list and *fits'
         NARGs = NARGs + 1
         FILENAME_N(NARGs) = FILENAME
   77    continue
         enddo
 
 
      write(*,'(''# '')')
      write(*,'(''# INPUT-LIST NARGs: '',i4,1x,a80)')
     .              NARGs,FILENAME(1:80)
      write(*,'(''# '')')
    3 continue
 
      if (DEBUG) then
          print*,'           '
          print*,'SHOW_REG : ',SHOW_REG
          print*,'SHOW_APER: ',SHOW_APER
          print*,'SHOW_PTRG: ',SHOW_POSTARG
          print*,'  DO_LIST: ',DO_LIST
          print*,'    NARGs: ',NARGs
          print*,'SHOW_KEYW: ',SHOW_KEYWORD
          print*,'     KEYW: ',KEYWORD
          print*,'           '
          endif
 
 
      NIM = 0
      do NARG = 1, NARGs
 
         KEYSTREAM = 'NOT_FOUND'
         APERNAME  = 'NONE'
         DATESTR   = 'DATELESS'
         TIMESTR   = 'TIMELESS'
         ROOTSTR   = 'ROOTLESS'
         FILTSTR   = 'FLTLSS'
         rEXPT     = 0
         iEXPT     = 0
         rRA  = 0.
         rDEC = 0.
         iPID = 0
         iPA  = 0
         rPA  = 0
         rSUNA = 0.
         rTEMP = 0.
         vaFAC = 0.
         rPT1  = 99.99
         rPT2  = 99.99
         FILENAME = FILENAME_N(NARG)
         NIM = NIM+1
 
c-----------------------------------------------
 
c        print*,'open10: ',FILENAME
         open(10,file=FILENAME,status='old',
     .     err=900,recl=2880,form='UNFORMATTED',
     .     access='DIRECT')
 
         NENDs = 1
         NEND  = 0
         NAXIS = 0
 
         i = 0
         nread = 0
 100     continue
         i = i + 1
         read(10,rec=i,iostat=ios) buff
         if (ios.lt.0) goto 900
         do k = 0, 35, 1
            field  = buff(k*80+01:k*80+08)
            stream = buff(k*80+11:k*80+31)
            if (DEBUG)
     .      write(*,'(i2.2,1x,i2.2,'':'',a8,'' | '',a20)')
     .               i,k,field,stream
            if (stream(1:1).eq.'''') stream    = stream(2:20)
            if (field.eq.'EXTEND  '.and.
     .          stream(20:20).eq.'T') NENDs = 2
            if (field.eq.'NAXIS   ') read(stream(01:20),*) NAXIS
            if (field.eq.'TIME-OBS') TIMESTR   = stream(01:08)
            if (field.eq.'DATE-OBS') DATESTRIN = stream(01:10)
            if (field.eq.'ROOTNAME') ROOTSTR   = stream(01:10)
            if (field.eq.'FILTNAM1'.and.
     .          stream(1:1).eq.'F')  FILTSTR   = stream(01:06)
            if (field.eq.'FILTNAM1'.and.
     .          stream(1:1).eq.'F')  FILTSTR   = stream(01:06)
            if (field.eq.'FILTER1 '.and.
     .          stream(1:1).eq.'F')  FILTSTR = stream(01:06)
            if (field.eq.'FILTER2 '.and.
     .          stream(1:1).eq.'F')  FILTSTR = stream(01:06)
            if (field.eq.'FILTER  '.and.
     .          stream(1:1).eq.'F')  FILTSTR = stream(01:06)
            if (field.eq.'EXPTIME ') read(stream,*,end=109) rEXPT
            if (field.eq.'RA_TARG ') read(stream,*,end=109) rRA_TARG
            if (field.eq.'DEC_TARG') read(stream,*,end=109) rDEC_TARG
            if (field.eq.'PROPOSID') read(stream,*,end=109) iPID
            if (field.eq.'PA_V3   ') read(stream,*,end=109) rPA
            if (field.eq.'SUN_ALT ') read(stream,*,end=109) rSUNA
            if (field.eq.'UBAY3TMP') read(stream,*,end=109) rTEMP
            if (field.eq.'VAFACTOR') read(stream,*,end=109) vaFAC
            if (field.eq.'CRVAL1  '.and.rRA.eq.0)
     .                               read(stream,*,end=109) rRA
            if (field.eq.'CRVAL2  '.and.rDEC.eq.0)
     .                               read(stream,*,end=109) rDEC
            if (field.eq.'POSTARG1') read(stream,*,end=109) rPT1
            if (field.eq.'POSTARG2') read(stream,*,end=109) rPT2
            if (field.eq.'RA_APER ') read(stream,*,end=109) rRA
            if (field.eq.'DEC_APER') read(stream,*,end=109) rDEC
            if (field.eq.'CRVAL1  ') read(stream,*,end=109) rRA
            if (field.eq.'CRVAL2  ') read(stream,*,end=109) rDEC
            if (field.eq.'APERTURE') APERNAME = stream(1:16)
            if (field.eq.KEYWORD)    KEYSTREAM = buff(k*80+11:k*80+30)
 
            if (field.eq.'END     ') then
               NEND = NEND + 1
               if (DEBUG) then
                   print*,' '
                   print*,'---> END CARD...'
                   print*,'---> NEND : ',NEND
                   print*,'---> NENDs: ',NENDs
                   print*,'---> NAXIS: ',NAXIS
                   print*,' '
                   endif
               if (NEND.eq.NENDs.or.NAXIS.ne.0) goto 101
               endif
 109        continue
            enddo
         goto 100
 101     continue
         close(10)
 
         if (abs(rRA).lt.0.01) then
            rRA  = rRA_TARG
            rDEC = rDEC_TARG
            endif
 
         iPA = rPA
         iEXPT = int(rEXPT+0.99) ! round up
 
         if (DEBUG) print*,'DEBUG: FIND DATESTRING: ',DATESTRIN
 
         if (DATESTRIN(3:3).eq.'/') DATESTR(1:8) = DATESTRIN(1:8)
         if (DATESTRIN(5:5).eq.'-') then
            DATESTR(1:2) = DATESTRIN(9:10)
            DATESTR(3:3) = '/'
            DATESTR(4:5) = DATESTRIN(6:7)
            DATESTR(6:6) = '/'
            DATESTR(7:8) = DATESTRIN(3:4)
            endif
 
         if (DATESTR.ne.'DATELESS') then
            read(DATESTR(1:2),*) iDAY
            read(DATESTR(4:5),*) iMON
            read(DATESTR(7:8),*) iYIR
            endif
 
         if (TIMESTR.ne.'TIMELESS') then
            read(TIMESTR(1:2),*) iHRS
            read(TIMESTR(4:5),*) iMIN
            read(TIMESTR(7:8),*) iSEC
            endif
 
         if (DEBUG) then
            print*,' iYIR: ',iYIR
            print*,' iDAY: ',iDAY,DOY(iMON)
            print*,' iMON: ',iMON
            endif
 
         rDAT = iYIR + (DOY(iMON)+iDAY)/365.0
         rDAT = rDAT + (iHRS + iMIN/60. + iSEC/60./60.)/8766.
 
         if (DEBUG) print*,' rDAT: ',rDAT
 
 
 
c
c this routine will figure out which interesting object
c is closest to the observed position and will compute
c where this field is located with respect to the object.
c
         if (DEBUG) print*,'DEBUG: FIND_OBJ...'
         call find_obj(rRA,rDEC,OBJNAME,rRA0,rDC0)
 
         dRA = (rRA -rRA0)*60.0*cos(rDC0*3.14159/180)  ! in arcmin (W is +ve)
         dDC = (rDEC-rDC0)*60.0                        ! in arcmin
 
         hRA = rRA*24/360
         iRAH = int(hRA)
         iRAM = (hRA-iRAH)*60
         iRAS = (hRA-iRAH-iRAM/60.)*3600
         iFRA = (hRA-iRAH-iRAM/60.-iRAS/3600.)*36000
 
         aDEC = abs(rDEC)
         iDEG = int(aDEC)
         iMIN = (aDEC-int(aDEC))*60
         iSEC = (aDEC-int(aDEC)-iMIN/60.)*3600
         iFDC = (aDEC-int(aDEC)-iMIN/60.-iSEC/3600.)*36000
 
         decsgn = '+'
         if (rDEC.lt.0) decsgn = '-'
 
c
c find end of "fits"
c
         NCHAR = 99
         do i = 5, 99
            if (filename(i:i).eq.'s') NCHAR = i
            if (filename(i:i).eq.'S') NCHAR = i
            enddo
 
         if (DEBUG) print*,'DEBUG: OUTPUT...'
         if (NIM.eq.1.and.SHOW_REG)     write(*,111)
         if (NIM.eq.1.and.SHOW_POSTARG) write(*,121)
         if (NIM.eq.1.and.SHOW_APER)    write(*,131)
         if (NIM.eq.1.and.SHOW_KEYWORD) write(*,141) KEYWORD
 111     format('#ROOTNAME   PID  DATE_OBS TIME_OBS ',
     .          'RA_TARG    DEC_TARG',
     .          '     PAV3    FILT   EXPT ',
     .          'RDATE    VAFACTOR    dRA('')  dDEC('') OBJNAME')
 121     format('#ROOTNAME   PID  DATE_OBS TIME_OBS ',
     .          'RA_TARG    DEC_TARG',
     .          '     PAV3    FILT   EXPT ',
     .          'RDATE    VAFACTOR    POSTARG1   POSTARG2 FILENAME')
 131     format('#ROOTNAME   PID  DATE_OBS TIME_OBS ',
     .          'RA_TARG    DEC_TARG',
     .          '     PAV3    FILT   EXPT ',
     .          'RDATE    VAFACTOR  APERNAME       FILENAME')
 141     format('#ROOTNAME   PID  DATE_OBS TIME_OBS ',
     .          'RA_TARG    DEC_TARG',
     .          '     PAV3    FILT   EXPT ',
     .          'RDATE    VAFACTOR   KEYWORD=',a08)
         if (SHOW_REG)
     .       write(*,112) ROOTSTR,iPID,
     .                    DATESTR,TIMESTR,
     .                    iRAH,iRAM,iRAS,iFRA,decsgn,
     .                    iDEG,iMIN,iSEC,iFDC,
     .                    rPA, FILTSTR, iEXPT,
     .                    rDAT,vaFAC,dRA,dDC,OBJNAME,
     .                    FILENAME(1:NCHAR)
         if (SHOW_POSTARG)
     .       write(*,122) ROOTSTR,iPID,
     .                    DATESTR,TIMESTR,
     .                    iRAH,iRAM,iRAS,iFRA,decsgn,
     .                    iDEG,iMIN,iSEC,iFDC,
     .                    rPA, FILTSTR, iEXPT,
     .                    rDAT,vaFAC,rPT1,rPT2,
     .                    FILENAME(1:NCHAR)
         if (SHOW_APER)
     .       write(*,132) ROOTSTR,iPID,
     .                    DATESTR,TIMESTR,
     .                    iRAH,iRAM,iRAS,iFRA,decsgn,
     .                    iDEG,iMIN,iSEC,iFDC,
     .                    rPA, FILTSTR, iEXPT,
     .                    rDAT,vaFAC,APERNAME,
     .                    FILENAME(1:NCHAR)
         if (SHOW_KEYWORD)
     .       write(*,142) ROOTSTR,iPID,
     .                    DATESTR,TIMESTR,
     .                    iRAH,iRAM,iRAS,iFRA,decsgn,
     .                    iDEG,iMIN,iSEC,iFDC,
     .                    rPA, FILTSTR, iEXPT,
     .                    rDAT,vaFAC,KEYSTREAM,
     .                    FILENAME(1:NCHAR)
 
 112     format(1x,a9,1x,i5,1x,a8,1x,a8,
     .          1x,i2.2,'d',i2.2,'m',i2.2,'s',i1,
     .          1x,a1,i2.2,':',i2.2,':',i2.2,'.',i1,
     .          1x,f8.4,1x,a6,1x,i4,1x,
     .          f8.5,1x,f9.7,1x,f8.2,1x,f8.2,1x,a11,1x,30a)
 122     format(1x,a9,1x,i5,1x,a8,1x,a8,
     .          1x,i2.2,'d',i2.2,'m',i2.2,'s',i1,
     .          1x,a1,i2.2,':',i2.2,':',i2.2,'.',i1,
     .          1x,f8.4,1x,a6,1x,i4,1x,
     .          f8.5,1x,f9.7,1x,f10.5,1x,f10.5,1x,30a)
 132     format(1x,a9,1x,i5,1x,a8,1x,a8,
     .          1x,i2.2,'d',i2.2,'m',i2.2,'s',i1,
     .          1x,a1,i2.2,':',i2.2,':',i2.2,'.',i1,
     .          1x,f8.4,1x,a6,1x,i4,1x,
     .          f8.5,1x,f9.7,1x,a16,1x,30a)
 142     format(1x,a9,1x,i5,1x,a8,1x,a8,
     .          1x,i2.2,'d',i2.2,'m',i2.2,'s',i1,
     .          1x,a1,i2.2,':',i2.2,':',i2.2,'.',i1,
     .          1x,f8.4,1x,a6,1x,i4,1x,
     .          f8.5,1x,f9.7,1x,a20,1x,30a)
 
 777     continue
         enddo!NARGs
      stop
 
  900 continue
      print*,' '
      print*,'readheader.e ERROR EXIT. '
      print*,' '
      print*,'ONE OF THE IMAGES WAS NOT IN STANDARD'
      print*,'HST FITS FORMAT.'
      print*,' '
      write(*,'(''PROBLEM FILE: '',a80)') FILENAME
      print*,' '
      end
 
 
 
c*********************************************
c****
c**** #include "/user/jayander/FORTRAN/WCS/find_obj.f"
c****
c*********************************************
 
 
c----------------------------------------------------
c
c find the closest object of interest to this
c location; info from HARRIS' on-line catalog
c is contained within data statements.
c
      subroutine find_obj(rRA,rDC,NAME,rRA0,rDC0)
      implicit none
      real*8 rRA , rDC      ! input
      character*11 NAME     ! output (nearest obj)
      real*8 rRA0, rDC0     ! output (locn)
 
      real*8 rRAn, rDCn
      real*8 dRA,  dDC
      real*8 dDD, dmin
      integer N
      character*10 RSTRING
      character*09 DSTRING
 
      integer iRAH, iRAM
      real    rRAS
      integer iDEG, iMIN, iSEC
 
      character*35 INFO(191)
      common /INFO/INFO ! need to do this in LINUX to preserve
      data INFO /
     .'NGC0104 47T 00 24 05.2 -72 04 51.00',
     .'NGC0288     00 52 47.5 -26 35 24.00',
     .'NGC0362     01 03 14.3 -70 50 54.00',
     .'NGC1261     03 12 15.3 -55 13 01.00',
     .'Pal1        03 33 23.0 +79 34 50.00',
     .'AM1     E01 03 55 02.7 -49 36 52.00',
     .'Eridanus    04 24 44.5 -21 11 13.00',
     .'Pal2        04 46 05.9 +31 22 51.00',
     .'NGC1851     05 14 06.3 -40 02 50.00',
     .'NGC1904 M79 05 24 10.6 -24 31 27.00',
     .'NGC2298     06 48 59.2 -36 00 19.00',
     .'NGC2419     07 38 08.5 +38 52 55.00',
     .'Pyxis       09 07 57.8 -37 13 17.00',
     .'NGC2808     09 12 02.6 -64 51 47.00',
     .'E3          09 20 59.3 -77 16 57.00',
     .'Pal3        10 05 31.4 +00 04 17.00',
     .'NGC3201     10 17 36.8 -46 24 40.00',
     .'Pal4        11 29 16.8 +28 58 25.00',
     .'NGC4147     12 10 06.2 +18 32 31.00',
     .'NGC4372     12 25 45.4 -72 39 33.00',
     .'Rup106      12 38 40.2 -51 09 01.00',
     .'NGC4590 M68 12 39 28.0 -26 44 34.00',
     .'NGC4833     12 59 35.0 -70 52 29.00',
     .'NGC5024 M53 13 12 55.3 +18 10 09.00',
     .'NGC5053     13 16 27.0 +17 41 53.00',
     .'NGC5139 OMC 13 26 45.9 -47 28 37.00',
     .'NGC5272 M03 13 42 11.2 +28 22 32.00',
     .'NGC5286     13 46 26.5 -51 22 24.00',
     .'AM4         13 55 50.1 -27 10 22.00',
     .'NGC5466     14 05 27.3 +28 32 04.00',
     .'NGC5634     14 29 37.3 -05 58 35.00',
     .'NGC5694     14 39 36.5 -26 32 18.00',
     .'IC4499      15 00 18.5 -82 12 49.00',
     .'NGC5824     15 03 58.5 -33 04 04.00',
     .'Pal5        15 16 05.3 -00 06 41.00',
     .'NGC5897     15 17 24.5 -21 00 37.00',
     .'NGC5904 M05 15 18 33.8 +02 04 58.00',
     .'NGC5927     15 28 00.5 -50 40 22.00',
     .'NGC5946     15 35 28.5 -50 39 34.00',
     .'BH176       15 39 07.3 -50 03 02.00',
     .'NGC5986     15 46 03.5 -37 47 10.00',
     .'Lynga7      16 11 03.0 -55 18 52.00',
     .'Pal14       16 11 04.9 +14 57 29.00',
     .'NGC6093 M80 16 17 02.5 -22 58 30.00',
     .'NGC6121 M04 16 23 35.5 -26 31 31.00',
     .'NGC6101     16 25 48.6 -72 12 06.00',
     .'NGC6144     16 27 14.1 -26 01 29.00',
     .'NGC6139     16 27 40.4 -38 50 56.00',
     .'Terzan 3    16 28 40.1 -35 21 13.00',
     .'NGC6171     16 32 31.9 -13 03 13.00',
     .'1636-283    16 39 25.5 -28 23 52.00',
     .'NGC6205 M13 16 41 41.5 +36 27 37.00',
     .'NGC6229     16 46 58.9 +47 31 40.00',
     .'NGC6218 M12 16 47 14.5 -01 56 52.00',
     .'NGC6235     16 53 25.4 -22 10 38.00',
     .'NGC6254 M10 16 57 08.9 -04 05 58.00',
     .'NGC6256     16 59 32.6 -37 07 17.00',
     .'Pal15       17 00 02.4 -00 32 31.00',
     .'NGC6266 M62 17 01 12.8 -30 06 49.00',
     .'NGC6273 M19 17 02 37.8 -26 16 05.00',
     .'NGC6284     17 04 28.8 -24 45 53.00',
     .'NGC6287     17 05 09.4 -22 42 29.00',
     .'NGC6293     17 10 10.2 -26 34 55.00',
     .'NGC6304     17 14 32.1 -29 27 44.00',
     .'NGC6316     17 16 37.3 -28 08 24.00',
     .'NGC6341 M92 17 17 07.3 +43 08 11.00',
     .'NGC6325     17 17 59.2 -23 45 57.00',
     .'NGC6333 M09 17 19 11.8 -18 30 59.00',
     .'NGC6342     17 21 10.2 -19 35 14.00',
     .'NGC6356     17 23 35.0 -17 48 47.00',
     .'NGC6355     17 23 58.6 -26 21 13.00',
     .'NGC6352     17 25 29.2 -48 25 22.00',
     .'IC1257      17 27 08.5 -07 05 35.00',
     .'Terzan2 HP3 17 27 33.1 -30 48 08.00',
     .'NGC6366     17 27 44.3 -05 04 36.00',
     .'Terzan4 HP4 17 30 39.0 -31 35 44.00',
     .'HP1         17 31 05.2 -29 58 54.00',
     .'NGC6362     17 31 54.8 -67 02 53.00',
     .'Liller1     17 33 24.5 -33 23 20.00',
     .'NGC6380 Ton 17 34 28.0 -39 04 09.00',
     .'Terzan1     17 35 47.2 -30 28 54.00',
     .'Ton2        17 36 10.5 -38 33 12.00',
     .'NGC6388     17 36 17.0 -44 44 06.00',
     .'NGC6402 M14 17 37 36.1 -03 14 45.00',
     .'NGC6401     17 38 36.6 -23 54 34.00',
     .'NGC6397     17 40 41.3 -53 40 25.00',
     .'Pal6        17 43 42.2 -26 13 21.00',
     .'NGC6426     17 44 54.7 +03 10 13.00',
     .'Djorg1      17 47 28.3 -33 03 56.00',
     .'Terzan5 Trz 17 48 04.9 -24 46 45.00',
     .'NGC6440     17 48 52.7 -20 21 37.00',
     .'NGC6441     17 50 12.9 -37 03 05.00',
     .'Terzan6 HP5 17 50 46.4 -31 16 31.00',
     .'NGC6453     17 50 51.7 -34 35 57.00',
     .'UKS1        17 54 27.2 -24 08 43.00',
     .'NGC6496     17 59 02.0 -44 15 54.00',
     .'Terzan9     18 01 38.8 -26 50 23.00',
     .'Djorg2      18 01 49.1 -27 49 33.00',
     .'NGC6517     18 01 50.6 -08 57 32.00',
     .'Terzan10    18 02 57.4 -26 04 00.00',
     .'NGC6522     18 03 34.1 -30 02 02.00',
     .'NGC6535     18 03 50.7 -00 17 49.00',
     .'NGC6528     18 04 49.6 -30 03 21.00',
     .'NGC6539     18 04 49.8 -07 35 09.00',
     .'NGC6540 Dj  18 06 08.6 -27 45 55.00',
     .'NGC6544     18 07 20.6 -24 59 51.00',
     .'NGC6541     18 08 02.2 -43 30 00.00',
     .'2MSGC01     18 08 21.8 -19 49 47.00',
     .'ESOSC06     18 09 06.0 -46 25 23.00',
     .'NGC6553     18 09 17.6 -25 54 31.00',
     .'2MSGC02     18 09 36.5 -20 46 44.00',
     .'NGC6558     18 10 17.6 -31 45 50.00',
     .'IC1276  Pal 18 10 44.2 -07 12 27.00',
     .'Terzan12    18 12 15.8 -22 44 31.00',
     .'NGC6569     18 13 38.8 -31 49 37.00',
     .'NGC6584     18 18 37.7 -52 12 54.00',
     .'NGC6624     18 23 40.5 -30 21 40.00',
     .'NGC6626 M28 18 24 32.9 -24 52 12.00',
     .'NGC6638     18 30 56.1 -25 29 51.00',
     .'NGC6637 M69 18 31 23.2 -32 20 53.00',
     .'NGC6642     18 31 54.1 -23 28 31.00',
     .'NGC6652     18 35 45.7 -32 59 25.00',
     .'NGC6656 M22 18 36 24.2 -23 54 12.00',
     .'Pal8        18 41 29.9 -19 49 33.00',
     .'NGC6681 M70 18 43 12.7 -32 17 31.00',
     .'NGC6712     18 53 04.3 -08 42 22.00',
     .'NGC6715 M54 18 55 03.3 -30 28 42.00',
     .'NGC6717 Pal 18 55 06.2 -22 42 03.00',
     .'NGC6723     18 59 33.2 -36 37 54.00',
     .'NGC6749     19 05 15.3 +01 54 03.00',
     .'NGC6752     19 10 52.0 -59 59 05.00',
     .'NGC6760     19 11 12.1 +01 01 50.00',
     .'NGC6779 M56 19 16 35.5 +30 11 05.00',
     .'Terzan7     19 17 43.7 -34 39 27.00',
     .'Pal10       19 18 02.1 +18 34 18.00',
     .'Arp2        19 28 44.1 -30 21 14.00',
     .'NGC6809 M55 19 39 59.4 -30 57 44.00',
     .'Terzan8     19 41 45.0 -34 00 01.00',
     .'Pal11       19 45 14.4 -08 00 26.00',
     .'NGC6838 M71 19 53 46.1 +18 46 42.00',
     .'NGC6864 M75 20 06 04.8 -21 55 17.00',
     .'NGC6934     20 34 11.6 +07 24 15.00',
     .'NGC6981 M72 20 53 27.9 -12 32 13.00',
     .'NGC7006     21 01 29.5 +16 11 15.00',
     .'NGC7078 M15 21 29 58.3 +12 10 01.00',
     .'NGC7089 M02 21 33 29.3 -00 49 23.00',
     .'NGC7099 M30 21 40 22.0 -23 10 45.00',
     .'Pal12       21 46 38.8 -21 15 03.00',
     .'Pal13       23 06 44.4 +12 46 19.00',
     .'NGC7492     23 08 26.7 -15 36 41.00',
     .'ANDROMEDM31 00 42 44.3 +41 16 09.00',
     .'LMC         05 23 34.6 -69 45 22.00', ! for now replace LMC with JWST-CAL
     .'JWST-CALIB  05 21 57.0 -69 29 54.00',
     .'SMC         00 52 42.0 -72 49 00.00',
     .'BPTAU       04 19 15.8 +29 06 27.00',
     .'DFTAU       04 27 02.8 +25 42 22.00',
     .'ORION-PAR   05 35 25.0 -05 30 59.00',
     .'UDF         03 32 39.0 -27 47 29.00',
     .'JWST-CALIB  05 21 57.0 -69 29 54.00',
     .'BULGE-SWEEP 17 58 59.1 -29 12 17.73',
     .'BULGE-STANK 17 54 41.7 -29 49 30.00',
     .'BULGE-OGL29 17 48 15.3 -37 09 01.40',
     .'BULGE-BAADE 18 03 10.2 -29 56 33.60',
     .'TBROWN-DISK 18 59 45.0 -04 26 08.40',
     .'NGC6791     19 20 53.0 +37 46 30.00', ! <---
     .'SDSSJ090744 09 07 44.0 +02 45 06.40',
     .'SDSSJ093320 09 33 20.9 +44 17 05.50',
     .'HE0437-5439 04 38 12.8 -53 33 11.80',
     .'SDSSJ091301 09 13 01.0 +30 51 20.20',
     .'SDSSJ091759 09 17 59.6 +67 22 38.80',
     .'SDSSJ110557 11 05 57.5 +09 34 39.47',
     .'SDSSJ113312 11 33 12.1 +01 08 24.87',
     .'SDSSJ094214 09 42 14.1 +20 03 22.07',
     .'SDSSJ102137 10 21 37.1 -00 52 34.77',
     .'SDSSJ120337 12 03 37.6 +18 02 50.35',
     .'SDSSJ105009 10 50 09.6 +03 15 50.67',
     .'SDSSJ105248 10 52 48.3 -00 01 33.94',
     .'SDSSJ144955 14 49 55.6 +31 03 51.37',
     .'SDSSJ101018 10 10 18.8 +30 20 28.22',
     .'SDSSJ110224 11 02 24.4 +02 50 02.77',
     .'MO9-BLG-195 18 03 01.2 -28 23 26.50',
     .'MO9-BLG-046 17 54 43.4 -35 03 07.50',
     .'MO9-BLG-260 17 58 28.6 -26 50 21.10',
     .'OGL9-BULG-D 17 30 0.00 -29 00 00.00',
     .'OGL9-BULG-E 17 30 0.00 -30 00 00.00',
     .'SDSSJ142001 14 20 01.9 +12 44 04.70',
     .'PLEIADES_01 03 47 17.0 +24 14 50.60',
     .'PLEIADES_02 03 48 12.2 +23 59 07.30',
     .'PLEIADES_03 03 48 37.6 +24 13 11.40',
     .'30__DORADUS 05 38 40.0 -69 06 00.00',
     .'GALACTICCEN 17 45 40.0 -29 00 28.10'/
 
 
      dmin = 30.0       ! 30 arcmin finding radius
      NAME = 'NO ID'
      rRA0 = rRA
      rDC0 = rDC
      do N = 1, 191
         RSTRING = INFO(N)(13:22)
         DSTRING = INFO(N)(24:32)
         read(RSTRING,*) iRAH,iRAM,rRAS
         read(DSTRING,*) iDEG,iMIN,iSEC
         rRAn = (iRAH+iRAM/60.+rRAS/60./60.)*360/24  ! in degs
         rDCn = abs(iDEG)+iMIN/60.+iSEC/60./60.      ! in degs
c        if (iDEG.lt.0) rDCn = -rDCn
         if (DSTRING(01:01).eq.'-') rDCn = -rDCn
         dRA = (rRA-rRAn)*60.0*cos(rDC*3.14159/180)  ! in arcmin (W+)
         dDC = (rDC-rDCn)*60.0                       ! in arcmin
         dDD = sqrt(dRA**2+dDC**2)
         if (dDD.lt.dmin) then
            dmin = dDD
            rRA0 = rRAn
            rDC0 = rDCn
            NAME = INFO(N)(01:11)
            endif
         enddo
 
 
      return
      end
 
 
