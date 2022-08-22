
      program imginfo
      implicit none

      integer i
      character filename*80
      integer NCHAR

      integer iargc, nargc
      integer q

      integer iEXPT
      real    rEXPT
      integer iPA
      real    rPA
      integer iPID
      character*6 FILTSTR
      character*12 DATESTRIN
      character*9 ROOTSTR
      real rRA, hRA
      real rDEC, aDEC

      character*8 DATESTR 
      character*8 TIMESTR

      integer iRAH, iRAM, iRAS, iFRA
      integer iDEG, iMIN, iSEC, iFDC
      character decsgn
      real    rRA0, rDC0
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


      character*11 OBJNAME

      integer DOY(12)
      data DOY / 000, 031, 059, 090, 120, 151, 
     .           181, 212, 243, 273, 304, 334/ 


      character*8 CHINJVAL
      integer NAXIS, NAXIS1, NAXIS2, NAXIS3
      character*12 DATATYPE
      integer BITPIX
      integer iii, iiii

      integer NEXTEND, NEXTENDs
      integer NBPER

      if (iargc().eq.0) then
         print*,' '
         print*,'TO USE THIS ROUTINE, YOU NEED TO'
         print*,'GIVE IT SOME FITS IMAGES ON THE'
         print*,'COMMAND LINE, USING ASTERISKS IF'
         print*,'DESIRED.'
         print*,' '
         print*,'FOR EXAMPLE, '
         print*,' '
         print*,'% imginfo.e *.fits '
         print*,' '
         print*,'WILL TABULATE THE INFO FOR ALL THE '
         print*,'FITS FILES IN THE CURRENT DIRECTORY. '
         print*,' '
         endif


      write(*,110) 
      write(*,111) 
      write(*,110) 
 110  format('#----      ---  -----     -- x   -- x -- ------',
     .       ' --------     ----')
 111  format('#NEXT      REC  NAXES     NX x   NY x NZ BITPIX',
     .       ' DATATYPE     FILE')

      do nargc = 1, iargc()

         CHINJVAL = '        '
         DATESTR  = 'DATELESS'
         TIMESTR  = 'TIMELESS'
         ROOTSTR  = 'ROOTLESS'
         FILTSTR  = 'FLTLSS'
         rEXPT    = 0
         iEXPT    = 0
         rRA  = 0.
         rDEC = 0.
         iPID = 0
         iPA  = 0
         rPA  = 0
         rSUNA = 0.
         rTEMP = 0.
         vaFAC = 0.
         call getarg(nargc,filename)

c-----------------------------------------------

         open(10,file=FILENAME,status='old',
     .     err=900,recl=2880,form='UNFORMATTED',
     .     access='DIRECT')

         NAXIS  = 0
         NAXIS1 = 0
         NAXIS2 = 1
         NAXIS3 = 1
         DATATYPE = 'NONE'


         NEXTEND  = 0
         NEXTENDs = 0

         i = 0
         nread = 0
 099     continue
         NAXIS = 0
         NAXIS1 = 0
         NAXIS2 = 1
         NAXIS3 = 1
 100     continue
         i = i + 1
         read(10,rec=i,iostat=ios) buff
         if (ios.eq.5002) goto 003
         if (ios.lt.0000) goto 900
         do k = 0, 35, 1
            field  = buff(k*80+01:k*80+08)
            stream = buff(k*80+11:k*80+31)
            if (stream(1:1).eq.'''') stream    = stream(2:20)
            if (field.eq.'CHINJECT') CHINJVAL = stream(1:8)
            if (field.eq.'BITPIX  ') read(stream,*) BITPIX
            if (field.eq.'NAXIS   ') read(stream,*) NAXIS
            if (field.eq.'NAXIS1  ') read(stream,*) NAXIS1
            if (field.eq.'NAXIS2  ') read(stream,*) NAXIS2
            if (field.eq.'NAXIS3  ') read(stream,*) NAXIS3
            if (field.eq.'NEXTEND ') read(stream,*) NEXTENDs
            if (field.eq.'DATATYPE') then
               DATATYPE = '          '
               iiii = 0
               do iii = 1,  12
                  if (stream(iii:iii).ne.' '.and.
     .                stream(iii:iii).ne.'''') then
                      iiii = iiii + 1
                      DATATYPE(iiii:iiii) = stream(iii:iii)
                      endif
                  enddo
               endif
c           write(*,'(i2.2,1x,i2.2,'':'',a8,''---'',a20,5x,a6)') 
c    .               i,k,field,stream,FILTSTR
            if (field.eq.'END     ') goto 101
 109        continue
            enddo 
         goto 100
 101     continue
         
         if (BITPIX.eq.  8) DATATYPE = 'CHAR*1      '
         if (BITPIX.eq. 16) DATATYPE = 'INTEGER*2   '
         if (BITPIX.eq. 32) DATATYPE = 'INTEGER*4   '
         if (BITPIX.eq.-32) DATATYPE = 'REAL*4      '
         if (BITPIX.eq.-64) DATATYPE = 'REAL*8      '

         write(*,112) NEXTEND, i, NAXIS, NAXIS1, NAXIS2, NAXIS3,
     .                            BITPIX, DATATYPE,
     .                            FILENAME(1:80)
 112     format(1x,i4,1x,i8,1x,i6,1x,i6,'x',i6,'x',i4,
     .          3x,i4,1x,a12,
     .          1x,80a)         

         NBPER = abs(BITPIX)/8*NAXIS1*NAXIS2*NAXIS3/2880

         NEXTEND = NEXTEND + 1
         i = i + NBPER

         goto 099

c        if (NEXTEND.lt.NEXTENDs) then
c           NEXTEND = NEXTEND + 1
c           goto 099
c           endif

         print*,'CLOSE(10)' 
    3    close(10)

         enddo!NARGs
      stop
 
  900 continue
      print*,' '
      print*,'probe_fits.e ERROR EXIT. '
      print*,' '
      print*,'ONE OF THE IMAGES WAS NOT IN STANDARD'
      print*,'HST FITS FORMAT.'
      print*,' '
      write(*,'(''PROBLEM FILE: '',a80)') FILENAME
      print*,' '
      end


c----------------------------------------------------
c
c find the closest object of interest to this
c location; info from HARRIS' on-line catalog
c is contained within data statements.
c
      subroutine find_obj(rRA,rDC,NAME,rRA0,rDC0)
      implicit none
      real rRA , rDC      ! input
      character*11 NAME   ! output (nearest obj)
      real rRA0, rDC0     ! output (locn)
  
      real rRAn, rDCn
      real dRA,  dDC
      real dDD, dmin
      integer N
      character*10 RSTRING
      character*09 DSTRING

      integer iRAH, iRAM
      real    rRAS
      integer iDEG, iMIN, iSEC

      character*32 INFO(163)
      common /INFO/INFO ! need to do this in LINUX to preserve
      data INFO / 
     .'NGC0104 47T 00 24 05.2 -72 04 51',
     .'NGC0288     00 52 47.5 -26 35 24',
     .'NGC0362     01 03 14.3 -70 50 54',
     .'NGC1261     03 12 15.3 -55 13 01',
     .'Pal1        03 33 23.0 +79 34 50',
     .'AM1     E01 03 55 02.7 -49 36 52',
     .'Eridanus    04 24 44.5 -21 11 13',
     .'Pal2        04 46 05.9 +31 22 51',
     .'NGC1851     05 14 06.3 -40 02 50',
     .'NGC1904 M79 05 24 10.6 -24 31 27',
     .'NGC2298     06 48 59.2 -36 00 19',
     .'NGC2419     07 38 08.5 +38 52 55',
     .'Pyxis       09 07 57.8 -37 13 17',
     .'NGC2808     09 12 02.6 -64 51 47',
     .'E3          09 20 59.3 -77 16 57',
     .'Pal3        10 05 31.4 +00 04 17',
     .'NGC3201     10 17 36.8 -46 24 40',
     .'Pal4        11 29 16.8 +28 58 25',
     .'NGC4147     12 10 06.2 +18 32 31',
     .'NGC4372     12 25 45.4 -72 39 33',
     .'Rup106      12 38 40.2 -51 09 01',
     .'NGC4590 M68 12 39 28.0 -26 44 34',
     .'NGC4833     12 59 35.0 -70 52 29',
     .'NGC5024 M53 13 12 55.3 +18 10 09',
     .'NGC5053     13 16 27.0 +17 41 53',
     .'NGC5139 OMC 13 26 45.9 -47 28 37',
     .'NGC5272 M03 13 42 11.2 +28 22 32',
     .'NGC5286     13 46 26.5 -51 22 24',
     .'AM4         13 55 50.1 -27 10 22',
     .'NGC5466     14 05 27.3 +28 32 04',
     .'NGC5634     14 29 37.3 -05 58 35',
     .'NGC5694     14 39 36.5 -26 32 18',
     .'IC4499      15 00 18.5 -82 12 49',
     .'NGC5824     15 03 58.5 -33 04 04',
     .'Pal5        15 16 05.3 -00 06 41',
     .'NGC5897     15 17 24.5 -21 00 37',
     .'NGC5904 M05 15 18 33.8 +02 04 58',
     .'NGC5927     15 28 00.5 -50 40 22',
     .'NGC5946     15 35 28.5 -50 39 34',
     .'BH176       15 39 07.3 -50 03 02',
     .'NGC5986     15 46 03.5 -37 47 10',
     .'Lynga7      16 11 03.0 -55 18 52',
     .'Pal14       16 11 04.9 +14 57 29',
     .'NGC6093 M80 16 17 02.5 -22 58 30',
     .'NGC6121 M04 16 23 35.5 -26 31 31',
     .'NGC6101     16 25 48.6 -72 12 06',
     .'NGC6144     16 27 14.1 -26 01 29',
     .'NGC6139     16 27 40.4 -38 50 56',
     .'Terzan 3    16 28 40.1 -35 21 13',
     .'NGC6171     16 32 31.9 -13 03 13',
     .'1636-283    16 39 25.5 -28 23 52',
     .'NGC6205 M13 16 41 41.5 +36 27 37',
     .'NGC6229     16 46 58.9 +47 31 40',
     .'NGC6218 M12 16 47 14.5 -01 56 52',
     .'NGC6235     16 53 25.4 -22 10 38',
     .'NGC6254 M10 16 57 08.9 -04 05 58',
     .'NGC6256     16 59 32.6 -37 07 17',
     .'Pal15       17 00 02.4 -00 32 31',
     .'NGC6266 M62 17 01 12.8 -30 06 49',
     .'NGC6273 M19 17 02 37.8 -26 16 05',
     .'NGC6284     17 04 28.8 -24 45 53',
     .'NGC6287     17 05 09.4 -22 42 29',
     .'NGC6293     17 10 10.2 -26 34 55',
     .'NGC6304     17 14 32.1 -29 27 44',
     .'NGC6316     17 16 37.3 -28 08 24',
     .'NGC6341 M92 17 17 07.3 +43 08 11',
     .'NGC6325     17 17 59.2 -23 45 57',
     .'NGC6333 M09 17 19 11.8 -18 30 59',
     .'NGC6342     17 21 10.2 -19 35 14',
     .'NGC6356     17 23 35.0 -17 48 47',
     .'NGC6355     17 23 58.6 -26 21 13',
     .'NGC6352     17 25 29.2 -48 25 22',
     .'IC1257      17 27 08.5 -07 05 35',
     .'Terzan2 HP3 17 27 33.1 -30 48 08',
     .'NGC6366     17 27 44.3 -05 04 36',
     .'Terzan4 HP4 17 30 39.0 -31 35 44',
     .'HP1         17 31 05.2 -29 58 54',
     .'NGC6362     17 31 54.8 -67 02 53',
     .'Liller1     17 33 24.5 -33 23 20',
     .'NGC6380 Ton 17 34 28.0 -39 04 09',
     .'Terzan1     17 35 47.2 -30 28 54',
     .'Ton2        17 36 10.5 -38 33 12',
     .'NGC6388     17 36 17.0 -44 44 06',
     .'NGC6402 M14 17 37 36.1 -03 14 45',
     .'NGC6401     17 38 36.6 -23 54 34',
     .'NGC6397     17 40 41.3 -53 40 25',
     .'Pal6        17 43 42.2 -26 13 21',
     .'NGC6426     17 44 54.7 +03 10 13',
     .'Djorg1      17 47 28.3 -33 03 56',
     .'Terzan5 Trz 17 48 04.9 -24 46 45',
     .'NGC6440     17 48 52.7 -20 21 37',
     .'NGC6441     17 50 12.9 -37 03 05',
     .'Terzan6 HP5 17 50 46.4 -31 16 31',
     .'NGC6453     17 50 51.7 -34 35 57',
     .'UKS1        17 54 27.2 -24 08 43',
     .'NGC6496     17 59 02.0 -44 15 54',
     .'Terzan9     18 01 38.8 -26 50 23',
     .'Djorg2      18 01 49.1 -27 49 33',
     .'NGC6517     18 01 50.6 -08 57 32',
     .'Terzan10    18 02 57.4 -26 04 00',
     .'NGC6522     18 03 34.1 -30 02 02',
     .'NGC6535     18 03 50.7 -00 17 49',
     .'NGC6528     18 04 49.6 -30 03 21',
     .'NGC6539     18 04 49.8 -07 35 09',
     .'NGC6540 Dj  18 06 08.6 -27 45 55',
     .'NGC6544     18 07 20.6 -24 59 51',
     .'NGC6541     18 08 02.2 -43 30 00',
     .'2MSGC01     18 08 21.8 -19 49 47',
     .'ESOSC06     18 09 06.0 -46 25 23',
     .'NGC6553     18 09 17.6 -25 54 31',
     .'2MSGC02     18 09 36.5 -20 46 44',
     .'NGC6558     18 10 17.6 -31 45 50',
     .'IC1276  Pal 18 10 44.2 -07 12 27',
     .'Terzan12    18 12 15.8 -22 44 31',
     .'NGC6569     18 13 38.8 -31 49 37',
     .'NGC6584     18 18 37.7 -52 12 54',
     .'NGC6624     18 23 40.5 -30 21 40',
     .'NGC6626 M28 18 24 32.9 -24 52 12',
     .'NGC6638     18 30 56.1 -25 29 51',
     .'NGC6637 M69 18 31 23.2 -32 20 53',
     .'NGC6642     18 31 54.1 -23 28 31',
     .'NGC6652     18 35 45.7 -32 59 25',
     .'NGC6656 M22 18 36 24.2 -23 54 12',
     .'Pal8        18 41 29.9 -19 49 33',
     .'NGC6681 M70 18 43 12.7 -32 17 31',
     .'NGC6712     18 53 04.3 -08 42 22',
     .'NGC6715 M54 18 55 03.3 -30 28 42',
     .'NGC6717 Pal 18 55 06.2 -22 42 03',
     .'NGC6723     18 59 33.2 -36 37 54',
     .'NGC6749     19 05 15.3 +01 54 03',
     .'NGC6752     19 10 52.0 -59 59 05',
     .'NGC6760     19 11 12.1 +01 01 50',
     .'NGC6779 M56 19 16 35.5 +30 11 05',
     .'Terzan7     19 17 43.7 -34 39 27',
     .'Pal10       19 18 02.1 +18 34 18',
     .'Arp2        19 28 44.1 -30 21 14',
     .'NGC6809 M55 19 39 59.4 -30 57 44',
     .'Terzan8     19 41 45.0 -34 00 01',
     .'Pal11       19 45 14.4 -08 00 26',
     .'NGC6838 M71 19 53 46.1 +18 46 42',
     .'NGC6864 M75 20 06 04.8 -21 55 17',
     .'NGC6934     20 34 11.6 +07 24 15',
     .'NGC6981 M72 20 53 27.9 -12 32 13',
     .'NGC7006     21 01 29.5 +16 11 15',
     .'NGC7078 M15 21 29 58.3 +12 10 01',
     .'NGC7089 M02 21 33 29.3 -00 49 23',
     .'NGC7099 M30 21 40 22.0 -23 10 45',
     .'Pal12       21 46 38.8 -21 15 03',
     .'Pal13       23 06 44.4 +12 46 19',
     .'NGC7492     23 08 26.7 -15 36 41',
     .'PLEIADESM45 03 47 00.0 +24 07 00',
     .'ANDROMEDM31 00 42 44.3 +41 16 09',
c    .'LMC         05 23 34.6 -69 45 22', ! for now replace LMC with JWST-CAL
     .'JWST-CALIB  05 21 57.0 -69 29 54',
     .'SMC         00 52 42.0 -72 49 00',
     .'BPTAU       04 19 15.8 +29 06 27',
     .'DFTAU       04 27 02.8 +25 42 22',
     .'ORION-PAR   05 35 25.0 -05 30 59',
     .'UDF         03 32 39.0 -27 47 29',
     .'JWST-CALIB  05 21 57.0 -69 29 54',
     .'ARCHES      17 45 50.5 -28 49 28',
     .'QUINTUPLET  17 46 13.9 -28 49 48',
     .'SAG-A*      17 45 40.0 -29 00 28',
     .'30-DOR      05 38 38.0 -69 05 40'/

      dmin = 120.0
      NAME = 'NO ID'
      rRA0 = rRA
      rDC0 = rDC
      do N = 1, 163
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

