

      program catheader
      implicit none

C     Print out all the header keywords in all extensions of a FITS file

      integer status,readwrite,blocksize,nkeys,nspace,i
      character filename*80,record*80

      integer iargc, nargc

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
      real    rRA0, rDC0
      real    dRA , dDC
      real    rSUNA
      real    rTEMP

      integer iDAY
      integer iMON
      integer iYIR
      integer iHRS

      integer iq
      integer nimg
      integer naxes
      integer laxis(3)

      character*8  field
      character*20 stream 
      
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer ios, k
      integer j

      character*2880 buff
      common /sneaky/ buff

      integer ifirst, i1, i2

      integer np1, np2, npt
      integer NEXTEND
      logical EXTEND
      integer nread
      real bscale, bzero
      integer FLEN
      logical DIAG
      data DIAG/.false./

      integer BSKIP
      integer ISKIP
      integer r2i
      real    i2r

      integer b

      integer DOY(12)
      data DOY / 000, 031, 059, 090, 120, 151, 
     .           181, 212, 243, 273, 304, 303/ 

      integer BITPIX
      integer NAXIS
      integer NAXIS1
      integer NAXIS2
      integer NAXIS3

      logical BYTEDUMP

      if (iargc().lt.1) stop 'need at least one fits image as agrument'

      BYTEDUMP = .false.

      do nargc = 1, iargc()


         rRA0 = 201.6910 ! OMEGA CEN 
         rDC0 = -47.4769 ! OMEGA CEN 

         DATESTR  = 'DATELESS'
         TIMESTR  = 'TIMELESS'
         ROOTSTR  = 'ROOTLESS'
         FILTSTR  = 'FILTLESS'
         rEXPT    = 0
         iEXPT    = 0
         rRA  = 0.
         rDEC = 0.
         iPID = 0
         iPA  = 0
         rPA  = 0
         rSUNA = 0.
         rTEMP = 0.
         call getarg(nargc,FILENAME)

         if (FILENAME(1:4).eq.'BYTEDUMP') then
             BYTEDUMP = .true.
             goto 9999
             endif
 
         FLEN = 80
         do i = 1, 80
            if (FILENAME(i:i).eq.' ') then
               FLEN = i-1
               goto 5
               endif
            enddo
    5    continue 
c
c-----------------------------------------------
c
         print*,' '
         write(*,'(''FILE: '',80a)') FILENAME

         open(10,file=FILENAME,status='old',
     .        err=900,recl=2880,form='UNFORMATTED',
     .        access='DIRECT')

         EXTEND  = .false.
         NEXTEND = 0
         NREAD   = 0
         i = 0
 100     continue
         BITPIX = 0
         NAXIS  = 0
         NAXIS1 = 0
         NAXIS2 = 0
         NAXIS3 = 0
         write(*,'(10x,30(''*''),''  EXTENSION: '',i2,3x,32(''*''))') 
     .                               NREAD 
 102     continue
         i = i + 1
         read(10,rec=i,iostat=ios) buff
         if (ios.ne.0) goto 9999
         do b = 0001, 2880                    ! purge of junk characters...
            if (BYTEDUMP) then
                do k = 0, 35
                   enddo
                endif
            if (buff(b:b).ne.'.'.and.
     .          buff(b:b).ne.' '.and.
     .          buff(b:b).ne.'+'.and.
     .          buff(b:b).ne.'_'.and.
     .          buff(b:b).ne.'-'.and.
     .          buff(b:b).ne.'/'.and.
     .          buff(b:b).ne.'='.and.
     .          buff(b:b).ne.'*'.and.
     .          buff(b:b).ne.'?'.and.
     .          buff(b:b).ne.'!'.and.
     .          buff(b:b).ne.'\\'.and.
     .          buff(b:b).ne.''''.and.
     .          buff(b:b).ne.':'.and.
     .          buff(b:b).ne.'@'.and.
     .          buff(b:b).ne.'$'.and.
     .          buff(b:b).ne.'^'.and.
     .          buff(b:b).ne.'&'.and.
     .          buff(b:b).ne.'('.and.
     .          buff(b:b).ne.')'.and.
     .          buff(b:b).ne.'['.and.
     .          buff(b:b).ne.']'.and.
     .          buff(b:b).ne.'{'.and.
     .          buff(b:b).ne.'}'.and.
     .          buff(b:b).ne.'~'.and.
     .          buff(b:b).ne.'#'.and.
     .          buff(b:b).ne.';'.and.
     .          buff(b:b).ne.'>'.and.
     .          buff(b:b).ne.'<'.and.
     .          buff(b:b).ne.'|'.and.
     .          buff(b:b).ne.','.and.
     .          (.not.(buff(b:b).ge.'A'.and.buff(b:b).le.'Z')).and.
     .          (.not.(buff(b:b).ge.'a'.and.buff(b:b).le.'z')).and.
     .          (.not.(buff(b:b).ge.'0'.and.buff(b:b).le.'9'))) then
                print*,'FIX: ',b,' ',iachar(buff(b:b)),buff(b:b)
                buff(b:b) = ' '
                endif
            enddo
         do k = 0, 35, 1
            write(*,122) i,k,NREAD,buff(k*80+1:k*80+80),
     .                   NARGC,FILENAME(1:FLEN)
  122       format('P'i6.6,1x,'L'i2.2,1x,'E',i2.2,1x,a80,1x,i5.5,1x,80a)
            field  = buff(k*80+01:k*80+08)
            stream = buff(k*80+11:k*80+31)
            if (field.eq.'EXTEND  ') read(stream,*) EXTEND
            if (field.eq.'NEXTEND ') read(stream,*) NEXTEND
            if (field.eq.'BITPIX  ') read(stream,*) BITPIX
            if (field.eq.'NAXIS   ') read(stream,*) NAXIS
            if (field.eq.'NAXIS1  ') read(stream,*) NAXIS1
            if (field.eq.'NAXIS2  ') read(stream,*) NAXIS2
            if (field.eq.'NAXIS3  ') read(stream,*) NAXIS3
            if (field.eq.'END     ') goto 101
            enddo 
         write(*,'(10x,80(''-''))') 
         goto 102
 101     continue
         if (DIAG) then
            print*,'          '
            print*,'  EXTEND: ',EXTEND
            print*,' NEXTEND: ',NEXTEND
            print*,'  NAXIS : ',NAXIS
            print*,'  NAXIS1: ',NAXIS1
            print*,'  NAXIS2: ',NAXIS2
            print*,'  NAXIS3: ',NAXIS3
            print*,'  BITPIX: ',BITPIX,abs(BITPIX)/8
            endif
         BSKIP = 0
         if (NAXIS.eq.1) BSKIP = abs(BITPIX)/8*NAXIS1
         if (NAXIS.eq.2) BSKIP = abs(BITPIX)/8*NAXIS1*NAXIS2
         if (NAXIS.eq.3) BSKIP = abs(BITPIX)/8*NAXIS1*NAXIS2*NAXIS3
         ISKIP = (BSKIP+2879)/2880
         NREAD = NREAD + 1
         i = i + ISKIP
         read(10,rec=i+1,iostat=ios) buff
         if (ios.ne.0) goto 9999
         if (EXTEND) goto 100
 9999    continue

         enddo

      stop 'END OF LIST'

  900 continue
      stop 'ERROR OPENING FILE'

      end


