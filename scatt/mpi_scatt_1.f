!====================================================================
	subroutine scatt_start_mpi
	use decdatamd
	use decarraymd
	implicit none

	include 'mpif.h'

	integer  :: ierr
	integer  :: provided, namelen
      character*(MPI_MAX_PROCESSOR_NAME) :: processor_name

	!!call mpi_init(ierr)
	call mpi_init_thread(mpi_thread_multiple, provided, ierr)

	call mpi_comm_rank(mpi_comm_world, tasknum, ierr)
	call mpi_comm_size(mpi_comm_world, nproc, ierr)
      call MPI_GET_PROCESSOR_NAME(processor_name, namelen, ierr) !!

	if(tasknum == 0) then !!
	  print*,'provided=',provided
	  print*,'mpi_thread_multiple=',mpi_thread_multiple 
	endif   

	if(tasknum == 0) write(9000,*) 'nproc=',nproc

      print*,'node ',taskNum,' @ ',trim(processor_name)

      if(kAddr /= MPI_ADDRESS_KIND) then !!
        if(taskNum == 0) then
          print*,'Err: kAddr /= MPI_ADDRESS_KIND'
          print*,'kAddr=',kAddr,' MPI_ADDRESS_KIND=',MPI_ADDRESS_KIND
          print*,'Value of kAddr should be larger, modify the code'   
        end if
        stop
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

	return
	end subroutine scatt_start_mpi
!====================================================================
	subroutine rowdecomp(n5d)
	use decdatamd
	use decarraymd 
	implicit none

	include 'mpif.h'

	integer  :: n5d
	integer  :: i, ierr 

	nrow=n5d
	!fullbufsz  = 2*kcol*int(ceiling(dble(nrow)/dble(nproc)))
	fullbufsz  = kcol*int(ceiling(dble(nrow)/dble(nproc)))
	blockcolsz = kcol*nproc

	do i=0, nproc-1
	   rowdec(1,i)=i*nrow/nproc + 1
	   rowdec(2,i)=(i+1)*nrow/nproc
	   rowdec(3,i)=rowdec(2,i)-rowdec(1,i)+1
	   rowdec(4,i)=rowdec(3,i)*kcol
	   write(9000,*)i,':',rowdec(1:4,i)
	enddo

	NAp=rowdec(3,tasknum)
	print*,'node:',tasknum,': NAp=', NAp

!	allocate(recvbuf(fullbufsz, 1:nproc-1), sendbuf(fullbufsz, 1:nproc-1),
!     &           sendrecvhandle(1:2*(nproc-1)), sendbufsize(0:nproc-1),
!     &           handlestat(MPI_STATUS_SIZE,1:2*(nproc-1)),
!     &           workTmp(nrow,kcol), stat=ierr)
      allocate(rawBuf(fullbufsz,0:nProc-1), 
     &         rstBuf(fullbufsz,0:nProc-1), 
     &         workTmp(nrow,kcol),         stat=ierr)
	if(ierr /= 0) STOP 'Failed: allocate in rowdecomp'

	workTmp(:,:)=(0.0,0.0)
!	recvbuf(:,:)=(0.0,0.0)
!	sendbuf(:,:)=(0.0,0.0) 
	rawBuf(:,:)=(0.0,0.0) 
	rstBuf(:,:)=(0.0,0.0) 

	return
	end subroutine rowdecomp
!====================================================================
!!	subroutine SendRecv(iMpiLoop,nMpiLoop, ColDecLoop,wfMatrix)
!!	use decDatamd
!!	use decArraymd
!!	implicit none
!!	include 'mpif.h'
!!
!!! Dummy Variables
!!	integer :: iMpiLoop, nMpiLoop
!!	integer :: ColDecLoop(3,0:nProc-1,1)
!!	complex(4) :: wfMatrix(1)
!!
!!      print *,'SendRecv: no more available'
!!      stop
!!
!!      end subroutine SendRecv
!====================================================================
      subroutine PreAllExchange()
      use decDatamd
      use nh4scattmd, only: wasy, wint
      use comparamd, only : zzty,r1ty,row1ty,row2ty
      implicit none
      include 'mpif.h'

      ! local variables
      integer(kind=MPI_ADDRESS_KIND) :: nsize1,nsize2
      integer :: szoftype,disp, ierr, ns0

      call MPI_TYPE_SIZE(MPI_COMPLEX8, szoftype, ierr)

      ns0   =1_kAddr*zzty%nzint*r1ty%nvb*row1ty%nvb*row2ty%nvb*nap
!      nsize1=1_kAddr*zzty%nz   *r1ty%nvb*row1ty%nvb*row2ty%nvb*nap*szoftype
      nsize1=1_kAddr*(zzty%nz-zzty%nzint)*r1ty%nvb*row1ty%nvb*row2ty%nvb*nap*szoftype
!      call MPI_Win_Create(wasy,nsize1,szoftype,
      call MPI_Win_Create(wasy(ns0+1),nsize1,szoftype,
     &                    MPI_INFO_NULL,MPI_COMM_WORLD, winAsy, ierr)

      nsize2=1_kAddr*zzty%nzint*r1ty%nrb*row1ty%nvb*row2ty%nvb*nap*szoftype
      call MPI_Win_Create(wint,       nsize2,szoftype,
     &                    MPI_INFO_NULL,MPI_COMM_WORLD, winInt, ierr)

      print*, taskNum, " : winAsy=",winAsy,", winInt=",winInt
      call MPI_Win_Fence(0, winAsy, ierr)
      call MPI_Win_Fence(0, winInt, ierr)
        
      end subroutine PreAllExchange
!====================================================================
      subroutine PostAllExchange
      use decDatamd, only: winAsy,winInt
      implicit none
      include 'mpif.h'

      integer :: ierr

      call MPI_Win_Free(winAsy, ierr)
      call MPI_Win_Free(winInt, ierr)
      end subroutine PostAllExchange
!====================================================================
      subroutine FinishExchange(cPart)
      use decDatamd, only: winAsy,winInt, taskNum
      implicit none
      include 'mpif.h'
      ! dummy variables
      character*3 :: cPart
      ! local variables
      integer :: ierr

!      print*,"FinishExchange-0",taskNum
      select case (cPart)
        case('asy')
          call MPI_Win_Fence(0, winAsy, ierr)
        case('int')
          call MPI_Win_Fence(0, winInt, ierr)
        case default
          print *,'FinishExchange: invalid value for cPart...'
          stop
      end select
!      print*,"FinishExchange-1",taskNum, ierr
      end subroutine finishExchange
!====================================================================
      subroutine getRawData(cPart,iLoop,nLoop,colDecLoop)
      use decDatamd
      use decArraymd
      implicit none
      include 'mpif.h'
      ! dummy variables
      character*3 :: cPart
      integer :: iLoop, nLoop
      integer :: colDecLoop(3, 0:nProc-1, nLoop)

      ! local variables
      integer(kind=MPI_ADDRESS_KIND) :: disp ! disp: remote offset in the window
      integer :: nsize !!
      integer :: iProc, nRowx,iCol0,nColx, ierr

      iCol0=colDecLoop(1,taskNum,iLoop)
      nColx=colDecLoop(3,taskNum,iLoop)
!      print*,taskNum,' :nColx, ',iCol0,nColx,iLoop      
      do iProc=0,nProc-1
        nRowx=rowDec(3,iProc)
        disp =1_kAddr*nRowx*(iCol0-1)
        nsize=nRowx*nColx

        select case (cPart)
          case('asy')
            call MPI_Get(rawBuf(1,iProc),nsize,MPI_COMPLEX8, iProc,
     &                   disp,nsize,MPI_COMPLEX8, winAsy, ierr)
!            print*,taskNum,' dispA=',disp,nsize," from:",iProc
          case('int')
            call MPI_Get(rawBuf(1,iProc),nsize,MPI_COMPLEX8, iProc,
     &                   disp,nsize,MPI_COMPLEX8, winInt, ierr)
!            print*,taskNum,' dispI=',disp,nsize," from:",iProc
          case default
            print *,'getRawData: invalid value for cPart...'
            stop
        end select
      end do      
      end subroutine getRawData
!====================================================================
      subroutine unpackRawData(cPart,iLoop,nLoop,colDecLoop)
      use decDatamd
      use decArraymd
      implicit none
      include 'mpif.h'
      ! dummy variables
      character*3 :: cPart
      integer :: iLoop, nLoop
      integer :: colDecLoop(3, 0:nProc-1, nLoop)

      ! local variables
      integer :: iProc, iRow0,nRowx,iRowx, nColx,iColx, ierr
      integer :: nsize, idata

      nColx=colDecLoop(3,taskNum,iLoop)
      do iProc=0,nProc-1
        iRow0=rowDec(1,iProc)
        nRowx=rowDec(3,iProc)
c$OMP parallel default(shared)
c$OMP&    private(idata,iColx,iRowx)
c$OMP do 
        do idata=1,nRowx*nColx ! all data gotten
          iColx=(idata-1)/nRowx +1 !! for some Cta
          !iRowx=mod(idata-1,nRowx)+1 + (iRow0-1)
          iRowx=mod(idata-1,nRowx)+ iRow0 ! global index
          workTmp(iRowx, iColx)=rawBuf(idata,iProc)
        end do
c$OMP end do NOWAIT
c$OMP end parallel
      end do
      end subroutine unpackRawData
!====================================================================
      subroutine packRstData(cPart,iLoop,nLoop,colDecLoop)
      use decDatamd
      use decArraymd
      implicit none
      include 'mpif.h'
      ! dummy variables
      character*3 :: cPart
      integer :: iLoop, nLoop
      integer :: colDecLoop(3, 0:nProc-1, nLoop)

      ! local variables
      integer :: iProc, iRow0,nRowx,iRowx, nColx,iColx, ierr
      integer :: nsize, idata

      nColx=colDecLoop(3,taskNum,iLoop)
      do iProc=0,nProc-1
        iRow0=rowDec(1,iProc)
        nRowx=rowDec(3,iProc)
!$OMP parallel default(shared)
!$OMP&    private(idata,iColx,iRowx)
!$OMP do 
        do idata=1,nRowx*nColx ! all data gotten
          iColx=(idata-1)/nRowx +1 !! for some Cta
          !iRowx=mod(idata-1,nRowx)+1 + (iRow0-1)
          iRowx=mod(idata-1,nRowx)+ iRow0 ! global index
          rstBuf(idata,iProc)=workTmp(iRowx, iColx)
        end do
!$OMP end do NOWAIT
!$OMP end parallel
      end do

      end subroutine packRstData
!====================================================================
      subroutine putRstData(cPart,iLoop,nLoop,colDecLoop)
      use decDatamd
      use decArraymd
      implicit none
      include 'mpif.h'
      ! dummy variables
      character*3 :: cPart
      integer :: iLoop, nLoop
      integer :: colDecLoop(3, 0:nProc-1, nLoop)

      ! local variables
      integer(kind=MPI_ADDRESS_KIND) :: disp ! disp: remote offset in the window
      integer :: nsize !!
      integer :: iProc, nRowx,iCol0,nColx, ierr      

      iCol0=colDecLoop(1,taskNum,iLoop)
      nColx=colDecLoop(3,taskNum,iLoop)
!      print*,taskNum,' :nColx, ',iCol0,nColx,iLoop      
      do iProc=0,nProc-1
        nRowx=rowDec(3,iProc)
        disp =1_kAddr*nRowx*(iCol0-1)
        nsize=nRowx*nColx
        select case (cPart)
          case('asy')
            call MPI_Put(rstBuf(1,iProc),nsize,MPI_COMPLEX8, iProc,
     &                   disp,nsize,MPI_COMPLEX8, winAsy, ierr)
!            print*,taskNum,' dispA=',disp,nsize," to:",iProc
          case('int')
            call MPI_Put(rstBuf(1,iProc),nsize,MPI_COMPLEX8, iProc,
     &                   disp,nsize,MPI_COMPLEX8, winInt, ierr)
!            print*,taskNum,' dispI=',disp,nsize," to:",iProc
          case default
            print *,'putRstData: invalid value for cPart...'
            stop
        end select
      end do

      end subroutine putRstData
!====================================================================
      subroutine ExchangeData(cPart,iLoop,nLoop,colDecLoop,ifputrst)
      use decDatamd
      use decArraymd
      implicit none
      include 'mpif.h'
      ! dummy variables
      character*3 :: cPart
      integer :: iLoop, nLoop
      integer :: colDecLoop(3, 0:nProc-1, nLoop)
      integer :: ifputrst
      integer :: itest

!      if(taskNum==0) print*, "iLoop=",iLoop,"@",cPart
      if(iLoop == 1) then ! First Loop

        call getRawData(cPart,iLoop,nLoop,colDecLoop)
        call FinishExchange(cPart) ! new Raw data come
        call unpackRawData(cPart,iLoop,nLoop,colDecLoop)

        call getRawData(cPart,iLoop+1,nLoop,colDecLoop) ! the 2nd Loop
        return

      else if(iLoop /= nLoop+1) then ! inner Loops

        call FinishExchange(cPart) ! new Raw data come & old Rst data written.
        
        !// RstBuf is available, pack Rst data to buf. 
        !// this step frees the workTmp array.
        if(ifputrst /= 0) then
          call packRstData(cPart,iLoop-1,nLoop,colDecLoop) ! belong to the former Loop.
          call putRstData(cPart,iLoop-1,nLoop,colDecLoop)
        end if

        !// workTmp array is available, unpack new Raw data.
        !// this step frees the RawBuf.
        call unpackRawData(cPart,iLoop,nLoop,colDecLoop) ! for use in current Loop.

        !// RawBuf is available, pack Raw data for next Loop.
        if(iLoop+1 <= nLoop) then
          call getRawData(cPart,iLoop+1,nLoop,colDecLoop) ! for next Loop
        end if

        return
      else ! Last Loop, i.e. nLoop+1

        call FinishExchange(cPart) ! old Rst data written.

        !// RstBuf is available, pack Rst data to buf. 
        !// this step frees the workTmp array.
        if(ifputrst /= 0)  then
          call packRstData(cPart,iLoop-1,nLoop,colDecLoop) ! belong to the former Loop.
          call putRstData(cPart,iLoop-1,nLoop,colDecLoop)
        end if

        call FinishExchange(cPart) ! old Rst data written.
        return
      end if      

      end subroutine ExchangeData
!====================================================================
!====================================================================

