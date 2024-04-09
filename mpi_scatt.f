!====================================================================
	subroutine scatt_start_mpi
	use decdatamd
	use decarraymd
	implicit none

	include 'mpif.h'

	integer  :: ierr
	integer  :: provided

	 call mpi_init(ierr)
!	 call mpi_init_thread(mpi_thread_multiple, provided, ierr)

	 call mpi_comm_rank(mpi_comm_world, tasknum, ierr)
	 call mpi_comm_size(mpi_comm_world, nproc, ierr)

	 if(tasknum == 0) then
!	   print*,'provided=',provided
!	   print*,'mpi_thread_multiple=',mpi_thread_multiple 
	 endif   

	 if(tasknum == 0) write(9000,*) 'nproc=',nproc

	  print*, 'node:', tasknum 

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
	fullbufsz  = 2*kcol*int(ceiling(dble(nrow)/dble(nproc)))
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

	allocate(recvbuf(fullbufsz, 1:nproc-1), sendbuf(fullbufsz, 1:nproc-1),
     &           sendrecvhandle(1:2*(nproc-1)), sendbufsize(0:nproc-1),
     &           handlestat(MPI_STATUS_SIZE,1:2*(nproc-1)),
     &           workTmp(nrow,kcol), stat=ierr)
	if(ierr /= 0) STOP 'Failed: allocate in rowdecomp'

	workTmp(:,:)=(0.0,0.0)
	recvbuf(:,:)=(0.0,0.0)
	sendbuf(:,:)=(0.0,0.0) 

	return
	end subroutine rowdecomp
!====================================================================
	subroutine SendRecv(iMpiLoop,nMpiLoop, ColDecLoop,wfMatrix)
	use decDatamd
	use decArraymd
	implicit none
	include 'mpif.h'

! Dummy Variables
	integer :: iMpiLoop, nMpiLoop
	integer :: ColDecLoop(3,0:nProc-1,1)
	complex(4) :: wfMatrix(1)

! Local variables
	integer :: j, k,L, ierr
        integer :: M
	integer :: iCol, nnCol

	if(iMpiLoop == 1) then  ! First Loop

	sendBufSize(:) = 0
	call packRaw(wfMatrix(1), colDecLoop(1:3,0:nProc-1,1))
	call initSendRecv(MPI_COMPLEX)
	call MPI_WaitAll(2*(nProc-1), sendRecvHandle, handleStat, ierr)


	call unpackRaw(iMpiLoop, colDecLoop(3,taskNum,1))

	!copy local wf to workTmp
	M = (NAp) * (colDecLoop(1,taskNum,1)-1) 
	do k=1,colDecLoop(3,taskNum,1)                   
	  do j=rowDec(1,taskNum),rowDec(2,taskNum)        
	    M = M + 1                                
	    workTmp(j,k) = wfMatrix(M)                  
	  end do                                         
	end do                                           


	sendBufSize(:) = 0
	call packRaw(wfMatrix(1), colDecLoop(1:3,0:nProc-1,2))
	call initSendRecv(MPI_COMPLEX)
	return

	else if(iMpiLoop /= nMpiLoop+1) then

	   call MPI_WaitAll(2*(nProc-1), sendRecvHandle, handleStat, ierr)
	   sendBufSize(:) = 0

        if(iMpiLoop >= 3) then
	  call unpackComp(wfMatrix(1), colDecLoop(1:3,0:nProc-1,iMpiLoop-2))
        end if
     

	call packComp(colDecLoop(3,taskNum,iMpiLoop-1))

        ! save comp data of LAST loop to local WF
        M = (NAp) * (colDecLoop(1,taskNum,iMpiLoop-1) - 1) 
        do j=1, colDecLoop(3,taskNum,iMpiLoop-1)                        
	  do k=rowDec(1,taskNum), rowDec(2,taskNum)                     
	    M = M + 1                                               
	    wfMatrix(M) = workTmp(k,j)                             
          end do                                                    
        end do                                                      


        call unpackRaw(iMpiLoop, colDecLoop(3,taskNum,iMpiLoop))

        !copy local wf to workTmp
        M = (NAp) * (colDecLoop(1,taskNum,iMpiLoop) - 1)   
        do j=1,colDecLoop(3,taskNum,iMpiLoop)                       
          do k=rowDec(1,taskNum),rowDec(2,taskNum)                  
            M = M + 1                                           
            workTmp(k,j) = wfMatrix(M)                             
          end do                                                    
        end do                                                      


        if(iMpiLoop <= nMpiLoop-1) then
          call packRaw(wfMatrix(1), colDecLoop(1:3,0:nProc-1,iMpiLoop+1))
	end if

	call initSendRecv(MPI_COMPLEX)

      else  ! Last Loop

        call MPI_WaitAll(2*(nProc-1), sendRecvHandle, handleStat, ierr)
        sendBufSize(:)=0

        call unpackComp(wfMatrix(1), colDecLoop(1:3,0:nProc-1,iMpiLoop-2))


        call packComp( colDecLoop(3,taskNum,iMpiLoop-1) )

        ! save comp data of LAST loop to local WF
        M = (NAp)*(colDecLoop(1,taskNum,iMpiLoop-1) - 1)   
        do j=1, colDecLoop(3,taskNum,iMpiLoop-1)                    
          do k=rowDec(1,taskNum), rowDec(2,taskNum)                 
            M = M + 1                                           
            wfMatrix(M) = workTmp(k,j)                             
          end do                                                    
        end do                                                      


        call initSendRecv(MPI_COMPLEX)
        call MPI_WaitAll(2*(nProc-1), sendRecvHandle, handleStat, ierr)
        call unpackComp(wfMatrix(1), colDecLoop(1:3,0:nProc-1,iMpiLoop-1))

      end if

      end subroutine SendRecv
!====================================================================
	subroutine packRaw(localWF, colDec)
	use decDatamd
	use decArraymd
	implicit none

! Dummy variables
      integer :: colDec(3,0:nProc-1)
      complex :: localWF(1)

      ! Local variables
      integer :: j,k, size, n,p

      k=0
      do j=0,nProc-1
        if(j == taskNum) cycle

        k = k + 1
        size = colDec(3,j) * NAp

        p = (NAp) * (colDec(1,j)-1)            
        do n=sendBufSize(j)+1, sendBufSize(j)+size      
          p = p + 1                                 
          sendBuf(n,k) = localWF(p)                    
        end do
        sendBufSize(j) = sendBufSize(j) + size
      end do

      end subroutine packRaw
!====================================================================
	SUBROUTINE initSendRecv(dataType)
	USE decDatamd
	USE decArraymd
	IMPLICIT NONE
	include 'mpif.h'

!  dummy variables
	INTEGER, INTENT(IN) :: dataType

!  subroutine variables
	INTEGER :: j, k,L, ierr

! initiate sending/receiving of data
	k = 0; L = 0
	do j = 0, nProc-1
	   if (j == taskNum) CYCLE

	   k = k + 1
	   call MPI_Isend(sendBuf(1, k), sendBufSize(j), dataType,
     &                    j, 1, MPI_COMM_WORLD,
     &                    sendRecvHandle(k+L), ierr)

	   L = L + 1
	   call MPI_Irecv(recvBuf(1, L), fullBufsz, dataType,
     &                    j, 1, MPI_COMM_WORLD,
     &                    sendRecvHandle(k+L), ierr)
	end do

	return
	END SUBROUTINE initSendRecv
!====================================================================
      subroutine unpackRaw(iMpiLoop, nCurrCol)
      use decDatamd
      use decArraymd
      implicit none

      ! Dummy variables
      integer :: iMpiLoop, nCurrCol

      ! Local variables
      integer :: j, L,M, n,p

      L = 0
      do j=0, nProc-1
        if(j == taskNum) cycle

        M=rowDec(4,taskNum)
        if(iMpiLoop <= 2) M = 0  ! when there's NO comp data 

        L = L + 1
        do n=1,nCurrCol
          do p=rowDec(1,j), rowDec(2,j)
            M = M + 1
            workTmp(p,n) = recvBuf(M,L)
          end do
        end do
      end do !... do j=0,nProc-1

      end subroutine unpackRaw
!====================================================================
      subroutine unpackComp(localWF, colDec)
      use decDatamd
      use decArraymd
      implicit none

      ! Dummy variables
      integer    :: colDec(3,0:nProc-1)
      complex(4) :: localWF(1)

      ! Local variables
      integer :: j, L, n, k

      L=0
      do j=0,nProc-1
        if(j == taskNum) cycle

        L = L + 1
        n = (NAp) * (colDec(1,j)-1)                    

        localWF( n+1 : n+(NAp)*(colDec(3,j)) ) =    
     &          recvBuf(1 : (NAp)*(colDec(3,j)), L)   
      end do

      end subroutine unpackComp
!====================================================================
      subroutine packComp(nCurrCol)
      use decDatamd
      use decArraymd
      implicit none

      ! Dummy variables
      integer :: nCurrCol

      ! Local variables
      integer :: j, L,M, n,p
      integer :: size0

      L=0
      do j=0, nProc-1
        if(j == taskNum) cycle

        L = L + 1
        size0 = nCurrCol*rowDec(3,j)

        M = sendBufSize(j)
        do n=1,nCurrCol
          do p=rowDec(1,j), rowDec(2,j)
            M = M + 1
            sendBuf(M,L) = workTmp(p,n)
          end do
        end do

        sendBufSize(j) = sendBufSize(j) + size0
      end do

      end subroutine packComp
!====================================================================

