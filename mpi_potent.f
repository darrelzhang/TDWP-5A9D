        subroutine prepotent(NZ,NZINT,NDIM2asy,NDIM2int,NPOT)
        use decdatamd
        use decarraymd
        implicit none
        include 'mpif.h'

        integer   :: ierr
        integer, intent(in) :: NZ, NZINT, NDIM2asy, NDIM2int, NPOT

        allocate( asyColDistInfo((NZ-NZINT)*NDIM2asy) )
        allocate( asyColSizePerTask(0:nProc-1) )
        allocate( intColDistInfo(NZINT*NDIM2int) )
        allocate( intColSizePerTask(0:nProc-1) )
       asyColDistInfo(:)=0 ; asyColSizePerTask(:)=0
       intColDistInfo(:)=0 ; intColSizePerTask(:)=0

        nAsyMpiLoop=CEILING(DBLE((NZ-NZINT)*NDIM2asy)/DBLE(blockColsz))
        nIntMpiLoop=CEILING(DBLE(NZINT*NDIM2int)/DBLE(blockColsz))
        allocate( AsyColDecLoop(3, 0:nProc-1, nAsyMpiLoop+1), stat=ierr)
        allocate( IntColDecLoop(3, 0:nProc-1, nIntMpiLoop+1), stat=ierr)
        AsyColDecLoop(:,:,:)=0
        IntColDecLoop(:,:,:)=0
        if(ierr /= 0) STOP 'allocate coldecloop error in prepotent'

        call mpi_barrier(mpi_comm_world, ierr)
        call colDistri((NZ-NZINT)*NDIM2asy,asyColDistInfo,asyColSizePerTask,
     &                 AsyColDecLoop)
        call colDistri(NZINT*NDIM2int,     intColDistInfo,intColSizePerTask,
     &                 IntColDecLoop)
        call mpi_barrier(mpi_comm_world, ierr)

        asyCol=asyColSizePerTask(taskNum)
        intCol=intColSizePerTask(taskNum)
        print*, 'asycol=', asycol, 'intcol=', intcol

        allocate( asyPotSize(0:nProc-1,0:NPOT-1), 
     &            asyColG2L((NZ-NZINT)*NDIM2asy), 
     &            intPotSize(0:nProc-1,0:NPOT-1),
     &            intColG2L(NZINT*NDIM2int), stat=ierr )
        print*, 'ierr=', ierr
       asyPotSize(:,:)=0 ; asyColG2L(:)=0
       intPotSize(:,:)=0 ; intColG2L(:)=0 
        call colDistri2( (NZ-NZINT)*NDIM2asy, asyPotSize,asyColG2L,
     &                   asyColDistInfo,NPOT)
        call colDistri2( NZINT*NDIM2int,      intPotSize,intColG2L,
     &                   intColDistInfo,NPOT)

        return
        end subroutine prepotent
!================================================
        subroutine colDistri(numCol, colDistInfo, colSizePerTask,colDecLoop)
        use decDatamd
        implicit none

      ! Dummy variables
        integer :: numCol, colDistInfo(numCol),colSizePerTask(0:nProc-1),
     &             colDecLoop(3,0:nProc-1,1)

      ! Local variables
        integer :: L, i,j,k, numColLeft,L0,L1,L2

        colDistInfo(:) = -1
        colSizePerTask(:) = 0

        L=0
        do i=1,numCol/blockColsz
        do j=0,nProc-1
           colDecLoop(1,j,i)=L+1
           colDecLoop(2,j,i)=L+kCol
           colDecLoop(3,j,i)=kCol
           do k=1,kCol
              L=L+1
              colDistInfo(L)=j
              colSizePerTask(j)=colSizePerTask(j)+1
           end do
        end do
        end do

        L0 = L
        print*, L, 'numcol=', numcol
        if(L .LT. numCol) then
        numColLeft=numCol-L
        do j=0,nProc-1
           L1=numColLeft*j/nProc + 1 + L0
           L2=numColLeft*(j+1)/nProc + L0

           colDecLoop(1,j,i)=L1
           colDecLoop(2,j,i)=L2
           colDecLoop(3,j,i)=L2-L1+1
           do k=L1,L2
              L=L+1
              colDistInfo(L)=j
              colSizePerTask(j)=colSizePerTask(j)+1
           end do
        end do
        end if

        i = i + 1
        do j=0,nProc-1
           colDecLoop(1,j,i)=-1
           colDecLoop(2,j,i)=-1
           colDecLoop(3,j,i)=0
        end do

        if(taskNum == 0) then

        write(9001,*) '=========================='
        write(9001,*) 'numCol=',numCol , i

        do k=1,i
        do j=0,nProc-1
           write(9001,*) k,'@',j, colDecLoop(1,j,k),'<->',
     &                            colDecLoop(2,j,k),
     &                   'tot=',  colDecLoop(3,j,k)
        end do
        end do
        write(9001,*) '=========================='

        end if

        return
        end subroutine colDistri
!==========================================================
        subroutine colDistri2(numCol,PotSize,ColG2L,
     &                        ColDistInfo,NPOT)
        use decDatamd
        use decArraymd
        implicit none

        ! Dummy variables
        integer :: numCol,  NPOT
        integer :: PotSize(0:nProc-1,0:NPOT-1),ColG2L(numCol),ColDistInfo(numCol)

        ! Local variables
        integer :: counter(0:nProc-1)
        integer :: IPOT
        integer :: N1,N2, j,k, L,cLoc

        potSize(:,:)=0
        ColG2L(:)=-1
        counter(:)=0

       do IPOT=0,NPOT-1
           N1=numcol*IPOT/NPOT+1    
           N2=numcol*(IPOT+1)/NPOT  
           do L=N1, N2              

              cLoc=ColDistInfo(L)

              PotSize(cLoc,IPOT)=PotSize(cLoc,IPOT)+1

              counter(cLoc)=counter(cLoc)+1
              ColG2L(L)=counter(cLoc)
          end do
        end do

        if(taskNum == 0) then
           write(9000,*)'===== potsize ====='
           do IPOT=0,NPOT-1
           do j=0,nProc-1
              write(9000,*) 'IPOT=',IPOT, '@',j,'total:',PotSize(j,IPOT)
           end do
           end do

           write(9000,*) '===== PotSize per task ====='
           do j=0,nProc-1
              write(9000,*) '@',j,':',sum( PotSize(j,0:NPOT-1) )
           end do

           write(9000,*) '===== colG2L ====='
           do L=1,numCol
              write(9000,*) 'col index=',L,'@',ColDistInfo(L),
     &                      'local index=',ColG2L(L)
           end do
        end if

        return
        end subroutine colDistri2
!==========================================================
!==========================================================
        subroutine read_nh4pot(idout, nz1, nz2, nr, vref, 
     &                         zq, rq, npot, vm,
     &                         ColDistInfo, ColG2L, potSize)
	use comparamd, only : row1ty, row2ty, thetaty,
     &                        betasty, gammasty 
        use nh4rotmd,  only : na1node, nb1node
        use smallmd
        use decDatamd
        use decArraymd
        implicit none

        include 'mpif.h'

! Dummy variables
        integer :: idout
        integer, intent(in) :: nz1,nz2, nr, npot
        integer :: ColDistInfo(1), ColG2L(1),
     &             potSize(0:nProc-1,0:NPOT-1)
        real(8) :: vref(1:nr),zq(1:nz2), rq(1:nr)
        real(4) :: vm(1)

! Local variables
        real(4) :: vvmin(1:NR,1:NZ2-NZ1+1),
     &             tvvmin(1:NR,1:NZ2-NZ1+1), vv
        integer :: ipot, ios,ierr, ms, target0
        integer :: i, j, k, l, ia
        integer(8) :: ibase
        integer :: n1, n2, n5d 
        integer :: i1, i2, i3, i4, nleft
        integer :: potpos(0:npot-1)
        integer :: stat(mpi_status_size)
        real(4), allocatable :: vmtmp(:)

        n5d=betasty%nvb*gammasty%nvb*thetaty%nvb*na1node*nb1node
        vvmin(:,:)=1.d+5
        print*, 'idout=', idout
        print*, 'nz1=',nz1,'nz2=',nz2

        if(idout == 31 .or. idout == 32) then
            if(idout == 31) then
            open(300+tasknum,file='result/'//'VM'//itoc4(tasknum)//'.ASY',
     &           form='binary',status='old')
              l=asycol
            else if(idout == 32) then
            open(300+tasknum,file='result/'//'VM'//itoc4(tasknum)//'.INT',
     &           form='binary',status='old')
              l=intcol
            end if

            print*,'Node',taskNum,'Reading In ','VM'//ITOC4(taskNum),
     &                                          1_8*L*n5d

            do i=1,L
            read(300+taskNum) VM(1_8*(i-1)*n5d+1_8 : 1_8*i*n5d)
            end do

            vv=minval( VM(1_8:1_8*L*n5d) )
            print*,'@',taskNum,': vmin of VM=',vv
            close(300+taskNum)
            
            return
        endif   

!read POT...{ASY,INT}
        allocate(vmtmp(1:n5d), stat=ierr)
        potPos(:) = -1
        do ipot=0,npot-1

          if(tasknum == 0) print*, 'ipot=', ipot
          if(idout == 1) then
            open(200+ipot, file='result/'//'POT'//itoc4(ipot)//'.ASY', 
     &                     status='old',form='binary', iostat=ios)   !###
          else if(idout == 2) then
            open(200+ipot, file='result/'//'POT'//itoc4(ipot)//'.INT',
     &                     status='old',form='binary', iostat=ios)   !###
          else
             stop 'invalid value for idout'
          end if

        if(ios == 0) then
          potPos(ipot)=taskNum ! file exists
          do i=0,nProc-1
            if(i == taskNum) cycle
            call mpi_send(potpos(ipot), 1, mpi_integer,
     &                    i, ipot, mpi_comm_world, ierr)
          end do
        else
            call mpi_recv(potpos(ipot), 1, mpi_integer,
     &             mpi_any_source, ipot, mpi_comm_world,
     &             stat, ierr)
        end if
        end do !...do ipot=0,npot-1

        do ipot=0, npot-1
           if(potpos(ipot) == -1) then
              print*,'potPos',ipot,' : Error Value'
              stop
            end if
            if(tasknum == 0)
     &        print*,'POT'//itoc4(ipot), potpos(ipot)
        end do

        do ipot=0,npot-1
           if(potpos(ipot) == tasknum) then ! i have this POT file
           print*, 'ipot:', ipot

             n1=(nz2-nz1+1)*nr*row1ty%nvb*row2ty%nvb*ipot/npot+1
             n2=(nz2-nz1+1)*nr*row1ty%nvb*row2ty%nvb*(ipot+1)/npot      

             do ms=n1,n2
                i1=(ms-1)/(nr*row1ty%nvb*row2ty%nvb)+1
                nleft=mod(ms-1,nr*row1ty%nvb*row2ty%nvb)+1
                i2=(nleft-1)/(row1ty%nvb*row2ty%nvb)+1
                nleft=mod(nleft-1,row1ty%nvb*row2ty%nvb)+1
                i3=(nleft-1)/row2ty%nvb+1
                nleft=mod(nleft-1,row2ty%nvb)+1
                i4=nleft
                
               
                read(200+ipot)(vmtmp(ia),ia=1,n5d)
                vv=minval(vmtmp)
                vvmin(i2,i1)=min(vv, vvmin(i2,i1))

                vmtmp(:)=vmtmp(:)-vref(i2)-row1ty%vref(i3)-row2ty%vref(i4)

                target0=ColDistInfo(ms)

                print*,ms,'==>',target0,'@',ColG2L(ms),' size=',n5d

                if(target0 == taskNum) then
                  do ia=1,n5d
                     vm( 1_8*(colg2l(ms)-1)*n5d+ia)=vmtmp(ia)
                  end do
                else
                  call mpi_send(vmtmp, n5d, mpi_real4,
     $                      target0, ms, mpi_comm_world, ierr)
                end if
             end do 
           close(200+ipot)

           else  ! if(potPos(IPOT) /= taskNum)

            do i=1,potsize(tasknum, ipot)
               call mpi_recv(vmtmp, n5d, mpi_real4,
     $                potpos(ipot), mpi_any_tag, mpi_comm_world,
     $                stat, ierr)
               ms=stat(mpi_tag)
               do ia=1,n5d
                  vm( 1_8*(colg2l(ms)-1)*n5d + ia ) = vmtmp(ia)
               end do
            end do

           end if  ! if(potPos(IPOT) /= taskNum) 

            call mpi_barrier(mpi_comm_world, ierr)
        end do !...do IPOT=0,NPOT-1
        
        call mpi_reduce(vvmin,tvvmin, (nz2-nz1+1)*nr, mpi_real4,
     &                   mpi_min, 0, mpi_comm_world, ierr)
        if(tasknum == 0) then
        do i=1,nz2-nz1+1
           do j=1,nr
              write(30+idout,*) zq(i+nz1-1),rq(j), tvvmin(j,i)*27.211386d0
           end do
           write(30+idout,*)
              write(40+idout,*) zq(i+nz1-1), minval(tvvmin(:,i))*27.211386d0
        end do
        close(30+idout)
        end if

c      return

        if(idout == 1) then
        open(300+taskNum,file='result/'//'VM'//ITOC4(taskNum)//'.ASY',
     &       form='binary',status='UNKNOWN')
          L=asyCol
        else if(idout == 2) then
        open(300+taskNum,file='result/'//'VM'//ITOC4(taskNum)//'.INT',
     &       form='binary',status='UNKNOWN')
          L=intCol
        end if

        print*,'Node',taskNum,'Writing Out ','VM'//ITOC4(taskNum), L*n5d

        do i=1,L
          write(300+taskNum) VM(1_8*(i-1)*n5d+1 :1_8*i*n5d)
        end do

        close(300+taskNum)

        vv=minval( VM(1_8:1_8*L*n5d) )
        print*,'@',taskNum,': vmin of VM=',vv

        deallocate(vmtmp, stat=ierr)
        return
        end subroutine read_nh4pot
!==========================================================



!==========================================================
