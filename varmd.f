        module decdatamd
        implicit none

        integer,parameter :: kAddr=8
        integer              :: tasknum, nproc
        integer              :: ncol, nrow, kcol, asycol, intcol
        integer              :: fullbufsz, blockcolsz
        integer              :: NAp
        integer, allocatable :: rowdec(:,:)

        integer              :: nasympiloop, nintmpiloop
        integer,allocatable  :: asycoldecloop(:,:,:),
     &                          intcoldecloop(:,:,:)

        integer, allocatable :: asyColDistInfo(:), intColDistInfo(:)
        integer, allocatable :: asyColSizePerTask(:),
     &                          intColSizePerTask(:)
        integer, allocatable :: asyPotSize(:,:), intPotSize(:,:)
        integer, allocatable :: asyColG2L(:), intColG2L(:)

        real(8), allocatable :: tCour2D(:,:), Cour2D(:,:)
        integer :: winAsy, winInt

        end module decdatamd

        module decarraymd
        implicit none

!        integer, allocatable :: handlestat(:,:)
!        integer, allocatable :: sendrecvhandle(:)
!        integer, allocatable :: sendbufsize(:)

        complex, allocatable :: workTmp(:,:)
!        complex, allocatable :: sendbuf(:,:), recvbuf(:,:)

        complex(4), allocatable :: rawBuf(:,:), rstBuf(:,:)

        end module decarraymd


