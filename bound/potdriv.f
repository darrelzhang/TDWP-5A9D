	module potdrivmd
	implicit none

	public pot_zero, potdriv

	contains

	subroutine pot_zero
	use massmd,   only : amass
	use comparamd,only : cart0,
     &                       zg0, betanh30,alphas0, 
     &                       r10, betas0,gammas0,  
     &                       row10, row20, theta0,
     &                       vinf
	implicit none
        integer :: i, j
	real(8) :: ev, vx
	real(8) :: zz
        real(8) :: cart(3,5)
	real(8), external :: fzhangxr


	 ev=27.211386d0
	 call cart_hnh3_to_radau(cart0,zg0,betanh30,alphas0,
     &                           r10,betas0,gammas0,
     &                           row10,row20,theta0)

!nh4 pes lijun
!When the potential energy surface needs to initialize some variables, 
!we will need to call this program. Otherwise, there is no need to call it.
         call pes_init

	 call potdriv(zg0,betanh30,alphas0,
     &                     r10,betas0,gammas0, 
     &                     row10,row20,theta0,   vx)
	 vinf=vx

	 write(*,*)'===================================================' 
	 write(*,*)'initial radau coordinate (in bohr):'
	 write(*,*)'---------------------------------------------------' 
	 write(*,"(1x,'zg0=',f10.6)")zg0
	 write(*,"(1x,'betanh30=',f10.6)")betanh30
	 write(*,"(1x,'alphas0=',f10.6)")alphas0
	 write(*,*)'---------------------------------------------------' 
	 write(*,"(1x,'r10=',f10.6)")r10
	 write(*,"(1x,'betas0=',f10.6)")betas0
	 write(*,"(1x,'gammas0=',f10.6)")gammas0
	 write(*,*)'---------------------------------------------------' 
	 write(*,"(1x,'row10=',f10.6)")row10
	 write(*,"(1x,'row20=',f10.6)")row20
	 write(*,"(1x,'theta0=',f10.6)")theta0   
	 write(*,*)'---------------------------------------------------' 
	 write(*,*)'initial potential energy of H+NH3 (in eV):',vinf*ev 
	 write(*,*)'===================================================' 

	return
	end subroutine pot_zero


	subroutine nh3_pot_eng_suf
	use comparamd, only : row1ty, row2ty, thetaty,
     &                        betasty, gammasty, 
     &                        r1ty,  zg0, betanh30,alphas0
	use nh3_bound_eigmd, only : nh3pes
        use nh3_bound_eigmd, only : nbbound
	implicit none
        integer :: irow2, irow1, ith, ig, ib, ir1, ibas
        integer :: iz, ibnh3, ianh3
        integer :: ms, ierr
        real(8) :: row2, row1, theta, gammas, betas, r1
        real(8) :: zz, betanh3, alphanh3
        real(8) :: vx, vrmin(row2ty%nvb,row1ty%nvb), vtmp, vmin
        
        allocate(nh3pes(row2ty%nvb * row1ty%nvb * thetaty%nvb *
     &               gammasty%nvb * betasty%nvb *r1ty%nvb ) ,
     &          stat=ierr)

        ms=0
        do ir1=1, r1ty%nvb
            r1=r1ty%rqa(ir1)
        do ib =1, betasty%nvb
            betas=betasty%rqa(ib)
        do ig =1, gammasty%nvb
            gammas=gammasty%rqa(ig)
        do ith=1, thetaty%nvb 
            theta=thetaty%rqa(ith)
        do irow1=1, row1ty%nvb
            row1=row1ty%rqa(irow1)
        do irow2=1, row2ty%nvb
            row2=row2ty%rqa(irow2)
         
            ms=ms+1

	    call potdriv(zg0,betanh30,alphas0,
     &                     r1,betas,gammas, 
     &                     row1,row2,theta,   vx)
            nh3pes(ms)=vx-r1ty%vref(ir1)-row1ty%vref(irow1)-row2ty%vref(irow2)

	enddo 
	enddo 
	enddo 
	enddo 

	enddo 
	enddo

        vmin=minval(nh3pes)
        print*, 'minimum value of nh3pes:', vmin*27.211386d0

	return
	end subroutine nh3_pot_eng_suf


	subroutine potdriv(zz,betanh3,alphas,
     &                     r1,betas,gammas, 
     &                     row1,row2,theta, vx)
	use comparamd, only : vinf, zg0
	implicit none
        integer, parameter  :: num_cart=15
        real(8)             :: zz1  
	real(8), intent(in) :: zz, betanh3,alphas,
     &                         r1,betas,gammas, 
     &                         row1,row2,theta
	real(8), intent(out):: vx
	real(8) :: cart(3,5), cart1(3,5), der_cart(num_cart)
        real(8) :: r_in(6) 
        real(8) :: vxa, vxb, vxc

        zz1=zz
!N H H H H
	cart=0.d0
	 call radau_hnh3_to_cart(zz1,betanh3,alphas,
     &                           r1,betas,gammas, 
     &                           row1,row2,theta, cart)

!nh4 pes lijun
!atoms given in cart is N H1 H2   H3     H4, and cart is in bohr
!the output energy is vx, which is in Hatree.
!one can change the sequence of the atoms to call the PES
        call rearrange_atom(cart,cart1)
        cart1=cart1*0.529177d0
        call nh4pipNN(cart1,vx, vxa, vxb, vxc)
        vx=vx/27.211386d0
        vx=vx-vinf

	return
	end subroutine potdriv

!==========================================================
!N H1 H2   H3     H4
        subroutine  radau_hnh3_to_cart(rs,beta_nh3,alpha_s,
     &                                 r1,beta_s, gamma_s,
     &                                 row1,row2,theta,cart0)
        use massmd, only : amass
        implicit none
        integer  :: i, j
        real(8), intent(in)   :: rs,   beta_nh3, alpha_s, 
     &                          r1,   beta_s, gamma_s, 
     &                          row1, row2,   theta
        real(8), intent(out)  :: cart0(3,5)
        real(8)   :: cent_nh2(3), cent_nh3, cent_hnh3

!N H1 H2
        call radau_cart(amass(1:3), row1, row2, theta, cart0(1,1))

!COM of NH2        
         do i=1,3
            cent_nh2(i)=0.d0
            do j=1,3
               cent_nh2(i)=cent_nh2(i)+cart0(i,j)*amass(j)
            enddo
               cent_nh2(i)=cent_nh2(i)/(sum(amass(1:3)))
         enddo
 
!move COM of NH2 to origin
         do j=1, 3
            do i=1,3
               cart0(i,j)=cart0(i,j)-cent_nh2(i)
            enddo
         enddo 

        call roteular(3,cart0,0.d0,beta_s,gamma_s) 
        
!H3
        cart0(1,4)=0.d0
        cart0(2,4)=0.d0
        cart0(3,4)=r1

!move center of mass of NH3 to origin
        cent_nh3=cart0(3,4)*amass(4)/(sum(amass(1:4)))
        do j=1, 4
           cart0(3,j)=cart0(3,j)-cent_nh3
        enddo 
         
        call roteular(4,cart0,0.d0,beta_nh3,alpha_s) 

!H4
        cart0(1,5)=0.d0
        cart0(2,5)=0.d0 
        cart0(3,5)=rs 
        
!!move COM of HNH3 to origin 
!        cent_hnh3=cart0(3,5)*amass(5)/(sum(amass(1:5))) 
!        do j=1,5
!           cart0(3,j)=cart0(3,j)-cent_hnh3
!        enddo
      
        return
        end subroutine radau_hnh3_to_cart
!==========================================================
        subroutine cart_hnh3_to_radau(cart0,rs,beta_nh3,alpha_s,
     &                        r1,beta_s,gamma_s,
     &                        row1,row2,theta)
        use massmd, only : amass,atom
        implicit none
        integer :: i, j
        real(8), intent(in)  :: cart0(3,5)
        real(8)              :: cart1(3,5), cart2(3,5), cart3(3,5),
     &                          cart4(3,5)
        real(8), intent(out) :: rs,   beta_nh3, alpha_s,   
     &                          r1,   beta_s, gamma_s,    
     &                          row1, row2,   theta 
        real(8)              :: alpha, beta, phi
        real(8)              :: rh4, rh3, rh2, rh1
        real(8)              :: GHNH3(3), GNH3(3), GNH2(3), Ap(3)
        real(8)              :: r1_vec(3), r2_vec(3)
        
!COM of HNH3
        do i=1, 3
           GHNH3(i)=0.d0
           do j=1,5
              GHNH3(i)=GHNH3(i)+amass(j)*cart0(i,j) 
           enddo
              GHNH3(i)=GHNH3(i)/(sum(amass(1:5)))
        enddo

!move COM of HNH3 to origin
        do j=1,5
           do i=1,3
              cart1(i,j)=cart0(i,j)-GHNH3(i)
           enddo
        enddo 
        rh4=sqrt(sum(cart1(:,5)**2))
        beta=acos(cart1(3,5)/rh4)
        alpha=acos(cart1(1,5)/rh4/sin(beta))
        if(cart1(2,5) < 0.d0) alpha=2.d0*acos(-1.d0)-alpha

        call roteular(5,cart1,0.d0,-beta,-alpha) 
        
!@        write(41,*)5
!@        write(41,*)
!@        do j=1, 5
!@           write(41,"(a1,3f10.5)")atom(j),(cart1(i,j),i=1,3)
!@        enddo


!COM of NH3
        do i=1, 3
           GNH3(i)=0.d0
           do j=1,4
              GNH3(i)=GNH3(i)+amass(j)*cart1(i,j) 
           enddo
              GNH3(i)=GNH3(i)/(sum(amass(1:4)))
        enddo

!move NH3 COM to origin
        do j=1,5
           do i=1,3        
             cart2(i,j)=cart1(i,j)-GNH3(i)
           enddo
        enddo 
        rs=sqrt(sum(cart2(:,5)**2))

        rh3=sqrt(sum(cart2(:,4)**2))
        beta_nh3=acos(cart2(3,4)/rh3)
        phi=acos(cart2(1,4)/rh3/sin(beta_nh3))
        if(cart2(2,4) < 0.d0) phi=2.d0*acos(-1.d0)-phi

        call roteular(5,cart2,0.d0,-beta_nh3,-phi) 
        
!@        write(51,*)5
!@        write(51,*)
!@        do j=1, 5
!@           write(51,"(a1,3f10.5)")atom(j),(cart2(i,j),i=1,3)
!@        enddo

!COM of NH2
        do i=1, 3
           GNH2(i)=0.d0
           do j=1, 3
              GNH2(i)=GNH2(i)+amass(j)*cart2(i,j)
           enddo
              GNH2(i)=GNH2(i)/sum(amass(1:3))
        enddo            

!move NH2 COM to origin, H3 is in Z-axis
        do j=1,5
           do i=1,3
              cart3(i,j)=cart2(i,j)-GNH2(i)
           enddo
        enddo 
        r1=sqrt(sum(cart3(:,4)**2))

        call cart_radau(amass, cart3(1,1), row1, row2, theta, Ap)

        do i=1,3
          r2_vec(i)=cart3(i,3)-Ap(i)
        enddo

        beta_s=(sum(r2_vec(:)*cart3(:,4)))/( row2*abs(cart3(3,4)) )
        beta_s=acos(beta_s)
        alpha_s=acos(r2_vec(1)/(sqrt(sum(r2_vec(:)**2)))/sin(beta_s))
        if(r2_vec(2) < 0.d0) alpha_s=2.d0*acos(-1.d0)-alpha_s
        call roteular(5,cart3,0.d0,-beta_s, -alpha_s) 
        call cart_radau(amass, cart3(1,1), row1, row2, theta, Ap)
!@        print*, "H2:", (cart3(i,3),i=1,3)
!@        print*, "A:",(Ap(i),i=1,3)

        do j=1, 5
           do i=1, 3
              cart4(i,j)=cart3(i,j)-Ap(i)
           enddo
        enddo
        
        do i=1, 3
           r1_vec(i)=cart4(i,2)
        enddo
!@        print*, 'row1=', row1, sqrt(sum(r1_vec(:)**2)) 

        gamma_s=acos(r1_vec(1)/(sqrt(sum(r1_vec(:)**2)))/sin(theta))
        if(r1_vec(2) < 0.d0)gamma_s=2.d0*acos(-1.d0)-gamma_s
!@        call roteular(5,cart4,0.d0, 0.d0, -gamma_s) 
!@      
!@
!@        write(61,"('5')")
!@        write(61,*)
!@        do j=1, 5
!@           write(61,"(a1,3f10.5)")atom(j),(cart4(i,j),i=1,3)
!@        enddo
!@        call cart_radau(amass, cart4(1,1), row1, row2, theta, Ap)
!@        write(61,"(a1,3f10.5)")"O",(Ap(i),i=1,3)
        

        return
        end subroutine cart_hnh3_to_radau 
!==========================================================
!A point : A-GH2=sqrt(N-GH2 * GNH2-GH2)
        subroutine cart_radau(amass, cart0, row1, row2, theta, Ap)
        implicit none
        integer                 :: i, j
        real(8), intent(in)      :: cart0(3,3)
        real(8), intent(in)      :: amass(3)
        real(8), intent(out)     :: row1, row2, theta, Ap(3)
        real(8)                  :: GNH2(3), GH2(3)
        real(8)                  :: facx, facy, facz
        real(8)                  :: r1, r2
        real(8)                  :: r1_vec(3),r2_vec(3)
        real(8)                  :: tmp

!mass center of NH2
        do i=1, 3
           tmp=0.d0
           do j=1,3
              tmp=tmp+cart0(i,j)*amass(j)
           enddo 
           GNH2(i)=tmp/sum(amass(:))
        enddo        

!mass center of H1-H2
        do i=1, 3
           tmp=0.d0
           do j=2, 3
              tmp=tmp+cart0(i,j)*amass(j)
           enddo
          GH2(i)=tmp/sum(amass(2:3)) 
        enddo

!N, A, GNH2, GH2 in the same line
        facx=cart0(1,1)-GH2(1)
        facy=cart0(2,1)-GH2(2)
        facz=cart0(3,1)-GH2(3)

!r1 == N-GH2 , r2 == GNH2-GH2
        r1=sqrt(sum((cart0(:,1)-GH2(:))**2))
        r2=sqrt(sum((GNH2(:)-GH2(:))**2))

!A(a,(a-GH2(1))*facy/facx-GH2(2),(a-GH2(1))*facz/facx-GH2(3))
!tmp == abs(a-GH2(1))
        tmp=sqrt(r1*r2/(1.d0+(facy/facx)**2+(facz/facx)**2))

!check A point in (N,GH2) or (GH2, N)
        if(cart0(1,1)>=GH2(1))then
!          if((tmp+GH2(1))>= cart0(1,1))then
!             Ap(1)=-tmp+GH2(1)
!          else
             Ap(1)=tmp+GH2(1)
!          endif
        else
!          if((tmp+GH2(1))>= GH2(1))then
             Ap(1)=-tmp+GH2(1)
!          else
!             Ap(1)=tmp+GH2(1)
!          endif
        endif
        Ap(2)=(Ap(1)-GH2(1))*facy/facx+GH2(2)
        Ap(3)=(Ap(1)-GH2(1))*facz/facx+GH2(3)

!R1, R2
        row1=sqrt(sum((Ap(:)-cart0(:,2))**2))
        row2=sqrt(sum((Ap(:)-cart0(:,3))**2))
        
        do i=1, 3
           r1_vec(i)=cart0(i,2)-Ap(i)
           r2_vec(i)=cart0(i,3)-Ap(i)
        enddo

!theta        
        theta=sum(r1_vec(:)*r2_vec(:))/(row1*row2)
        theta=acos(theta)


        return
        end subroutine cart_radau
c====================================================================
       subroutine radau_cart(amass, row1, row2, theta, cart0) 
       implicit none
       integer  :: i, j
       real(8), intent(in)  :: row1, row2, theta
       real(8), intent(in)  :: amass(3)
       real(8), intent(out) :: cart0(3,3) 
       real(8)              :: GNH2(3), GH2(3),Ap(3)
       real(8)              :: facx, facy, facz 
       real(8)              :: row1_vec(3), row2_vec(3) 
       real(8)              :: r1, r2 ,r3 
       real(8)              :: alpha
       real(8)              :: tmp

       
        Ap=0.d0

        cart0(1,2)=row1*sin(theta)
        cart0(2,2)=0.d0
        cart0(3,2)=row1*cos(theta)

        cart0(1,3)=0.d0
        cart0(2,3)=0.d0
        cart0(3,3)=row2
                

!mass center of H-H
        do i=1, 3
           tmp=0.d0
           do j=2, 3
              tmp=tmp+cart0(i,j)*amass(j)
           enddo
          GH2(i)=tmp/sum(amass(2:3)) 
        enddo

!N, A, GNH2, GH2 in the same line
        facx=Ap(1)-GH2(1)
        facy=Ap(2)-GH2(2)
        facz=Ap(3)-GH2(3)

!R1 -- A-H1, R2 -- A-H2
        do i=1, 3
           row1_vec(i)=cart0(i,2)-Ap(i)
           row2_vec(i)=cart0(i,3)-Ap(i)
        enddo

        alpha=(1.d0-sqrt(sum(amass(:))/amass(1)))/(amass(2)+amass(3))

!GNH2(a,(a-GH2(1))*facy/facx-GH2(2),(a-GH2(1))*facz/facx-GH2(3))
!row1_vec, x,==> a, (tmp is a)
        tmp=( row1_vec(1)-(1.d0-alpha*amass(2))*cart0(1,2)
     &      +alpha*amass(3)*cart0(1,3) )/
     &      (-1.d0+alpha*amass(2)+alpha*amass(3))
!        print*, 'r1, x', tmp
 
!GNH2(x, y, z)       
        GNH2(1)=tmp
        GNH2(2)=(GNH2(1)-GH2(1))*facy/facx+GH2(2)
        GNH2(3)=(GNH2(1)-GH2(1))*facz/facx+GH2(3)

!N atom, (x, y, z)
        do i=1,3
           cart0(i,1)=(GNH2(i)*(sum(amass(:)))
     &               -cart0(i,2)*amass(2)-cart0(i,3)*amass(3))/amass(1) 
        enddo

!r1 -- H1-GNH2 , r2 -- H2-GNH2       
!        do i=1, 3
!           grow1_vec(i)=cart0(i,2)-GNH2(i)
!           grow2_vec(i)=cart0(i,3)-GNH2(i)
!        enddo

!'row1_vec, z,' , tmp==a  ---for checking 
!        tmp=( (row1_vec(3)+alpha*amass(3)*(cart0(3,3)-GH2(3))-
!     &         (1.d0-alpha*amass(2))*(cart0(3,2)-GH2(3)))*facx/facz
!     &        +alpha*amass(3)*GH2(1)-(1.d0-alpha*amass(2))*GH2(1) )/
!     &      (alpha*amass(2)+alpha*amass(3)-1.d0)
!        print*, 'r1, z', tmp

!--------------------------------------

!'row2_vec, x,', tmp==a  ---for checking
!        tmp=( row2_vec(1)+alpha*amass(2)*cart0(1,2)-
!     &       (1.d0-alpha*amass(3))*cart0(1,3) )/
!     &      (alpha*(amass(2)+amass(3))-1.d0) 
!        print*, 'r2, x', tmp

!'row2_vec,z',  tmp==a  ---for checking
!        tmp=( (row2_vec(3)+alpha*amass(2)*(cart0(3,2)-GH2(3))-
!     &        (1.d0-alpha*amass(3))*(cart0(3,3)-GH2(3))  )*facx/facz 
!     &       +(alpha*(amass(2)+amass(3))-1.d0)*GH2(1) )/  
!     &       (alpha*(amass(2)+amass(3))-1.d0) 
!        print*, 'r2, z', tmp

!N, H1, H2
!        do j=1,3
!           write(22,*) (cart0(i,j),i=1,3)
!        enddo

        
!      r1=sqrt(sum((cart0(:,1)-GH2(:))**2))
!      r2=sqrt(sum((GNH2(:)-GH2(:))**2))  
!      r3=sqrt(sum((Ap(:)-GH2(:))**2))  
!      print*,'check radau point:', r1*r2, r3*r3  


       return 
       end subroutine radau_cart        
!********************************************************************
!*   Program    :  roteular                                         *
!*   Function   :  rotate a vector by eular angle (a,b,c)           *
!********************************************************************
	subroutine roteular(natom,cartes,a,b,c)
	implicit none
	integer :: iatom,ixyz,jxyz
	real(8) :: a,b,c,tmatrot(3,3)
	real(8) :: sina,cosa,sinb,cosb,sinc,cosc,temp(3)
	integer, intent(in)   :: natom
	real(8), intent(inout):: cartes(3,natom)

        sina=dsin(a); cosa=dcos(a)
        sinb=dsin(b); cosb=dcos(b)
        sinc=dsin(c); cosc=dcos(c)
 
        tmatrot(1,1)= cosa*cosb*cosc-sina*sinc
        tmatrot(2,1)= sina*cosb*cosc+cosa*sinc
        tmatrot(3,1)=-sinb*cosc
 
        tmatrot(1,2)=-cosa*cosb*sinc-sina*cosc
        tmatrot(2,2)=-sina*cosb*sinc+cosa*cosc
        tmatrot(3,2)= sinb*sinc
 
        tmatrot(1,3)= cosa*sinb
        tmatrot(2,3)= sina*sinb
        tmatrot(3,3)= cosb

	do iatom=1,natom
	   do ixyz=1,3
	      temp(ixyz)=0.0d0
	      do jxyz=1,3
	         temp(ixyz)=temp(ixyz)
     &                    +tmatrot(ixyz,jxyz)*cartes(jxyz,iatom)
	      enddo
	   enddo
	   do ixyz=1,3
	      cartes(ixyz,iatom)=temp(ixyz)
	   enddo
 	enddo

	return
	end subroutine roteular


        subroutine rearrange_atom(cartin, cartout)
        implicit none
        integer, parameter :: natom=5
        integer :: i, new(natom)
        real(8), intent(in)  :: cartin(3,natom)
        real(8), intent(out) :: cartout(3,natom)
        
        new(1:natom)=[5,1,2,3,4]
        do i=1, 5
           call dcopy1d(3,cartin(1,i),cartout(1,new(i)))
        enddo

        return
        end subroutine rearrange_atom
        subroutine dcopy1d(ndim,xin,xout)
        implicit none
        integer, intent(in) :: ndim
        real(8), intent(in) :: xin(ndim)
        real(8), intent(out):: xout(ndim)

        xout(1:ndim)=xin(1:ndim)

        return
        end subroutine dcopy1d

        subroutine cart_bond_chenjun(cart,rout)
        implicit none
        real(8), intent(in)  :: cart(3,4)
        real(8), intent(out) :: rout(6)
        real(8) :: vec1(3), vec2(3)
        real(8) :: radtoang

        rout(1)=sqrt(sum((cart(:,1)-cart(:,2))**2))
        rout(2)=sqrt(sum((cart(:,1)-cart(:,3))**2))
        rout(3)=sqrt(sum((cart(:,1)-cart(:,4))**2))

        rout(4)=sqrt(sum((cart(:,2)-cart(:,3))**2))
        rout(5)=sqrt(sum((cart(:,2)-cart(:,4))**2))
        rout(6)=sqrt(sum((cart(:,3)-cart(:,4))**2))

        return
        end subroutine cart_bond_chenjun

        subroutine cart_bond(cart,rout)
        implicit none
        real(8), intent(in)  :: cart(3,4)
        real(8), intent(out) :: rout(6)
        real(8) :: vec1(3), vec2(3)
        real(8) :: radtoang

        rout(1)=sqrt(sum((cart(:,1)-cart(:,2))**2))*0.529177d0
        rout(2)=sqrt(sum((cart(:,1)-cart(:,3))**2))*0.529177d0
        rout(3)=sqrt(sum((cart(:,1)-cart(:,4))**2))*0.529177d0

!        rout(4)=sqrt(sum((cart(:,2)-cart(:,3))**2))
!        rout(5)=sqrt(sum((cart(:,2)-cart(:,4))**2))
!        rout(6)=sqrt(sum((cart(:,3)-cart(:,4))**2))

        radtoang=180.d0/dacos(-1.d0)

        vec1(:)=(cart(:,2)-cart(:,1))*0.529177d0
        vec2(:)=(cart(:,3)-cart(:,1))*0.529177d0
        rout(4)=acos(sum(vec1*vec2)/(rout(1)*rout(2)))*radtoang

        vec1(:)=(cart(:,2)-cart(:,1))*0.529177d0
        vec2(:)=(cart(:,4)-cart(:,1))*0.529177d0
        rout(5)=acos(sum(vec1*vec2)/(rout(1)*rout(3)))*radtoang
        
        vec1(:)=(cart(:,3)-cart(:,1))*0.529177d0
        vec2(:)=(cart(:,4)-cart(:,1))*0.529177d0
        rout(6)=acos(sum(vec1*vec2)/(rout(2)*rout(3)))*radtoang


        return
        end subroutine cart_bond

	end module potdrivmd
