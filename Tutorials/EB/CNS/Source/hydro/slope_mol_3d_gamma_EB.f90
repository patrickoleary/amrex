module ebslope_module
  
  implicit none

  private

  public ebslopex, ebslopey, ebslopez

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ebslopex(q,qd_lo,qd_hi, &
                        dqxal,qpd_lo,qpd_hi, &
                        flag, fglo, fghi, &
                        ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,nv)

      use mempool_module, only : bl_allocate, bl_deallocate
      use meth_params_module
      use bl_constants_module

      implicit none

      integer          :: qd_lo(3), qd_hi(3)
      integer          :: qpd_lo(3),qpd_hi(3)
      integer, intent(in) :: fglo(3),fghi(3)
      integer          :: ilo1, ilo2, ihi1, ihi2, ilo3, ihi3, nv

      double precision :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),nv)
      double precision :: dqxal(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),nv)
      integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fglo(3))


      integer i, j, k, n

      double precision  slop, dq1
      double precision dm, dp, dc, ds, sl, dl, dfm, dfp

      integer ilo, ihi      
      
      double precision :: dsgn,dlim,dcen
      double precision, pointer :: dlft(:,:),drgt(:,:)

      integer :: nbr(-1:1,-1:1,-1:1)

      call bl_allocate (dlft( ilo1-2:ihi1+2,nv))
      call bl_allocate (drgt( ilo1-2:ihi1+2,nv))

      if(plm_iorder.eq.1) then

         do n = 1, nv
            do k = ilo3-1, ihi3+1
            do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1
                  dqxal(i,j,k,n) = ZERO
               enddo
            enddo
            enddo
         enddo

      else


            ! Compute slopes in first coordinate direction
         do k = ilo3-1, ihi3+1
            do j = ilo2-1, ihi2+1

               do i = ilo1-1, ihi1+1

                  call get_neighbors( flag(i,j,k), nbr )

                  if( nbr(0,0,0) .eq. 0 .or. nbr(-1,0,0) .eq. 0) then
                     dlft(i,:) = 0.d0
                  else

                     dlft(i,1) = 0.5d0*(q(i,j,k,QP)-q(i-1,j,k,QP))/qaux(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU))
                     dlft(i,2) = (q(i,j,k,QRHO)-q(i-1,j,k,QRHO))- (q(i,j,k,QP) - q(i-1,j,k,QP))/qaux(i,j,k,QC)**2
                     dlft(i,3) = 0.5d0*(q(i,j,k,QP)-q(i-1,j,k,QP))/qaux(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU))
                     dlft(i,4) = q(i,j,k,QV) - q(i-1,j,k,QV) 
                     dlft(i,5) = q(i,j,k,QW) - q(i-1,j,k,QW) 
                  endif

                  if( nbr(0,0,0) .eq. 0 .or. nbr(1,0,0) .eq. 0) then
                     drgt(i,:) = 0.d0
                  else

                     drgt(i,1) = 0.5d0*(q(i+1,j,k,QP)-q(i,j,k,QP))/qaux(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU))
                     drgt(i,2) = (q(i+1,j,k,QRHO)-q(i,j,QRHO))- (q(i,j,k,QP) - q(i,j,k,QP))/qaux(i,j,k,QC)**2
                     drgt(i,3) = 0.5d0*(q(i+1,j,k,QP)-q(i,j,k,QP))/qaux(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU))
                     drgt(i,4) = q(i+1,j,k,QV) - q(i,j,k,QV) 
                     drgt(i,5) = q(i+1,j,k,QW) - q(i,j,k,QW) 
                  endif

               enddo

               do n=1,nv

               ! First compute Fromm slopes
               do i = ilo1-1, ihi1+1
                  dcen = 0.5d0 * (dlft(i,n)+drgt(i,n))
                  dsgn = sign(ONE, dcen)
                  slop =2.d0* min( abs(dlft(i,n)), abs(drgt(i,n)) )
                  if (dlft(i,n)*drgt(i,n) .ge. ZERO) then
                     dlim = slop
                  else
                     dlim = ZERO
                  endif
                  dqxal(i,j,k,n) = dsgn*min( dlim, abs(dcen) )
               enddo


            enddo
         enddo

      call bl_deallocate (dlft)
      call bl_deallocate (drgt)

      end subroutine ebslopex

      subroutine ebslopey(q,qd_lo,qd_hi, &
                        dqyal,qpd_lo,qpd_hi, &
                        flag, fglo, fghi, &
                        ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,nv)

      use mempool_module, only : bl_allocate, bl_deallocate
      use meth_params_module
      use bl_constants_module

      implicit none

      integer          :: qd_lo(3), qd_hi(3)
      integer          :: qpd_lo(3),qpd_hi(3)
      integer, intent(in) :: fglo(3),fghi(3)
      integer          :: ilo1, ilo2, ihi1, ihi2, ilo3, ihi3, nv

      double precision :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),nv)
      double precision :: dqyal(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),nv)
      integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fglo(3))

      integer i, j, k, n

      double precision slop
      double precision dm, dp, dc, ds, sl, dl, dfm, dfp

      integer ilo, ihi      
      integer :: nbr(-1:1,-1:1,-1:1)
      
      double precision::dsgn,dlim,dcen
      double precision, pointer::dlft(:,:),drgt(:,:)

      call bl_allocate (dlft( ilo1-2:ihi1+2,nv))
      call bl_allocate (drgt( ilo1-2:ihi1+2,nv))


      if(plm_iorder.eq.1) then

         do n = 1, nv
            do k = ilo3-1, ihi3+1
            do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1
                  dqyal(i,j,k,n) = ZERO
               enddo
            enddo
            enddo
         enddo

      else


            ! Compute slopes in first coordinate direction
         do k = ilo3-1, ihi3+1
            do j = ilo2-1, ihi2+1

               do i = ilo1-1, ihi1+1

                  call get_neighbors( flag(i,j,k), nbr )

                  if( nbr(0,0,0) .eq. 0 .or. nbr(0,-1,0) .eq. 0) then
                     dlft(i,:) = 0.d0
                  else

                     dlft(i,1) = 0.5d0*(q(i,j,k,QP)-q(i,j-1,k,QP))/qaux(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV))
                     dlft(i,2) = (q(i,j,k,QRHO)-q(i,j-1,k,QRHO))- (q(i,j,k,QP) - q(i,j-1,k,QP))/qaux(i,j,k,QC)**2
                     dlft(i,3) = 0.5d0*(q(i,j,k,QP)-q(i,j-1,k,QP))/qaux(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV))
                     dlft(i,4) = q(i,j,k,QU) - q(i,j-1,k,QU) 
                     dlft(i,5) = q(i,j,k,QW) - q(i,j-1,k,QW) 
                  endif

                  if( nbr(0,0,0) .eq. 0 .or. nbr(0,1,0) .eq. 0) then
                     drgt(i,:) = 0.d0
                  else

                     drgt(i,1) = 0.5d0*(q(i,j+1,k,QP)-q(i,j,k,QP))/qaux(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV))
                     drgt(i,2) = (q(i,j+1,k,QRHO)-q(i,j,QRHO))- (q(i,j+1,k,QP) - q(i,j,k,QP))/qaux(i,j,k,QC)**2
                     drgt(i,3) = 0.5d0*(q(i+1,j,k,QP)-q(i,j,k,QP))/qaux(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV))
                     drgt(i,4) = q(i,j+1,k,QU) - q(i,j,k,QU) 
                     drgt(i,5) = q(i,j+1,k,QW) - q(i,j,k,QW) 
                   endif

               enddo

               do n=1,nv

               ! First compute Fromm slopes
               do i = ilo1-1, ihi1+1
                  dcen = 0.5d0 * (dlft(i,n)+drgt(i,n))
                  dsgn = sign(ONE, dcen)
                  slop = 2.d0* min( abs(dlft(i,n)), abs(drgt(i,n)) )
                  if (dlft(i,n)*drgt(i,n) .ge. ZERO) then
                     dlim = slop
                  else
                     dlim = ZERO
                  endif
                  dqyal(i,j,k,n) = dsgn*min( dlim, abs(dcen) )
               enddo

               enddo
            enddo
          enddo

      call bl_deallocate (dlft)
      call bl_deallocate (drgt)

      end subroutine ebslopey

      subroutine ebslopez(q,qd_lo,qd_hi, &
                        dqyal,qpd_lo,qpd_hi, &
                        flag, fglo, fghi, &
                        ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,nv)

      use mempool_module, only : bl_allocate, bl_deallocate
      use meth_params_module
      use bl_constants_module

      implicit none

      integer          :: qd_lo(3), qd_hi(3)
      integer          :: qpd_lo(3),qpd_hi(3)
      integer, intent(in) :: fglo(3),fghi(3)
      integer          :: ilo1, ilo2, ihi1, ihi2, ilo3, ihi3, nv

      double precision :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),nv)
      double precision :: dqzal(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),nv)
      integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fglo(3))

      integer i, j, k, n
      integer :: nbr(-1:1,-1:1,-1:1)

      double precision slop, dq1
      double precision dm, dp, dc, ds, sl, dl, dfm, dfp

      integer ilo, ihi      
      
      double precision::dsgn,dlim,dcen
      double precision, pointer::dlft(:,:),drgt(:,:)

      call bl_allocate (dlft( ilo1-2:ihi1+2,nv))
      call bl_allocate (drgt( ilo1-2:ihi1+2,nv))


      if(plm_iorder.eq.1) then

         do n = 1, nv
            do k = ilo3-1, ihi3+1
            do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1
                  dqzal(i,j,k,n) = ZERO
               enddo
            enddo
            enddo
         enddo

      else


            ! Compute slopes in first coordinate direction
         do k = ilo3-1, ihi3+1
            do j = ilo2-1, ihi2+1

               do i = ilo1-1, ihi1+1


                  call get_neighbors( flag(i,j,k), nbr )

                  if( nbr(0,0,0) .eq. 0 .or. nbr(0,0,-1) .eq. 0) then
                     dlft(i,:) = 0.d0
                  else

                     dlft(i,1) = 0.5d0*(q(i,j,k,QP)-q(i,j,k-1,QP))/qaux(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW))
                     dlft(i,2) = (q(i,j,k,QRHO)-q(i,j,k-1,QRHO))- (q(i,j,k,QP) - q(i,j,k-1,QP))/qaux(i,j,k,QC)**2
                     dlft(i,3) = 0.5d0*(q(i,j,k,QP)-q(i,j,k-1,QP))/qaux(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW))
                     dlft(i,4) = q(i,j,k,QU) - q(i,j,k-1,QU) 
                     dlft(i,5) = q(i,j,k,QV) - q(i,j,k-1,QV) 
                  endif

                  if( nbr(0,0,0) .eq. 0 .or. nbr(0,0,1) .eq. 0) then
                     drgt(i,:) = 0.d0
                  else

                     drgt(i,1) = 0.5d0*(q(i,j,k+1,QP)-q(i,j,k,QP))/qaux(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW))
                     drgt(i,2) = (q(i,j,k+1,QRHO)-q(i,j,QRHO))- (q(i,j+1,k,QP) - q(i,j,k,QP))/qaux(i,j,k,QC)**2
                     drgt(i,3) = 0.5d0*(q(i,j,k+1,QP)-q(i,j,k,QP))/qaux(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW))
                     drgt(i,4) = q(i,j,k+1,QU) - q(i,j,k,QU) 
                     drgt(i,5) = q(i,j,k+1,QV) - q(i,j,k,QV) 
                   endif

               enddo

               do n=1,nv

               ! First compute Fromm slopes
               do i = ilo1-1, ihi1+1
                  dcen = 0.5d0 * (dlft(i,n)+drgt(i,n))
                  dsgn = sign(ONE, dcen)
                  slop = 2.d0*min( abs(dlft(i,n), abs(drgt(i,n)) )
                  if (dlft(i,n)*drgt(i,n) .ge. ZERO) then
                     dlim = slop
                  else
                     dlim = ZERO
                  endif
                  dqzal(i,n) = dsgn*min( dlim, abs(dcen) )
               enddo

               enddo
            enddo
          enddo


      call bl_deallocate (dlft)
      call bl_deallocate (drgt)

      end subroutine ebslopez

end module ebslope_module