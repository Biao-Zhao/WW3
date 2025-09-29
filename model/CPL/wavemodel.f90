module wave_model_mod

  use mpp_mod, only: mpp_npes, mpp_pe, mpp_get_current_pelist
  use mpp_domains_mod, only: mpp_global_field, mpp_get_compute_domain
  use mpp_domains_mod, only: mpp_define_domains, mpp_define_layout

  use wave_type_mod, only: wave_data_type, atmos_wave_boundary_type, ice_wave_boundary_type

  use time_manager_mod, only: time_type, operator(+), get_date

  USE WMMDATMD, ONLY: MDSI, MDSO, MDSS, MDST, MDSE, &
                      NMPROC, IMPROC, NMPSCR, NRGRD, ETIME


  implicit none

  INTEGER, ALLOCATABLE :: TEND(:,:), TSTRT(:,:)
  INTEGER              :: MPI_COMM = -99

  contains


    !> wave_model_init initializes the wave model interface
    subroutine wave_model_init(Atm2Waves,Ice2Waves,Wav)
      ! Author: Brandon Reichl
      ! Origination: August 2019

      ! USE statements
      use wminitmd, only: wminit, wminitnml
      use w3gdatmd, only: NX, NY
      use w3gdatmd, only: usspf, ussp_wn

      ! Subroutine arguments
      type(atmos_wave_boundary_type), intent(inout)  :: Atm2Waves
      type(ice_wave_boundary_type), intent(inout)  :: Ice2Waves
      type(wave_data_type), intent(inout)            :: Wav

      ! Local parameters
      logical              :: Flag_NML
      INTEGER              :: Ierr_MPI
      integer :: layout(2)

      !----------------------------------------------------------------------
      !This block of code checks if the wave model has been initialized
      if (Wav%Waves_Is_Init) then
        write(*,*)'wave_model_init was called, but it is already registered'
        stop
      end if
      Wav%Waves_Is_Init = .true.
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      !This block of code sets up the MPI communicator for WW3
      ! The equivalent calls in ww3_multi are given in comments
      ! -> FMS has already done the MPI_INIT.
      ! CALL MPI_INIT      ( IERR_MPI )
      ! MPI_COMM = MPI_COMM_WORLD
      ! -> We need to access the commID from FMS.  This is the only
      !    way I can figure out using MPP libraries.
      call mpp_get_current_pelist(Wav%pelist,commID=MPI_COMM)
      ! -> FMS does provide an easy way to get the comm_size and rank
      ! CALL MPI_COMM_SIZE ( MPI_COMM, NMPROC, IERR_MPI )
      !CALL MPI_COMM_RANK ( MPI_COMM, IMPROC, IERR_MPI )
      NMPROC = mpp_npes()
      IMPROC = mpp_pe()
      IMPROC = IMPROC + 1
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      !This block of code calls the WW3 intializers
      INQUIRE(FILE="ww3_multi.nml", EXIST=Flag_NML)
      IF (Flag_NML) THEN
        CALL WMINITNML ( MDSI, MDSO, MDSS, 10, MDSE, 'ww3_multi.nml', MPI_COMM )
      ELSE
        CALL WMINIT ( MDSI, MDSO, MDSS, 10, MDSE, 'ww3_multi.inp', MPI_COMM )
      END IF
      ALLOCATE ( TEND(2,NRGRD), TSTRT(2,NRGRD) )
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      !This block of codes sets up the wave domains
      call mpp_define_layout((/1,NX,1,NY/),mpp_npes(),layout)
      call mpp_define_domains((/1,NX,1,NY/),layout, Wav%domain )
      !----------------------------------------------------------------------

      !
      !This block of code sets up the global coupler domains
      allocate(Atm2Waves%wavgrd_u10_glo(1:NX,1:NY,1))
      Atm2Waves%wavgrd_u10_glo(1:NX,1:NY,1) = 0.0
      allocate(Atm2Waves%wavgrd_v10_glo(1:NX,1:NY,1))
      Atm2Waves%wavgrd_v10_glo(1:NX,1:NY,1) = 0.0
      allocate(Ice2Waves%wavgrd_ucurr_glo(1:NX,1:NY,1))
      Ice2Waves%wavgrd_ucurr_glo(1:NX,1:NY,1) = 0.0
      allocate(Ice2Waves%wavgrd_vcurr_glo(1:NX,1:NY,1))
      Ice2Waves%wavgrd_vcurr_glo(1:NX,1:NY,1) = 0.0
      allocate(Ice2Waves%wavgrd_Ice_glo(1:NX,1:NY,1))
      Ice2Waves%wavgrd_Ice_glo(1:NX,1:NY,1) = 0.0

      if (usspf(1).gt.usspf(2)) then
        print*,'usspf(1): ',usspf(1)
        print*,'usspf(2): ',usspf(2)
        print*,'usspf(1) is greater than usspf(2) in mod_def.ww3. ',&
             ' Check ww3_grid.inp file, recreate mod_def.ww3 with ussp <= iussp to continue'
        stop
      endif
      Wav%num_Stk_bands = usspf(2)-usspf(1)+1
      allocate(Wav%stk_wavenumbers(1:Wav%num_Stk_Bands)) ; Wav%stk_wavenumbers(:) = 0.0
      Wav%stk_wavenumbers(:) = USSP_WN(usspf(1):usspf(2))

      allocate(Wav%ustkb_glo(1:NX,1:NY,Wav%num_stk_bands))
      allocate(Wav%vstkb_glo(1:NX,1:NY,Wav%num_stk_bands))
      allocate(Wav%hs_glo(1:NX,1:NY,1))
      allocate(Wav%ust_glo(1:NX,1:NY,1))
      allocate(Wav%ustdir_glo(1:NX,1:NY,1))
      allocate(Wav%charn_glo(1:NX,1:NY,1))
      allocate(Wav%tauox_glo(1:NX,1:NY,1))
      allocate(Wav%tauoy_glo(1:NX,1:NY,1))
      Wav%ustkb_glo(:,:,:)  = 0.0
      Wav%vstkb_glo(:,:,:)  = 0.0
      Wav%hs_glo(:,:,:)     = 0.0
      Wav%ust_glo(:,:,:)    = 0.0
      Wav%ustdir_glo(:,:,:) = 0.0
      Wav%charn_glo(:,:,:)  = 0.0
      Wav%tauox_glo(:,:,:)  = 0.0
      Wav%tauoy_glo(:,:,:)  = 0.0
      return
    end subroutine wave_model_init

    !> wave_model_timestep integrates the wave fields over one
    !!  coupling timestep
    subroutine update_wave_model(Atm2Waves, Ice2Waves, Wav, Time_start, Time_increment)
      ! Author: Brandon Reichl
      ! Origination: August 2019

      ! USE statements
      use wmwavemd, only: wmwave
      use w3gdatmd, only: NX, NY
      use w3idatmd, only: wx0, wxN, wy0, wyN, TW0, TWN, &
                          cx0, cxN, cy0, cyN, ICEI, TC0, TCN
      use w3adatmd, only: ussx, ussy, ussp, &
                          charn, hs, tauox, tauoy
      use w3wdatmd, only: UST, USTDIR
      use w3gdatmd, only: nseal, mapsf, NK, MAPSTA, MAPST2, &
                          gtype, ungtype, iclose, iclose_none
      use w3odatmd, only: iaproc, naproc
      ! Subroutine arguments
      type(atmos_wave_boundary_type), intent(in) :: Atm2Waves
      type(ice_wave_boundary_type),   intent(in) :: Ice2Waves
      type(wave_data_type),        intent(inout) :: Wav
      type(time_type),                intent(in) :: Time_start,&
                                                    Time_increment

      ! Local parameters
      integer :: I, yr, mo, da, hr, mi, se, is, ie, js, je
      integer :: isea, isea_g, ix, iy, jx, jy, pix, piy
      integer :: b

      integer :: glob_loc_x(NX,NY), glob_loc_y(NX,NY)
      real    :: mapout(NX,NY,1)
      logical :: is_west, is_east, is_south, is_north
      !----------------------------------------------------------------------
      !Convert the ending time of this call into WW3 time format, which
      ! is integer(2) :: (YYYYMMDD, HHMMSS)
      DO I=1, NRGRD
        call get_date(Time_start,&
             yr,mo,da,hr,mi,se)
        TSTRT(1,I) = yr*1e4+mo*1e2+da
        TSTRT(2,I) = hr*1e4+mi*1e2+se
        call get_date(Time_start+Time_increment,&
             yr,mo,da,hr,mi,se)
        TEND(1,I) = yr*1e4+mo*1e2+da
        TEND(2,I) = hr*1e4+mi*1e2+se
      END DO
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      !Call WW3 timestepper with an argument for the time to return
      ! back
      TW0(:) = TSTRT(:,1)
      TWN(:) = TEND(:,1)
      TC0(:) = TSTRT(:,1)
      TCN(:) = TEND(:,1)
      call mpp_global_field(Wav%domain,atm2waves%wavgrd_u10_mpp(:,:,:),atm2waves%wavgrd_u10_glo(:,:,:))
      wx0(:,:) = Atm2Waves%wavgrd_u10_glo(:,:,1)
      wxN(:,:) = wx0(:,:)
      call mpp_global_field(Wav%domain,atm2waves%wavgrd_v10_mpp(:,:,:),atm2waves%wavgrd_v10_glo(:,:,:))
      wy0(:,:) = Atm2Waves%wavgrd_v10_glo(:,:,1)
      wyN(:,:) = wy0(:,:)
      call mpp_global_field(Wav%domain,ice2waves%wavgrd_ucurr_mpp(:,:,:),ice2waves%wavgrd_ucurr_glo(:,:,:))
      cx0(:,:) = Ice2Waves%wavgrd_ucurr_glo(:,:,1)
      cxN(:,:) = cx0(:,:)
      call mpp_global_field(Wav%domain,ice2waves%wavgrd_vcurr_mpp(:,:,:),ice2waves%wavgrd_vcurr_glo(:,:,:))
      cy0(:,:) = Ice2Waves%wavgrd_vcurr_glo(:,:,1)
      cyN(:,:) = cy0(:,:)

      call mpp_global_field(Wav%domain,ice2waves%wavgrd_Ice_mpp(:,:,:),ice2waves%wavgrd_Ice_glo(:,:,:))
      ICEI(:,:) = Ice2Waves%wavgrd_Ice_glo(:,:,1)
      ! write(*,*)'Into wave model U10 max: ',maxval(atm2Waves%U_10_global),maxval(wx0)
      ! write(*,*)'Into wave model V10 max: ',maxval(atm2Waves%V_10_global),maxval(wy0)
      ! write(*,*)'Into wave model UO max: ',maxval(ice2Waves%Ucurr_global),maxval(cx0)
      ! write(*,*)'Into wave model VO max: ',maxval(ice2Waves%Vcurr_global),maxval(cy0)
      CALL WMWAVE ( TEND )

      Wav%glob_loc_X(:,:) = 0
      Wav%glob_loc_Y(:,:) = 0
      
      call mpp_get_compute_domain( Wav%domain, is, ie, js, je )
      isea = 1
      do ix=is,ie
        do iy=js,je
          if (isea<=nseal) then
            ISEA_G   = IAPROC + (ISEA-1)*NAPROC
            jx = MAPSF(ISEA_G,1)
            jy = MAPSF(ISEA_G,2)
            Wav%glob_loc_X(ix,iy) = jx
            Wav%glob_loc_Y(ix,iy) = jy
            do b=1,Wav%num_stk_bands
              Wav%ustkb_mpp(ix,iy,b) = USSP(isea,b)
              Wav%vstkb_mpp(ix,iy,b) = USSP(isea,NK+b)
            enddo
            ! send wave variables to coupler and atmospheric model, added by Biao
            Wav%landmask(ix,iy,1)    = MAPSTA(JY,JX)
            Wav%hs(ix,iy,1)          = HS(isea)
            Wav%ust_wav(ix,iy,1)     = UST(isea_g)
            Wav%ustdir_wav(ix,iy,1)  = USTDIR(isea_g)
            Wav%charn_wav(ix,iy,1)   = CHARN(isea)
            Wav%tauox_wav(ix,iy,1)   = TAUOX(isea)
            Wav%tauoy_wav(ix,iy,1)   = TAUOY(isea)
          endif
          isea = isea+1
        end do
      end do

      do b = 1,Wav%num_stk_bands
        call mpp_global_field(Wav%domain,Wav%ustkb_mpp(:,:,b),wav%ustkb_glo(:,:,b))
        call mpp_global_field(Wav%domain,Wav%vstkb_mpp(:,:,b),wav%vstkb_glo(:,:,b))
      enddo

      call mpp_global_field(Wav%domain,Wav%hs,wav%hs_glo)
      call mpp_global_field(Wav%domain,Wav%ust_wav,wav%ust_glo)
      call mpp_global_field(Wav%domain,Wav%ustdir_wav,wav%ustdir_glo)
      call mpp_global_field(Wav%domain,Wav%charn_wav,wav%charn_glo)
      call mpp_global_field(Wav%domain,Wav%tauox_wav,wav%tauox_glo)
      call mpp_global_field(Wav%domain,Wav%tauoy_wav,wav%tauoy_glo)

      call mpp_global_field(Wav%domain,Wav%glob_loc_X,glob_loc_X)
      call mpp_global_field(Wav%domain,Wav%glob_loc_Y,glob_loc_Y)
      call mpp_global_field(Wav%domain,Wav%landmask,mapout)

      Wav%ustkb_mpp(:,:,:)  = 0.0
      Wav%vstkb_mpp(:,:,:)  = 0.0
      Wav%hs(:,:,:)         = 0.0
      Wav%ust_wav(:,:,:)    = 0.0
      Wav%ustdir_wav(:,:,:) = 0.0
      Wav%charn_wav(:,:,:)  = 0.0
      Wav%tauox_wav(:,:,:)  = 0.0
      Wav%tauoy_wav(:,:,:)  = 0.0
      Wav%landmask(:,:,:)   = 3.0
      do ix=1,NX
         do iy=1,NY
           Pix = glob_loc_X(ix,iy)
           Piy = glob_loc_Y(ix,iy)
           if (Pix>=is .and. Pix<=ie .and. Piy>=js .and. Piy<=je) then
              do b = 1,Wav%num_stk_bands
                Wav%ustkb_mpp(Pix,Piy,b) = wav%ustkb_glo(ix,iy,b)
                Wav%vstkb_mpp(Pix,Piy,b) = wav%vstkb_glo(ix,iy,b)
              enddo
                Wav%hs(Pix,Piy,1) = wav%hs_glo(ix,iy,1)
                Wav%ust_wav(Pix,Piy,1) = wav%ust_glo(ix,iy,1)
                Wav%ustdir_wav(Pix,Piy,1) = wav%ustdir_glo(ix,iy,1)
                Wav%charn_wav(Pix,Piy,1) = wav%charn_glo(ix,iy,1)
                Wav%tauox_wav(Pix,Piy,1) = wav%tauox_glo(ix,iy,1)
                Wav%tauoy_wav(Pix,Piy,1) = wav%tauoy_glo(ix,iy,1)
                Wav%landmask(Pix,Piy,1)  = mapout(ix,iy,1)
           endif
         enddo
      enddo
     
     ! Added by Biao,  WW3 flags the points at four boudaries as excluded, so those cells stay 0.
     ! Before coupling to MOM6/SHiELD, call fill_boundary to copy the nearest
     ! interior row/column to these four lines, preventing zeros on domain edges. only work for 
     ! structured and none closure grid setup
     if (gtype.ne.ungtype .and. iclose .eq. iclose_none ) then
        is_west  = (is == 1)
        is_east  = (ie == NX)
        is_south = (js == 1)
        is_north = (je == NY)
        call fill_boundary(Wav%ustkb_mpp, Wav%landmask, is_west, is_east, is_south, is_north)
        call fill_boundary(Wav%vstkb_mpp, Wav%landmask, is_west, is_east, is_south, is_north) 
        call fill_boundary(Wav%hs,        Wav%landmask, is_west, is_east, is_south, is_north)
        call fill_boundary(Wav%ust_wav,   Wav%landmask, is_west, is_east, is_south, is_north)
        call fill_boundary(Wav%ustdir_wav,Wav%landmask, is_west, is_east, is_south, is_north)
        call fill_boundary(Wav%charn_wav, Wav%landmask, is_west, is_east, is_south, is_north)
        call fill_boundary(Wav%tauox_wav, Wav%landmask, is_west, is_east, is_south, is_north)
        call fill_boundary(Wav%tauoy_wav, Wav%landmask, is_west, is_east, is_south, is_north)
      end if
      !----------------------------------------------------------------------
      return
    end subroutine update_wave_model

    !>This subroutine finishes the wave routine and deallocates memory
    subroutine wave_model_end(Waves, Atm2Waves,Ice2Waves)
      ! Use statements
      use wmfinlmd, only: wmfinl

      ! Subroutine arguments
      type(wave_data_type), intent(inout) :: Waves
      type(atmos_wave_boundary_type), intent(inout) :: Atm2Waves
      type(ice_wave_boundary_type), intent(inout) :: Ice2Waves

      ! Local parameters
      INTEGER              :: IERR_MPI
      !----------------------------------------------------------------------
      !Finalize the driver
      CALL WMFINL
      DEALLOCATE ( TEND, TSTRT )
      deallocate(Atm2Waves%wavgrd_u10_glo)
      deallocate(Atm2Waves%wavgrd_v10_glo)
      deallocate(Ice2Waves%wavgrd_ucurr_glo)
      deallocate(Ice2Waves%wavgrd_vcurr_glo)
      deallocate(Ice2Waves%wavgrd_Ice_glo)
      deallocate(Waves%ustkb_glo)
      deallocate(Waves%vstkb_glo)
      deallocate(Waves%hs_glo)
      deallocate(Waves%ust_glo)
      deallocate(Waves%ustdir_glo)
      deallocate(Waves%charn_glo)
      deallocate(Waves%tauox_glo)
      deallocate(Waves%tauoy_glo)
      CALL MPI_BARRIER ( MPI_COMM, IERR_MPI ) !Do we need this?
      !----------------------------------------------------------------------

      return
    end subroutine wave_model_end
   
    !> This subroutine copies the nearest interior row/column onto the outer lines, added by Biao 
    subroutine fill_boundary(field, landmask, is_west, is_east, is_south, is_north)
      implicit none
      real, intent(inout) :: field(:,:,:)
      real, intent(in)    :: landmask(:,:,:)
      logical, intent(in) :: is_west, is_east, is_south, is_north
      !   local variables
      integer :: nxl,nyl,nzl,i,j

      nxl = size(field,1)
      nyl = size(field,2)
      nzl = size(field,3)
      
  if (is_west .and. nxl>=2) then
    do j=1,nyl
      if (landmask(1,j,1)==3 .and. landmask(2,j,1)==1) field(1,j,:) = field(2,j,:)
    end do
  end if

  if (is_east .and. nxl>=2) then
    do j=1,nyl
      if (landmask(nxl,j,1)==3 .and. landmask(nxl-1,j,1)==1) field(nxl,j,:) = field(nxl-1,j,:)
    end do
  end if

  if (is_south .and. nyl>=2) then
    do i=1,nxl
      if (landmask(i,1,1)==3 .and. landmask(i,2,1)==1) field(i,1,:) = field(i,2,:)
    end do
  end if

  if (is_north .and. nyl>=2) then
    do i=1,nxl
      if (landmask(i,nyl,1)==3 .and. landmask(i,nyl-1,1)==1) field(i,nyl,:) = field(i,nyl-1,:)
    end do
  end if 

      return
   end subroutine fill_boundary

end module wave_model_mod
