 program inc_gauss2tile

!---------------------------------------------------------------------
!
! Read a gaussian atmospheric increment file in netcdf.  Interpolate
! all fields to equivalent cubed sphere grid resolution.  Output the result
! in another netcdf file.
!
! Namelist variables:
! -------------------
! lev                   - Number of vertical levels.  Must be
!                         the same for the input and output grids.
! infile                - Path/name of input gaussian increment
!                         file (netcdf)
! outfile               - Path/name of output cubed sphere increment
!                         file (netcdf)
! gridpath              - Path to files oro_data.tile1.nc,...,oro_data.tile6.nc
!                         grid definition containing geolat/geolon
!
! 2023-07-07        Initial version.
!
!---------------------------------------------------------------------

 use netcdf
 use mpi

 implicit none

 integer, parameter :: num_recs = 9

 character(len=128) :: outfile, infile, gridpath, gridfile
 character(len=11)  :: records_in(num_recs), records_out(num_recs)
 character(len=1)   :: tilestr

 integer :: i, j, k, t
 integer :: lon_in, lat_in
 integer :: lev, ilev, lev_in
 integer :: ncid_in, id_var
 integer :: ncid_geo, res
 integer :: ncid_out, error
 integer :: dim_lon_out, dim_lat_out
 integer :: dim_lev_out, dim_ilev_out
 integer :: id_u_inc_out, id_v_inc_out
 integer :: id_lon_out, id_lat_out, id_lev_out
 integer :: id_pfull_out, id_ilev_out
 integer :: id_hyai_out, id_hybi_out
 integer :: id_delp_inc_out, id_delz_inc_out
 integer :: id_t_inc_out, id_sphum_inc_out
 integer :: id_liq_wat_inc_out, id_o3mr_inc_out
 integer :: id_icmr_inc_out, id_dim

 integer :: mpierr, mype, npes, mpistat(mpi_status_size)


 real, allocatable :: dummy_in(:,:,:)
 real, allocatable :: dummy_out(:,:,:,:)
 real, allocatable :: dummy_2d(:,:)

 real(8) :: rad2deg,dlondeg,deg2rad
 real(8), allocatable :: latitude_in(:), longitude_in(:)
 real(8), allocatable :: lat_out(:,:,:), lon_out(:,:,:)

 real, allocatable :: s2c(:,:,:,:)
 real, allocatable, dimension(:,:,:) :: id1, id2, jdc
 real, allocatable :: agrid(:,:,:,:)

 ! NOTE: u_inc,v_inc must be consecutive
 data records_in /'u_inc', 'v_inc', 'delp_inc', 'delz_inc', 'T_inc', &
                 'sphum_inc', 'liq_wat_inc', 'o3mr_inc', 'icmr_inc' /
 data records_out /'u', 'v', 'delp', 'delz', 'T', &
                  'sphum', 'liq_wat', 'o3mr', 'icmr' /

 namelist /setup/ outfile, infile, gridpath, lev


!-----------------------------------------------------------------
! MPI initialization
call mpi_init(mpierr)
call mpi_comm_rank(mpi_comm_world, mype, mpierr)
call mpi_comm_size(mpi_comm_world, npes, mpierr)
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Open and create output file records.  These will be filled
! with data below.
!-----------------------------------------------------------------

 if (mype == 0) print*,'- READ SETUP NAMELIST'
 open (43, file="./inc_gauss2tile.nml")
 read (43, nml=setup, iostat=error)
 if (error /= 0) then
   print*,"- FATAL ERROR READING NAMELIST. ISTAT IS: ", error
   stop 44
 endif
 close (43)

! Set constants
 rad2deg = 180.0_8 / (4.0_8 * atan(1.0_8))
 deg2rad = (4.0_8 * atan(1.0_8)) / 180.0_8

 ilev=lev+1

 call mpi_barrier(mpi_comm_world, mpierr)
! if (mype == 0) then
!   print*,'- OPEN OUTPUT FILE: ', trim(outfile)
!
!   error = nf90_create(outfile, cmode=IOR(NF90_CLOBBER,NF90_NETCDF4), ncid=ncid_out)
!   call netcdf_err(error, 'CREATING FILE='//trim(outfile) )
!
!   error = nf90_def_dim(ncid_out, 'lon', lon_out, dim_lon_out)
!   call netcdf_err(error, 'defining dimension lon for file='//trim(outfile) )
!
!   error = nf90_def_dim(ncid_out, 'lat', lat_out, dim_lat_out)
!   call netcdf_err(error, 'defining dimension lat for file='//trim(outfile) )
!
!   error = nf90_def_dim(ncid_out, 'lev', lev, dim_lev_out)
!   call netcdf_err(error, 'defining dimension lev for file='//trim(outfile) )
!
!   error = nf90_def_dim(ncid_out, 'ilev', ilev, dim_ilev_out)
!   call netcdf_err(error, 'defining dimension ilev for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'lon', nf90_double, (/dim_lon_out/), id_lon_out)
!   call netcdf_err(error, 'defining variable lon for file='//trim(outfile) )
!
!   error = nf90_put_att(ncid_out, id_lon_out, "units", "degrees_east")
!   call netcdf_err(error, 'define lon attribute for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'lat', nf90_double, (/dim_lat_out/), id_lat_out)
!   call netcdf_err(error, 'defining varable lat for file='//trim(outfile) )
!
!   error = nf90_put_att(ncid_out, id_lat_out, "units", "degrees_north")
!   call netcdf_err(error, 'defining lat att for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'lev', nf90_float, (/dim_lev_out/), id_lev_out)
!   call netcdf_err(error, 'defining variable lev for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'pfull', nf90_float, (/dim_lev_out/), id_pfull_out)
!   call netcdf_err(error, 'defining variable pfull for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'ilev', nf90_float, (/dim_ilev_out/), id_ilev_out)
!   call netcdf_err(error, 'defining variable ilev for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'hyai', nf90_float, (/dim_ilev_out/), id_hyai_out)
!   call netcdf_err(error, 'defining variable hyai for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'hybi', nf90_float, (/dim_ilev_out/), id_hybi_out)
!   call netcdf_err(error, 'defining variable hybi for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'u_inc', nf90_float, (/dim_lon_out,dim_lat_out,dim_lev_out/), id_u_inc_out)
!   call netcdf_err(error, 'defining variable u_inc for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'v_inc', nf90_float, (/dim_lon_out,dim_lat_out,dim_lev_out/), id_v_inc_out)
!   call netcdf_err(error, 'defining variable v_inc for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'delp_inc', nf90_float, (/dim_lon_out,dim_lat_out,dim_lev_out/), id_delp_inc_out)
!   call netcdf_err(error, 'defining variable delp_inc for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'delz_inc', nf90_float, (/dim_lon_out,dim_lat_out,dim_lev_out/), id_delz_inc_out)
!   call netcdf_err(error, 'defining variable delz_inc for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'T_inc', nf90_float, (/dim_lon_out,dim_lat_out,dim_lev_out/), id_t_inc_out)
!   call netcdf_err(error, 'defining variable t_inc for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'sphum_inc', nf90_float, (/dim_lon_out,dim_lat_out,dim_lev_out/), id_sphum_inc_out)
!   call netcdf_err(error, 'defining variable sphum_inc for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'liq_wat_inc', nf90_float, (/dim_lon_out,dim_lat_out,dim_lev_out/), id_liq_wat_inc_out)
!   call netcdf_err(error, 'defining variable liq_wat_inc for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'o3mr_inc', nf90_float, (/dim_lon_out,dim_lat_out,dim_lev_out/), id_o3mr_inc_out)
!   call netcdf_err(error, 'defining variable o3mr_inc for file='//trim(outfile) )
!
!   error = nf90_def_var(ncid_out, 'icmr_inc', nf90_float, (/dim_lon_out,dim_lat_out,dim_lev_out/), id_icmr_inc_out)
!   call netcdf_err(error, 'defining variable icmr_inc for file='//trim(outfile) )
!
!   error = nf90_put_att(ncid_out, nf90_global, 'source', 'GSI')
!   call netcdf_err(error, 'defining source attribute for file='//trim(outfile) )
!
!   error = nf90_put_att(ncid_out, nf90_global, 'comment', 'interpolated global analysis increment')
!   call netcdf_err(error, 'defining comment attribute for file='//trim(outfile) )
!
!   error = nf90_enddef(ncid_out, header_buffer_val, 4,0,4)
!   call netcdf_err(error, 'end meta define for file='//trim(outfile) )
! end if

!----------------------------------------------------
! Open and read input file dimensions and lat/lon arrays
!----------------------------------------------------

 if (mype == 0) print*,'- OPEN INPUT FILE: ', trim(infile)

 error = nf90_open(trim(infile), ior(nf90_nowrite, nf90_mpiio), &
                   comm=mpi_comm_world, info = mpi_info_null, ncid=ncid_in)
 call netcdf_err(error, 'opening file='//trim(infile) )

 error = nf90_inq_dimid(ncid_in, 'lon', id_dim)
 call netcdf_err(error, 'inquiring lon dimension for file='//trim(infile) )
 error = nf90_inquire_dimension(ncid_in, id_dim, len=lon_in)
 call netcdf_err(error, 'reading lon dimension for file='//trim(infile) )
 allocate(longitude_in(lon_in))
 error = nf90_inq_varid(ncid_in, 'lon', id_dim)
 call netcdf_err(error, 'inquiring var lon dimension for file='//trim(infile) )
 error = nf90_get_var(ncid_in, id_dim, longitude_in)
 call netcdf_err(error, 'reading longitude_in for file='//trim(infile) )

 error = nf90_inq_dimid(ncid_in, 'lat', id_dim)
 call netcdf_err(error, 'inquiring lat dimension for file='//trim(infile) )
 error = nf90_inquire_dimension(ncid_in, id_dim, len=lat_in)
 call netcdf_err(error, 'reading lat dimension for file='//trim(infile) )
 allocate(latitude_in(lat_in))
 error = nf90_inq_varid(ncid_in, 'lat', id_dim)
 call netcdf_err(error, 'inquiring var lat dimension for file='//trim(infile) )
 error = nf90_get_var(ncid_in, id_dim, latitude_in)
 call netcdf_err(error, 'reading latitude_in for file='//trim(infile) )

 error = nf90_inq_dimid(ncid_in, 'lev', id_dim)
 call netcdf_err(error, 'inquiring lev dimension for file='//trim(infile) )
 error = nf90_inquire_dimension(ncid_in, id_dim, len=lev_in)
 call netcdf_err(error, 'reading lev dimension for file='//trim(infile) )

!----------------------------------------------------
! Read in cubed sphere geometry
!----------------------------------------------------
! first allocate arrays based on the first tile
 t = 1
 write (tilestr, "(I1)") t
 gridfile=trim(gridpath)//'oro_data.tile'//trim(tilestr)//'.nc'
 error = nf90_open(trim(gridfile), ior(nf90_nowrite, nf90_mpiio), &
                   comm=mpi_comm_world, info = mpi_info_null, ncid=ncid_geo)
 call netcdf_err(error, 'opening file='//trim(gridfile) )
 error = nf90_inq_dimid(ncid_geo, 'lon', id_dim)
 call netcdf_err(error, 'inquiring lon dimension for file='//trim(gridfile) )
 error = nf90_inquire_dimension(ncid_geo, id_dim, len=res)
 call netcdf_err(error, 'getting case for file='//trim(gridfile) )
 print *, 'cubed sphere case is ',res
 allocate(lat_out(res,res,6),lon_out(res,res,6))
 allocate(dummy_2d(res,res))
 ! read in the geolat and geolon
 error = nf90_inq_varid(ncid_geo, 'geolat', id_dim)
 call netcdf_err(error, 'inquiring var geolat dimension for file='//trim(gridfile) )
 error = nf90_get_var(ncid_geo, id_dim, dummy_2d)
 call netcdf_err(error, 'reading geolat for file='//trim(gridfile) )
 lat_out(:,:,1) = dble(dummy_2d(:,:))
 error = nf90_inq_varid(ncid_geo, 'geolon', id_dim)
 call netcdf_err(error, 'inquiring var geolon dimension for file='//trim(gridfile) )
 error = nf90_get_var(ncid_geo, id_dim, dummy_2d)
 call netcdf_err(error, 'reading geolon for file='//trim(gridfile) )
 lon_out(:,:,1) = dble(dummy_2d(:,:))
 do t= 2, 6 ! loop through other tiles
   write (tilestr, "(I1)") t
   gridfile=trim(gridpath)//'oro_data.tile'//trim(tilestr)//'.nc'
   error = nf90_open(trim(gridfile), ior(nf90_nowrite, nf90_mpiio), &
                     comm=mpi_comm_world, info = mpi_info_null, ncid=ncid_geo)
   call netcdf_err(error, 'opening file='//trim(gridfile) )
   ! read in the geolat and geolon
   error = nf90_inq_varid(ncid_geo, 'geolat', id_dim)
   call netcdf_err(error, 'inquiring var geolat dimension for file='//trim(gridfile) )
   error = nf90_get_var(ncid_geo, id_dim, dummy_2d)
   call netcdf_err(error, 'reading geolat for file='//trim(gridfile) )
   lat_out(:,:,t) = dble(dummy_2d(:,:))
   error = nf90_inq_varid(ncid_geo, 'geolon', id_dim)
   call netcdf_err(error, 'inquiring var geolon dimension for file='//trim(gridfile) )
   error = nf90_get_var(ncid_geo, id_dim, dummy_2d)
   call netcdf_err(error, 'reading geolon for file='//trim(gridfile) )
   lon_out(:,:,t) = dble(dummy_2d(:,:))
 end do
!----------------------------------------------------
! Calculate interpolation weights
!----------------------------------------------------
do i=1,lon_in
  longitude_in(i) = longitude_in(i) * deg2rad
enddo
do i=1,lat_in
  latitude_in(i) = latitude_in(i) * deg2rad
enddo
lat_out(:,:,:) = lat_out(:,:,:) * deg2rad
lon_out(:,:,:) = lon_out(:,:,:) * deg2rad

allocate(s2c(res, res, 4, 6))
allocate(id1(res, res, 6))
allocate(id2(res, res, 6))
allocate(jdc(res, res, 6))
allocate(agrid(res+1, res+1, 2, 6))

do k=1,6
  do i=1,res
    do j=1,res
      agrid(i,j,1,k) = lon_out(i,j,k)
      agrid(i,j,2,k) = lat_out(i,j,k)
    enddo
  enddo
  call remap_coef( 1, res, 1, res, 1, res+1, 1, res+1, &
                   lon_in, lat_in, longitude_in, latitude_in, id1(:,:,k), &
                   id2(:,:,k), jdc(:,:,k), s2c(:,:,:,k), agrid(:,:,:,k))
enddo

!----------------------------------------------------
! Open output file for writing
!----------------------------------------------------

!----------------------------------------------------
! Loop through fields, interpolate, and write out
!----------------------------------------------------

 error = nf90_close(ncid_in)
 call mpi_barrier(mpi_comm_world, mpierr)

! if (mype == 0) then
!   print*,"- WRITE OUTPUT FILE: ", trim(outfile)
!
!  ! lev
!
!   allocate(levs(lev))
!   do j = 1, lev
!     levs(j) = j
!   enddo
!
!   error = nf90_put_var(ncid_out, id_lev_out, levs)
!   call netcdf_err(error, 'writing levs for file='//trim(outfile) )
!
!  ! pfull
!
!   error = nf90_put_var(ncid_out, id_pfull_out, levs)
!   call netcdf_err(error, 'writing pfull for file='//trim(outfile) )
!
!   deallocate (levs)
!   allocate (levs(ilev))
!   do j = 1, ilev
!     levs(j) = j
!   enddo
!
!  ! ilev
!
!   error = nf90_put_var(ncid_out, id_ilev_out, levs)
!   call netcdf_err(error, 'writing ilev for file='//trim(outfile) )
!
!  ! hyai
!
!   error = nf90_put_var(ncid_out, id_hyai_out, levs)
!   call netcdf_err(error, 'writing hyai for file='//trim(outfile) )
!
!  ! hybi
!
!   error = nf90_put_var(ncid_out, id_hybi_out, levs)
!   call netcdf_err(error, 'writing hybi for file='//trim(outfile) )
!
!  ! latitude
!
!   error = nf90_put_var(ncid_out, id_lat_out, latitude_out)
!   call netcdf_err(error, 'writing latitude for file='//trim(outfile) )
!
!  ! longitude
!
!   error = nf90_put_var(ncid_out, id_lon_out, longitude_out)
!   call netcdf_err(error, 'writing longitude for file='//trim(outfile) )
!
!   deallocate(levs)
!
!   error = nf90_close(ncid_out)
! end if

 call mpi_barrier(mpi_comm_world, mpierr)
 if (mype == 0) print*,'- NORMAL TERMINATION'

 call mpi_barrier(mpi_comm_world, mpierr)
 call mpi_finalize(mpierr)

 end program inc_gauss2tile

 subroutine netcdf_err( err, string )

 use netcdf

 implicit none
 integer, intent(in) :: err
 character(len=*), intent(in) :: string
 character(len=256) :: errmsg

 if( err.EQ.NF90_NOERR )return
 errmsg = NF90_STRERROR(err)
 print*,''
 print*,'FATAL ERROR: ', trim(string), ': ', trim(errmsg)
 print*,'STOP.'
 stop 999

 return
 end subroutine netcdf_err

  subroutine remap_coef( is, ie, js, je, isd, ied, jsd, jed, &
      im, jm, lon, lat, id1, id2, jdc, s2c, agrid )

    integer, intent(in):: is, ie, js, je, isd, ied, jsd, jed
    integer, intent(in):: im, jm
    real,    intent(in):: lon(im), lat(jm)
    real,    intent(out):: s2c(is:ie,js:je,4)
    integer, intent(out), dimension(is:ie,js:je):: id1, id2, jdc
    real,    intent(in):: agrid(isd:ied,jsd:jed,2)
    ! local:
    real :: rdlon(im)
    real :: rdlat(jm)
    real:: a1, b1
    integer i,j, i1, i2, jc, i0, j0
    do i=1,im-1
      rdlon(i) = 1. / (lon(i+1) - lon(i))
    enddo
    rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

    do j=1,jm-1
      rdlat(j) = 1. / (lat(j+1) - lat(j))
    enddo

    ! * Interpolate to cubed sphere cell center
    do 5000 j=js,je

      do i=is,ie

        if ( agrid(i,j,1)>lon(im) ) then
          i1 = im;     i2 = 1
          a1 = (agrid(i,j,1)-lon(im)) * rdlon(im)
        elseif ( agrid(i,j,1)<lon(1) ) then
          i1 = im;     i2 = 1
          a1 = (agrid(i,j,1)+2.*pi-lon(im)) * rdlon(im)
        else
          do i0=1,im-1
            if ( agrid(i,j,1)>=lon(i0) .and. agrid(i,j,1)<=lon(i0+1) ) then
              i1 = i0;  i2 = i0+1
              a1 = (agrid(i,j,1)-lon(i1)) * rdlon(i0)
              go to 111
            endif
          enddo
        endif
111     continue

        if ( agrid(i,j,2)<lat(1) ) then
          jc = 1
          b1 = 0.
        elseif ( agrid(i,j,2)>lat(jm) ) then
          jc = jm-1
          b1 = 1.
        else
          do j0=1,jm-1
            if ( agrid(i,j,2)>=lat(j0) .and. agrid(i,j,2)<=lat(j0+1) ) then
              jc = j0
              b1 = (agrid(i,j,2)-lat(jc)) * rdlat(jc)
              go to 222
            endif
          enddo
        endif
222     continue

        !if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
        !     write(*,*) 'gid=', i,j,a1, b1
        !endif

        s2c(i,j,1) = (1.-a1) * (1.-b1)
        s2c(i,j,2) =     a1  * (1.-b1)
        s2c(i,j,3) =     a1  *     b1
        s2c(i,j,4) = (1.-a1) *     b1
        id1(i,j) = i1
        id2(i,j) = i2
        jdc(i,j) = jc
      enddo   !i-loop
5000 continue   ! j-loop

  end subroutine remap_coef
