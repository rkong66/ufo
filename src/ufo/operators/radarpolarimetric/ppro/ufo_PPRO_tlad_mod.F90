! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_PPRO_tlad_mod

  use oops_variables_mod
  use obs_variables_mod
  use ufo_vars_mod
  use ufo_geovals_mod
  use vert_interp_mod
  use missing_values_mod
  use dualpol_op_mod
  use zhang21_tlad_mod
  use fckit_log_module, only: fckit_log
  use ufo_PPRO_mod, only: nopvar

! ------------------------------------------------------------------------------

  type, public :: ufo_PPRO_tlad
  private
    type(obs_variables), public :: obsvars
    type(oops_variables), public :: geovars
    integer, public :: obsvarindices(nopvar) ! Array that maps indices of opvars to indices of hofx columns
                                             ! -1 means not used
    character(len=MAXVARLEN), public :: v_coord ! GeoVaL to use to interpolate in vertical
    character(len=MAXVARLEN), public :: micro_option    ! Choice (enum) of microphysics option
    integer :: nval, nlocs
    real(kind_real), allocatable :: wf(:)
    integer, allocatable :: wi(:)
    integer, allocatable :: iband(:)

    ! trajectory trays, see trajectory variables in dualpol_op_tlad_mod.f90
    real(kind_real), allocatable :: qr_tray(:), qs_tray(:), qg_tray(:), qh_tray(:)
    real(kind_real), allocatable :: ntr_tray(:), nts_tray(:), ntg_tray(:), nth_tray(:)

    real(kind_real), allocatable :: qpr_tray(:), qms_tray(:), qmg_tray(:), qps_tray(:), qpg_tray(:), qmh_tray(:), qph_tray(:)
    real(kind_real), allocatable :: ntpr_tray(:), ntms_tray(:), ntmg_tray(:), ntps_tray(:), ntpg_tray(:), ntmh_tray(:), ntph_tray(:)

    real(kind_real), allocatable :: wr_tray(:), ws_tray(:), wg_tray(:), wps_tray(:), wpg_tray(:), wh_tray(:), wph_tray(:)
    real(kind_real), allocatable :: dmr_tray(:), dms_tray(:), dmg_tray(:), dmps_tray(:), dmpg_tray(:), dmh_tray(:), dmph_tray(:) 
    real(kind_real), allocatable :: zr_tray(:), zs_tray(:), zg_tray(:), zps_tray(:), zpg_tray(:), zmh_tray(:), zph_tray(:)
    real(kind_real), allocatable :: dens_tray(:), deng_tray(:), denh_tray(:)
    real(kind_real), allocatable :: rats_tray(:), ratg_tray(:), rath_tray(:)

    real(kind_real), allocatable :: zhr_tray(:), zhs_tray(:), zhg_tray(:), zhps_tray(:), zhpg_tray(:), zhh_tray(:), zhph_tray(:)
    real(kind_real), allocatable :: zdrr_tray(:), zdrs_tray(:), zdrg_tray(:), zdrps_tray(:), zdrpg_tray(:), zdrh_tray(:), zdrph_tray(:)
    real(kind_real), allocatable :: zbarr_tray(:), zbars_tray(:), zbarg_tray(:), zbarps_tray(:), zbarpg_tray(:), zbarh_tray(:), zbarph_tray(:)
    real(kind_real), allocatable :: zvr_tray(:), zvs_tray(:), zvg_tray(:), zvps_tray(:), zvpg_tray(:), zvh_tray(:), zvph_tray(:)
    real(kind_real), allocatable :: kdpr_tray(:), kdps_tray(:), kdpg_tray(:), kdpps_tray(:), kdppg_tray(:), kdph_tray(:), kdpph_tray(:)
    real(kind_real), allocatable :: phvr_tray(:), phvs_tray(:), phvg_tray(:), phvps_tray(:), phvpg_tray(:), phvh_tray(:), phvph_tray(:)

    real(kind_real), allocatable :: as_tray(:,:), ag_tray(:,:), ah_tray(:,:)

    real(kind_real), allocatable :: density_air_tray(:), temp_air_tray(:)
    real(kind_real), allocatable :: zh_tray(:), zdr_tray(:), kdp_tray(:), phv_tray(:)

    logical, allocatable :: badbg_tray(:)

  contains
    procedure :: setup => PPRO_tlad_setup_
    procedure :: cleanup => PPRO_tlad_cleanup_
    procedure :: settraj => PPRO_tlad_settraj_
    procedure :: simobs_tl => PPRO_simobs_tl_
    procedure :: simobs_ad => PPRO_simobs_ad_
    final :: destructor
  end type ufo_PPRO_tlad

  integer :: n_geovars
  character(len=maxvarlen), dimension(:), allocatable :: geovars_list
  character(len=maxvarlen):: varname_qr, varname_qs, varname_qg, varname_qh, varname_qnr, varname_qns, varname_qng, varname_qnh

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine PPRO_tlad_setup_(self, yaml_conf)
  use fckit_configuration_module, only: fckit_configuration
  use iso_c_binding
  use ufo_PPRO_mod, only: nopvar, opvars_list
  implicit none
  class(ufo_PPRO_tlad), intent(inout) :: self
  type(fckit_configuration), intent(in) :: yaml_conf

  character(kind=c_char,len=:), allocatable :: coord_name
  character(kind=c_char,len=:), allocatable :: micro_option
  character(kind=c_char,len=:), allocatable :: this_varname
  character(len=maxvarlen) :: var_string
  character(len=512) :: buffer
  integer :: i, j
  logical :: found

  ! Initialize all as not used variables
  self%obsvarindices(:) = -1

  DO i = 1, self%obsvars%nvars()
    write(buffer,*) 'PPRO TLAD obsvars(', i, ') = ', TRIM(self%obsvars%variable(i))
    call fckit_log%info(buffer); buffer = ''
  ENDDO

  DO i = 1, self%obsvars%nvars() ! loop of simulated variables
    found = .False.
    DO j = 1, nopvar ! loop of operator variables
      IF (TRIM(self%obsvars%variable(i)) .EQ. TRIM(opvars_list(j))) THEN
        self%obsvarindices(j) = i
        found = .True.
        EXIT
      ENDIF
    ENDDO
    IF (.not. found) THEN
      print*, ' Simulated variable: ', trim(self%obsvars%variable(i)), ' not found'
      call abor1_ftn("PPRO: simulated variable not found in operator, aborting")
    ENDIF
  ENDDO

  write(buffer,*) 'PPRO TLAD obsvarindices = ', self%obsvarindices(:)
  call fckit_log%info(buffer); buffer = ''

  ! YAML option for microphysics scheme
  call yaml_conf%get_or_die("microphysics option", micro_option)
  if(trim(micro_option) .eq. "WSM6") then
    self%micro_option = 'WSM6'
    n_geovars=6
  else if(trim(micro_option) .eq. "Thompson") then
    self%micro_option = 'THOMPSON'
    n_geovars=7
  else if(trim(micro_option) .eq. "NSSL") then
    self%micro_option = 'NSSL'
    n_geovars=11
  else
    print*, 'TLAD micro_option set to: ', trim(micro_option)
    call abor1_ftn("microphysics option not set or unsupported, aborting")
  endif

  if ( .not. allocated(geovars_list) ) allocate(geovars_list(n_geovars))
  !geovars_list(1) = var_airdens
  geovars_list(1) = var_prs
  geovars_list(2) = var_ts
  geovars_list(3) = var_q

  var_string="var_rain_mixing_ratio"
  if( yaml_conf%has(trim(var_string)) ) then
    call yaml_conf%get_or_die(trim(var_string), this_varname)
    geovars_list(4) = this_varname
    varname_qr = this_varname
  endif
  var_string="var_snow_mixing_ratio"
  if( yaml_conf%has(trim(var_string)) ) then
    call yaml_conf%get_or_die(trim(var_string), this_varname)
    geovars_list(5) = this_varname
    varname_qs = this_varname
  endif
  var_string="var_graupel_mixing_ratio"
  if( yaml_conf%has(trim(var_string)) ) then
    call yaml_conf%get_or_die(trim(var_string), this_varname)
    geovars_list(6) = this_varname
    varname_qg = this_varname
  endif
  
  if(trim(self%micro_option) .eq. "THOMPSON") then 
    var_string="var_rain_number_concentration"
    if( yaml_conf%has(trim(var_string)) ) then
      call yaml_conf%get_or_die(trim(var_string), this_varname)
      geovars_list(7) = this_varname
      varname_qnr = this_varname
    endif
  endif

  if(trim(self%micro_option) .eq. "NSSL") then
    var_string="var_hail_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
      call yaml_conf%get_or_die(trim(var_string), this_varname)
      geovars_list(7) = this_varname
      varname_qh = this_varname
    endif
    var_string="var_rain_number_concentration"
    if( yaml_conf%has(trim(var_string)) ) then
      call yaml_conf%get_or_die(trim(var_string), this_varname)
      geovars_list(8) = this_varname
      varname_qnr = this_varname
    endif
    var_string="var_snow_number_concentration"
    if( yaml_conf%has(trim(var_string)) ) then
      call yaml_conf%get_or_die(trim(var_string), this_varname)
      geovars_list(9) = this_varname
      varname_qns = this_varname
    endif
    var_string="var_grauple_number_concentration"
    if( yaml_conf%has(trim(var_string)) ) then
      call yaml_conf%get_or_die(trim(var_string), this_varname)
      geovars_list(10) = this_varname
      varname_qng = this_varname
    endif
    var_string="var_hail_number_concentration"
    if( yaml_conf%has(trim(var_string)) ) then
      call yaml_conf%get_or_die(trim(var_string), this_varname)
      geovars_list(11) = this_varname
      varname_qnh = this_varname
    endif
  endif

  ! YAML option for vertical coordinate name
  if( yaml_conf%has("VertCoord") ) then
    call yaml_conf%get_or_die("VertCoord",coord_name)
    self%v_coord = coord_name
    if( trim(self%v_coord) .ne. var_z .and. trim(self%v_coord) .ne. var_zm .and. &
        trim(self%v_coord) .ne. var_geomz ) then
        call abor1_ftn("ufo_PPRO: incorrect vertical coordinate specified")
    endif
  else ! default
    self%v_coord = var_z
  endif

  call self%geovars%push_back(geovars_list)
  call self%geovars%push_back(self%v_coord)

end subroutine PPRO_tlad_setup_

! ------------------------------------------------------------------------------

subroutine PPRO_tlad_settraj_(self, geovals, obss)
  use obsspace_mod
  use ufo_PPRO_mod, only: ufo_PPRO_sim1obs
  implicit none
  class(ufo_PPRO_tlad), intent(inout) :: self
  type(ufo_geovals),         intent(in)    :: geovals
  type(c_ptr), value,        intent(in)    :: obss

  !local variables
  logical :: variable_present
  integer :: iobs, ivar, nvars_geovars
  real(kind_real),  dimension(:), allocatable :: obsvcoord
  type(ufo_geoval), pointer :: vcoordprofile, profile  

  character(len=MAXVARLEN) :: geovar

  real(kind_real), allocatable :: tmp(:)
  real(kind_real) :: tmp2

  real(kind_real), allocatable :: fields(:,:)  ! background fields interplated vertically to obs height

  integer:: i
  integer:: nvars = 1 ! nvars = 1 for radar reflectivity only
  real (kind=8) :: rho ! kg/m^3
  real (kind=8) :: t ! temperature, Kelvin
  real (kind=8) :: qr, qs, qg, qh ! mixing ratio, g/kg
  real (kind=8) :: nr, ns, ng, nh
  real (kind=8) :: zh, zdr, kdp, phv

  real (kind=8) :: p, q 
  real(kind_real), parameter :: rd=287.04_kind_real
  real(kind_real), parameter :: one=1._kind_real
  real(kind_real), parameter :: D608=0.608_kind_real

  ! Make sure nothing already allocated
  call self%cleanup()

  ! Get height profiles from geovals
  call ufo_geovals_get_var(geovals, self%v_coord, vcoordprofile)
  self%nval = vcoordprofile%nval

  ! Get the observation vertical coordinates
  self%nlocs = obsspace_get_nlocs(obss)
  allocate(obsvcoord(self%nlocs))
  allocate(self%iband(self%nlocs))

  !call obsspace_get_db(obss, "MetaData", "geometric_height", obsvcoord)
  call obsspace_get_db(obss, "MetaData", "height", obsvcoord)
  if (variable_present) then
     call obsspace_get_db(obss, "MetaData", "iband", self%iband)
  else
     self%iband(:) = 1
  endif

  ! Allocate arrays for interpolation weights
  allocate(self%wi(self%nlocs))
  allocate(self%wf(self%nlocs))

  ! Calculate the interpolation weights
  allocate(tmp(vcoordprofile%nval))
  do iobs = 1, self%nlocs
    tmp = vcoordprofile%vals(:,iobs)
    tmp2 = obsvcoord(iobs)
    call vert_interp_weights(vcoordprofile%nval, tmp2, tmp, self%wi(iobs), self%wf(iobs))
  enddo

  ! Number of variables in geovars (without the vertical coordinate)
  nvars_geovars = self%geovars%nvars() - 1
  allocate(fields(nvars_geovars,self%nlocs))
  do ivar = 1, nvars_geovars
    ! Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)
    ! Get profiles for this variable from geovals, profiles at nlocs
    call ufo_geovals_get_var(geovals, geovar, profile)
    ! Interpolate from geovals to observational location into hofx
    do iobs = 1, self%nlocs
      call vert_interp_apply(profile%nval, profile%vals(:,iobs), &
                             & fields(ivar,iobs), self%wi(iobs), self%wf(iobs))
    enddo
  enddo

  call allocate_traj_trays(self)

  call init_coefs ! read all S-band and C-band coeffs from text files

  do iobs=1,self%nlocs

    p  = fields(1,iobs) ! pressure (Pa)
    t  = fields(2,iobs) ! temperature (K)
    q  = fields(3,iobs) ! specific humidity
    q = q/(one-q) !convert from specific humidity to mixing ratio
    rho = p/(rd*t*(one+D608*q))  !calculate dry air density, kg/m^3
    !rho  = fields(1,iobs) ! moist air density, kg/m^3
    !t    = fields(2,iobs) ! Kelvin
    qr   = 1000.0*fields(3,iobs) ! kg/kg -> g/kg
    qs   = 1000.0*fields(4,iobs) ! kg/kg -> g/kg
    qg   = 1000.0*fields(5,iobs) ! kg/kg -> g/kg
    if ( trim(self%micro_option) .eq. "THOMPSON" ) then
       nr = fields(6,iobs)
    else if ( trim(self%micro_option) .eq. "NSSL" ) then
       qh = 1000.0*fields(6,iobs)
       nr = 1000.0*fields(7,iobs)
       ns = 1000.0*fields(8,iobs)
       ng = 1000.0*fields(9,iobs)
       nh = 1000.0*fields(10,iobs)
    endif

    if (rho .LT. 0.0) then
      self%badbg_tray(iobs) = .True.
      cycle
    endif

    if ( trim(self%micro_option) .eq. "WSM6" ) then
        call ufo_PPRO_sim1obs(self%iband(iobs), rho, t, qr, qs, qg, zh, zdr, kdp, phv, &
             set_trajectory=.True.)
    else if ( trim(self%micro_option) .eq. "THOMPSON" ) then
        call ufo_PPRO_sim1obs(self%iband(iobs), rho, t, qr, qs, qg, zh, zdr, kdp, phv, &
             nr=nr, set_trajectory=.True.)
    else if ( trim(self%micro_option) .eq. "NSSL" ) then
        call ufo_PPRO_sim1obs(self%iband(iobs), rho, t, qr, qs, qg, zh, zdr, kdp, phv, &
             qh=qh, nr=nr, ns=ns, ng=ng, nh=nh, set_trajectory=.True.)
    endif
    
    ! save trajectory variables to trajectory trays
    call store_traj(self, iobs)
  enddo

! Cleanup memory
  deallocate(obsvcoord)
  deallocate(tmp)

end subroutine PPRO_tlad_settraj_

! ------------------------------------------------------------------------------

subroutine PPRO_simobs_tl_(self, geovals, obss, nvars, nlocs, hofx)
  implicit none
  class(ufo_PPRO_tlad), intent(in) :: self
  type(ufo_geovals),         intent(in) :: geovals
  integer,                   intent(in) :: nvars, nlocs
  real(c_double),         intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value,        intent(in) :: obss

  integer :: iobs, ivar, nvars_geovars
  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar

  real(kind_real), allocatable :: fields(:,:)  ! increment fields interplated vertically to obs height
  real(kind_real) :: qr_tl, qs_tl, qg_tl, qh_tl
  real(kind_real) :: nr_tl, ns_tl, ng_tl, nh_tl
  real(kind_real) :: zh_tl, zdr_tl, kdp_tl, phv_tl
  real(kind_real) :: zh_dBZ_tl, zdr_dB_tl
  integer       :: obsvarindex

  real(c_double) :: missing

! Number of variables in geovars (without the vertical coordinate)
  nvars_geovars = self%geovars%nvars() - 1
  allocate(fields(nvars_geovars,nlocs))
  fields=0.0

  do ivar = 1, nvars_geovars
    geovar = self%geovars%variable(ivar)

! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

! Interpolate from geovals to observational location into hofx
    do iobs = 1, nlocs
      call vert_interp_apply_tl(profile%nval, profile%vals(:,iobs), &
                             & fields(ivar,iobs), self%wi(iobs), self%wf(iobs))
    enddo
  enddo

  missing = missing_value(missing)

  do iobs=1,nlocs

    if (self%badbg_tray(iobs)) then
      hofx(:,iobs) = missing
      cycle
    endif

    qr_tl   = 1000.0*fields(3,iobs) ! kg/kg -> g/kg
    qs_tl   = 1000.0*fields(4,iobs) ! kg/kg -> g/kg
    qg_tl   = 1000.0*fields(5,iobs) ! kg/kg -> g/kg
    if ( trim(self%micro_option) .eq. "THOMPSON" ) then
       nr_tl = fields(6,iobs)
    else if ( trim(self%micro_option) .eq. "NSSL" ) then
       qh_tl = 1000.0*fields(6,iobs)
       nr_tl = 1000.0*fields(7,iobs)
       ns_tl = 1000.0*fields(8,iobs)
       ng_tl = 1000.0*fields(9,iobs)
       nh_tl = 1000.0*fields(10,iobs)
    endif

    ! restore trajectory variables from trajectory trays
    call restore_traj(self, iobs)

    if ( trim(self%micro_option) .eq. "WSM6" ) then
       call ufo_PPRO_sim1obs_tl(self%iband(iobs), & 
          qr_tl, qs_tl, qg_tl, &
          zh_tl, zdr_tl, kdp_tl, phv_tl) 
    else if ( trim(self%micro_option) .eq. "THOMPSON" ) then
       call ufo_PPRO_sim1obs_tl(self%iband(iobs), & 
          qr_tl, qs_tl, qg_tl, &
          zh_tl, zdr_tl, kdp_tl, phv_tl, &
          ntr_tl=nr_tl)
    else if ( trim(self%micro_option) .eq. "NSSL" ) then
       call ufo_PPRO_sim1obs_tl(self%iband(iobs), & 
          qr_tl, qs_tl, qg_tl, &
          zh_tl, zdr_tl, kdp_tl, phv_tl, &
          qh_tl=qh_tl, ntr_tl=nr_tl, nts_tl=ns_tl, ntg_tl=ng_tl, nth_tl=nh_tl)
    endif

    ! tangent linear of ZH = 10LOG10(zh)
    zh_dBZ_tl = 10.0 / LOG(10.0) * zh_tl / zh_traj

    ! tangent linear of ZDR = 10LOG10(zdr)
    zdr_dB_tl = 10.0 / LOG(10.0) * zdr_tl / zdr_traj

    do ivar=1, nopvar
      obsvarindex = self%obsvarindices(ivar)
      if (obsvarindex .EQ. -1) CYCLE ! VARIABLE NOT USED
      select case(ivar)
        case (1) ! zh
          hofx(obsvarindex,iobs) = zh_dBZ_tl
        case (2) ! zdr
          hofx(obsvarindex,iobs) = zdr_dB_tl
        case (3) ! kdp
          hofx(obsvarindex,iobs) = kdp_tl
        case (4) ! phv
          hofx(obsvarindex,iobs) = phv_tl
      end select
    enddo
  enddo

  deallocate(fields)

end subroutine PPRO_simobs_tl_

! ------------------------------------------------------------------------------

subroutine PPRO_simobs_ad_(self, geovals, obss, nvars, nlocs, hofx)
  implicit none
  class(ufo_PPRO_tlad), intent(in) :: self
  type(ufo_geovals),         intent(inout) :: geovals
  integer,                   intent(in)    :: nvars, nlocs
  real(c_double),            intent(in)    :: hofx(nvars, nlocs)
  type(c_ptr), value,        intent(in)    :: obss

  integer :: iobs, ivar, nvars_geovars
  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar
  real(c_double) :: missing
  real(kind_real), allocatable :: fields(:,:)  ! increment fields interplated vertically to obs height

  real(kind_real) :: qr_ad, qs_ad, qg_ad, qh_ad
  real(kind_real) :: nr_ad, ns_ad, ng_ad, nh_ad
  real(kind_real) :: zh_ad, zdr_ad, kdp_ad, phv_ad
  real(kind_real) :: zh_dBZ_ad, zdr_dB_ad
  integer       :: obsvarindex

  character(len=512) :: buffer

  missing = missing_value(missing)

! Number of variables in geovars (without the vertical coordinate)
  nvars_geovars = self%geovars%nvars() - 1
  allocate(fields(nvars_geovars,nlocs))

  fields=0.0
  do iobs=1,nlocs

    if (self%badbg_tray(iobs)) then
      cycle
    endif
    
    ! restore trajectory variables from trajectory trays
    call restore_traj(self, iobs)

    ! Initialize adjoint values
    zh_dBZ_ad=0.0;zdr_dB_ad=0.0;kdp_ad=0.0;phv_ad=0.0
    do ivar=1, nopvar
      obsvarindex = self%obsvarindices(ivar)
      if (obsvarindex .EQ. -1) CYCLE ! VARIABLE NOT USED
      if (hofx(obsvarindex,iobs) .EQ. missing) CYCLE
      select case(ivar)
        case (1) ! zh
          zh_dBZ_ad = hofx(obsvarindex,iobs) 
        case (2) ! zdr
          zdr_dB_ad = hofx(obsvarindex,iobs)
        case (3) ! kdp
          kdp_ad = hofx(obsvarindex,iobs) 
        case (4) ! phv
          phv_ad = hofx(obsvarindex,iobs) 
      end select
    enddo

    ! adjoint of ZH = 10LOG10(zh)
    zh_ad = 10.0 / LOG(10.0) / zh_traj * zh_dBZ_ad 

    ! adjoint of ZDR = 10LOG10(zdr)
    zdr_ad = 10.0 / LOG(10.0) / zdr_traj * zdr_dB_ad         
    
    qr_ad=0.0;qs_ad=0.0;qg_ad=0.0
    if ( trim(self%micro_option) .eq. "WSM6" ) then
       call ufo_PPRO_sim1obs_ad(self%iband(iobs), & 
           qr_ad, qs_ad, qg_ad, &
           zh_ad, zdr_ad, kdp_ad, phv_ad)
    else if ( trim(self%micro_option) .eq. "THOMPSON" ) then
        nr_ad=0.0
        call ufo_PPRO_sim1obs_ad(self%iband(iobs), &
           qr_ad, qs_ad, qg_ad, &
           zh_ad, zdr_ad, kdp_ad, phv_ad, &
           ntr_ad=nr_ad)
    else if ( trim(self%micro_option) .eq. "NSSL" ) then
        qh_ad=0.0;nr_ad=0.0;ns_ad=0.0;ng_ad=0.0;nh_ad=0.0
        call ufo_PPRO_sim1obs_ad(self%iband(iobs), &
           qr_ad, qs_ad, qg_ad, &
           zh_ad, zdr_ad, kdp_ad, phv_ad, &
           qh_ad=qh_ad, ntr_ad=nr_ad, nts_ad=ns_ad, ntg_ad=ng_ad, nth_ad=nh_ad)
    endif

    if (ISNAN(qr_ad) .OR. ISNAN(qs_ad) .OR. ISNAN(qg_ad)) then
      call print_traj
      buffer = ''
      write(buffer,*) 'zh_ad = ', zh_ad, ', qr_ad = ', qr_ad, ', qs_ad = ', qs_ad, 'qg_ad = ', qg_ad
      call fckit_log%info(buffer)
      STOP
    endif

    if (trim(self%micro_option) .eq. "THOMPSON") then
       if (ISNAN(nr_ad)) then
          call print_traj
          buffer = ''
          write(buffer,*) 'Thompson scheme: zh_ad = ', zh_ad, ', nr_ad = ', nr_ad
          call fckit_log%info(buffer)
          STOP
       endif
    endif

    if (trim(self%micro_option) .eq. "NSSL") then
       if (ISNAN(qh_ad) .or. ISNAN(nr_ad) .or. ISNAN(ns_ad) .or. ISNAN(ng_ad) .or. ISNAN(nh_ad)) then
          call print_traj
          buffer = ''
          write(buffer,*) 'NSSL scheme: zh_ad = ', zh_ad, ', qh_ad = ', qh_ad, &
                  ', nr_ad = ', nr_ad, ', ns_ad = ', ns_ad, ', ng_ad = ', ng_ad, ', nh_ad = ', nh_ad
          call fckit_log%info(buffer)
          STOP
       endif
    endif

    ! adjoint of unit conversion g/kg -> kg/kg
    fields(3,iobs) = fields(3,iobs) + qr_ad*1000.0
    fields(4,iobs) = fields(4,iobs) + qs_ad*1000.0
    fields(5,iobs) = fields(5,iobs) + qg_ad*1000.0
    if (trim(self%micro_option) .eq. "THOMPSON") then
      fields(6,iobs) = fields(6,iobs) + nr_ad
    endif
    if (trim(self%micro_option) .eq. "NSSL") then
      fields(6,iobs) = fields(6,iobs) + qh_ad*1000.0
      fields(7,iobs) = fields(7,iobs) + nr_ad
      fields(8,iobs) = fields(8,iobs) + ns_ad
      fields(9,iobs) = fields(9,iobs) + ng_ad
      fields(10,iobs) = fields(10,iobs) + nh_ad
    endif
    
  enddo

  do ivar = 1, nvars_geovars
! Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)

! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

! Interpolate from geovals to observational location into hofx
    do iobs = 1, nlocs
      call vert_interp_apply_ad(profile%nval, profile%vals(:,iobs), &
                             & fields(ivar,iobs), self%wi(iobs), self%wf(iobs))
    enddo
  enddo

  deallocate(fields)
end subroutine PPRO_simobs_ad_

! ------------------------------------------------------------------------------

subroutine PPRO_tlad_cleanup_(self)
  implicit none
  class(ufo_PPRO_tlad), intent(inout) :: self
  self%nval = 0
  self%nlocs = 0
  if (allocated(self%wi)) deallocate(self%wi)
  if (allocated(self%wf)) deallocate(self%wf)
  if (allocated(self%iband)) deallocate(self%iband)
  call deallocate_traj_trays(self)
end subroutine PPRO_tlad_cleanup_

! ------------------------------------------------------------------------------

subroutine  destructor(self)
  type(ufo_PPRO_tlad), intent(inout)  :: self

  call self%cleanup()

end subroutine destructor

subroutine ufo_PPRO_sim1obs_tl(iband, & 
  qr_tl, qs_tl, qg_tl, &
  zh_tl, zdr_tl, kdp_tl, phv_tl, &
  qh_tl, ntr_tl, nts_tl, ntg_tl, nth_tl)

  integer, intent(in) :: iband
  real(kind=8), intent(in) :: qr_tl, qs_tl, qg_tl
  real(kind=8), intent(in), optional :: qh_tl, ntr_tl, nts_tl, ntg_tl, nth_tl
  real(kind=8), intent(out) :: zh_tl, zdr_tl, kdp_tl, phv_tl
  
  real (kind=8), dimension(:,:),  pointer :: rain_coefs => NULL()
  real (kind=8), dimension(:,:),  pointer :: snow_coefs => NULL()
  real (kind=8), dimension(:,:),  pointer :: graupel_coefs => NULL()
  real (kind=8), dimension(:,:),  pointer :: hail_coefs => NULL()
  real (kind=8), dimension(:),    pointer :: snow_a => NULL()
  real (kind=8), dimension(:),    pointer :: graupel_a => NULL()
  real (kind=8), dimension(:),    pointer :: hail_a => NULL()

  integer :: i
  
  ! output of melting scheme
  real (kind=8) :: qpr_tl, qms_tl, qmg_tl, qps_tl, qpg_tl, qmh_tl, qph_tl
  real (kind=8) :: ntpr_tl, ntms_tl, ntmg_tl, ntps_tl, ntpg_tl, ntmh_tl, ntph_tl

  ! "r" is "pr", pure rain
  ! "s" is "ms", melting snow
  ! "g" is "mg", melting graupel
  real (kind=8) :: zhr_tl, zdrr_tl, kdpr_tl, phvr_tl
  real (kind=8) :: zhs_tl, zdrs_tl, kdps_tl, phvs_tl
  real (kind=8) :: zhg_tl, zdrg_tl, kdpg_tl, phvg_tl
  real (kind=8) :: zhps_tl, zdrps_tl, kdpps_tl, phvps_tl
  real (kind=8) :: zhpg_tl, zdrpg_tl, kdppg_tl, phvpg_tl
  real (kind=8) :: zhh_tl, zdrh_tl, kdph_tl, phvh_tl
  real (kind=8) :: zhph_tl, zdrph_tl, kdpph_tl, phvph_tl

  real (kind=8) :: wr_tl, ws_tl, wg_tl, wps_tl, wpg_tl, wh_tl, wph_tl
  real (kind=8) :: zr_tl, zs_tl, zg_tl, zps_tl, zpg_tl, zmh_tl, zph_tl
  real (kind=8) :: dmr_tl, dms_tl, dmg_tl, dmps_tl, dmpg_tl, dmh_tl, dmph_tl
  real (kind=8) :: dens_tl, deng_tl, denh_tl
  real (kind=8) :: rats_tl, ratg_tl, rath_tl
  real (kind=8) :: as_tl(0:12), ag_tl(0:12), ah_tl(0:12)
  real (kind=8) :: aps_tl(0:12), apg_tl(0:12), aph_tl(0:12)
  real (kind=8) :: denr_tl, denps_tl, denpg_tl,denph_tl

  aps_tl(:) = 0.0_8
  apg_tl(:) = 0.0_8
  denr_tl = 0.0_8
  denps_tl = 0.0_8
  denpg_tl = 0.0_8
  aph_tl(:) = 0.0_8
  denph_tl = 0.0_8

  if (iband == 1) then ! s-band radar
      rain_coefs => sband_rain_coefs
      snow_coefs => sband_snow_coefs
      graupel_coefs => sband_graupel_coefs
      hail_coefs => sband_hail_coefs
      snow_a => sband_snow_a
      graupel_a => sband_graupel_a
      hail_a => sband_hail_a
  else ! c-band radar
      rain_coefs => cband_rain_coefs
      snow_coefs => cband_snow_coefs
      graupel_coefs => cband_graupel_coefs
      hail_coefs => cband_hail_coefs
      snow_a => cband_snow_a
      graupel_a => cband_graupel_a
      hail_a => cband_hail_a
  endif

  stop  !Rong Kong temprarily added here
  ! melting scheme
  if (present(qh_tl) .and. present(ntr_tl) .and. present(nts_tl) .and. present(ntg_tl) .and. present(nth_tl)) then
     call melting_scheme_zhang24_tl( &
          qr_tl, qs_tl, qg_tl, qms_tl, qmg_tl, qpr_tl, qps_tl, qpg_tl, rats_tl, ratg_tl, &
          ntr_tl=ntr_tl, nts_tl=nts_tl, ntg_tl=ntg_tl, &
          ntms_tl=ntms_tl, ntmg_tl=ntmg_tl, ntpr_tl=ntpr_tl, ntps_tl=ntps_tl, ntpg_tl=ntps_tl, &
          qh_tl=qh_tl, nth_tl=nth_tl, rath_tl=rath_tl, &
          qmh_tl=qmh_tl, qph_tl=qph_tl, ntmh_tl=ntmh_tl, ntph_tl=ntph_tl )
  else !! WSM6 and Thompson are the same
     call melting_scheme_zhang24_tl( &
          qr_tl, qs_tl, qg_tl, qms_tl, qmg_tl, qpr_tl, qps_tl, qpg_tl, rats_tl, ratg_tl )
  endif

  ! 1. pure rain
  call set_traj_precip_type('rain')
  call watercontent_tl(qpr_tl, wr_tl)
  if (present(ntr_tl)) then
    call dm_z_2moment_ad(qpr_tl, ntpr_tl, denr_tl, dmr_tl, zr_tl)
  else
    call dm_z_wsm6_tl('rain', qpr_tl, denr_tl, dmr_tl, zr_tl)
  endif
  call dualpol_op_rain_tl(rain_coefs, wr_tl, dmr_tl, zhr_tl, zdrr_tl, kdpr_tl, phvr_tl)
  
  ! 2. melting snow
  call set_traj_precip_type('snow')
  call watercontent_tl(qms_tl, ws_tl)
  call weticephase_density_tl('snow', rats_tl, dens_tl)
  if (present(nts_tl)) then
    call dm_z_2moment_ad(qms_tl, ntms_tl, dens_tl, dms_tl, zs_tl)
  else
    call dm_z_wsm6_tl('snow', qms_tl, dens_tl, dms_tl, zs_tl)
  endif
  do i = 0, 12
    call coef_a_tl(snow_coefs(i,0:3), rats_tl, as_tl(i))
  end do
  call dualpol_op_icephase_tl(zs_tl, ws_tl, dens_tl, dms_tl, as_tl, &
                              zhs_tl, zdrs_tl, kdps_tl, phvs_tl)

  ! 3. melting graupel
  call set_traj_precip_type('graupel')
  call watercontent_tl(qmg_tl, wg_tl)
  call weticephase_density_tl('graupel', ratg_tl, deng_tl)
  if (present(ntg_tl)) then
    call dm_z_2moment_ad(qmg_tl, ntmg_tl, deng_tl, dmg_tl, zg_tl)
  else
    call dm_z_wsm6_tl('graupel', qmg_tl, deng_tl, dmg_tl, zg_tl)
    endif
  do i = 0, 12
    call coef_a_tl(graupel_coefs(i,0:3), ratg_tl, ag_tl(i))
  end do
  call dualpol_op_icephase_tl(zg_tl, wg_tl, deng_tl, dmg_tl, ag_tl, &
                              zhg_tl, zdrg_tl, kdpg_tl, phvg_tl)
  
  ! 4. pure snow
  call set_traj_precip_type('pure snow')
  call watercontent_tl(qps_tl, wps_tl)
  if (present(nts_tl)) then
    call dm_z_2moment_ad(qps_tl, ntps_tl, denps_tl, dmps_tl, zps_tl)
  else
    call dm_z_wsm6_tl('pure snow', qps_tl, denps_tl, dmps_tl, zps_tl)
  endif
  call dualpol_op_icephase_tl(zps_tl, wps_tl, denps_tl, dmps_tl, aps_tl, &
                              zhps_tl, zdrps_tl, kdpps_tl, phvps_tl)
          
  ! 5. pure graupel
  call set_traj_precip_type('pure graupel')
  call watercontent_tl(qpg_tl, wpg_tl)
  if (present(ntg_tl)) then
    call dm_z_2moment_ad(qpg_tl, ntpg_tl, denpg_tl, dmpg_tl, zpg_tl)
  else
    call dm_z_wsm6_tl('pure graupel', qpg_tl, denpg_tl, dmpg_tl, zpg_tl)
  endif
  call dualpol_op_icephase_tl(zpg_tl, wpg_tl, denpg_tl, dmpg_tl, apg_tl, &
                      zhpg_tl, zdrpg_tl, kdppg_tl, phvpg_tl)
  
  if (present(qh_tl) .and. present(nth_tl)) then
    ! 6. melting graupel
    call set_traj_precip_type('hail')
    call watercontent_tl(qmh_tl, wh_tl)
    call weticephase_density_tl('hail', rath_tl, denh_tl)
    call dm_z_2moment_ad(qmh_tl, ntmh_tl, denh_tl, dmh_tl, zmh_tl)
    do i = 0, 12
      call coef_a_tl(hail_coefs(i,0:3), rath_tl, ah_tl(i))
    end do
    call dualpol_op_icephase_tl(zmh_tl, wh_tl, denh_tl, dmh_tl, ah_tl, &
                              zhh_tl, zdrh_tl, kdph_tl, phvh_tl)

    ! 7. pure snow
    call set_traj_precip_type('pure hail')
    call watercontent_tl(qph_tl, wph_tl)
    call dm_z_2moment_ad(qph_tl, ntph_tl, denph_tl, dmph_tl, zph_tl)
    call dualpol_op_icephase_tl(zph_tl, wph_tl, denph_tl, dmps_tl, aph_tl, &
                              zhph_tl, zdrph_tl, kdpph_tl, phvph_tl)
  endif

  ! total
  if (present(qh_tl) .and. present(nth_tl)) then
    call dualpol_op_total_tl( &
      zhr_tl,  zhs_tl,  zhg_tl, &
      zdrr_tl, zdrs_tl, zdrg_tl, &
      kdpr_tl, kdps_tl, kdpg_tl, &
      phvr_tl, phvs_tl, phvg_tl, &
      zh_tl, zdr_tl, kdp_tl, phv_tl, &
      zhps_tl,  zhpg_tl, &
      zdrps_tl, zdrpg_tl, &
      kdpps_tl, kdppg_tl, &
      phvps_tl, phvpg_tl, &
      zhh_tl=zhh_tl, zdrh_tl=zdrh_tl, &
      kdph_tl=kdph_tl, phvh_tl=phvh_tl, &
      zhph_tl=zhph_tl, zdrph_tl=zdrph_tl, &
      kdpph_tl=kdpph_tl, phvph_tl=phvph_tl )
  else
    call dualpol_op_total_tl( &
      zhr_tl,  zhs_tl,  zhg_tl, &
      zdrr_tl, zdrs_tl, zdrg_tl, &
      kdpr_tl, kdps_tl, kdpg_tl, &
      phvr_tl, phvs_tl, phvg_tl, &
      zh_tl, zdr_tl, kdp_tl, phv_tl, &
      zhps_tl,  zhpg_tl, &
      zdrps_tl, zdrpg_tl, &
      kdpps_tl, kdppg_tl, &
      phvps_tl, phvpg_tl)
  endif

end subroutine ufo_PPRO_sim1obs_tl

subroutine ufo_PPRO_sim1obs_ad(iband, &
  qr_ad, qs_ad, qg_ad, &
  zh_ad, zdr_ad, kdp_ad, phv_ad, &
  qh_ad, ntr_ad, nts_ad, ntg_ad, nth_ad)

  integer, intent(in) :: iband
  real(kind=8), intent(inout) :: qr_ad, qs_ad, qg_ad
  real(kind=8), intent(inout), optional :: qh_ad, ntr_ad, nts_ad, ntg_ad, nth_ad
  real(kind=8), intent(in) :: zh_ad, zdr_ad, kdp_ad, phv_ad

  real (kind=8), dimension(:,:),  pointer :: rain_coefs => NULL()
  real (kind=8), dimension(:,:),  pointer :: snow_coefs => NULL()
  real (kind=8), dimension(:,:),  pointer :: graupel_coefs => NULL()
  real (kind=8), dimension(:,:),  pointer :: hail_coefs => NULL()
  real (kind=8), dimension(:),    pointer :: snow_a => NULL()
  real (kind=8), dimension(:),    pointer :: graupel_a => NULL()
  real (kind=8), dimension(:),    pointer :: hail_a => NULL()
  integer :: i

  ! output of melting scheme
  real (kind=8) :: qpr_ad, qms_ad, qmg_ad, qps_ad, qpg_ad, qmh_ad, qph_ad
  real (kind=8) :: ntpr_ad, ntms_ad, ntmg_ad, ntps_ad, ntpg_ad, ntmh_ad, ntph_ad

  ! "r" is "pr", pure rain
  ! "s" is "ms", melting snow
  ! "g" is "mg", melting graupel
  real (kind=8) :: zhr_ad, zdrr_ad, kdpr_ad, phvr_ad
  real (kind=8) :: zhs_ad, zdrs_ad, kdps_ad, phvs_ad
  real (kind=8) :: zhg_ad, zdrg_ad, kdpg_ad, phvg_ad
  real (kind=8) :: zhh_ad, zdrh_ad, kdph_ad, phvh_ad
  real (kind=8) :: zhps_ad, zdrps_ad, kdpps_ad, phvps_ad
  real (kind=8) :: zhpg_ad, zdrpg_ad, kdppg_ad, phvpg_ad
  real (kind=8) :: zhph_ad, zdrph_ad, kdpph_ad, phvph_ad

  real (kind=8) :: wr_ad, ws_ad, wg_ad, wps_ad, wpg_ad, wh_ad, wph_ad
  real (kind=8) :: zr_ad, zs_ad, zg_ad, zps_ad, zpg_ad, zmh_ad, zph_ad
  real (kind=8) :: dmr_ad, dms_ad, dmg_ad, dmps_ad, dmpg_ad, dmh_ad, dmph_ad
  real (kind=8) :: dens_ad, deng_ad, denh_ad
  real (kind=8) :: rats_ad, ratg_ad, rath_ad
  real (kind=8) :: as_ad(0:12), ag_ad(0:12), ah_ad(0:12)
  real (kind=8) :: aps_ad(0:12), apg_ad(0:12), aph_ad(0:12)
  real (kind=8) :: denr_ad, denps_ad, denpg_ad, denph_ad

  if (iband == 1) then ! s-band radar
      rain_coefs => sband_rain_coefs
      snow_coefs => sband_snow_coefs
      graupel_coefs => sband_graupel_coefs
      hail_coefs => sband_hail_coefs
      snow_a => sband_snow_a
      graupel_a => sband_graupel_a
      hail_a => sband_hail_a
  else ! c-band radar
      rain_coefs => cband_rain_coefs
      snow_coefs => cband_snow_coefs
      graupel_coefs => cband_graupel_coefs
      hail_coefs => cband_hail_coefs
      snow_a => cband_snow_a
      graupel_a => cband_graupel_a
      hail_a => cband_hail_a
  endif

  ! Initialization of temporary adjoint variables
  zhr_ad = 0.0_8; zdrr_ad = 0.0_8; kdpr_ad = 0.0_8; phvr_ad = 0.0_8
  zhs_ad = 0.0_8; zdrs_ad = 0.0_8; kdps_ad = 0.0_8; phvs_ad = 0.0_8
  zhg_ad = 0.0_8; zdrg_ad = 0.0_8; kdpg_ad = 0.0_8; phvg_ad = 0.0_8
  zhh_ad = 0.0_8; zdrh_ad = 0.0_8; kdph_ad = 0.0_8; phvh_ad = 0.0_8
  zhps_ad = 0.0_8; zdrps_ad = 0.0_8; kdpps_ad = 0.0_8; phvps_ad = 0.0_8
  zhpg_ad = 0.0_8; zdrpg_ad = 0.0_8; kdppg_ad = 0.0_8; phvpg_ad = 0.0_8
  zhph_ad = 0.0_8; zdrph_ad = 0.0_8; kdpph_ad = 0.0_8; phvph_ad = 0.0_8

  wr_ad = 0.0_8; ws_ad = 0.0_8; wg_ad = 0.0_8; wps_ad = 0.0_8; wpg_ad = 0.0_8; wh_ad = 0.0_8; wph_ad = 0.0_8
  zr_ad = 0.0_8; zs_ad = 0.0_8; zg_ad = 0.0_8; zps_ad = 0.0_8; zpg_ad = 0.0_8; zmh_ad = 0.0_8; zph_ad = 0.0_8
  dmr_ad = 0.0_8; dms_ad = 0.0_8; dmg_ad = 0.0_8; dmps_ad = 0.0_8; dmpg_ad = 0.0_8; dmh_ad = 0.0_8; dmph_ad = 0.0_8
  denr_ad = 0.0_8; dens_ad = 0.0_8; deng_ad = 0.0_8; denps_ad = 0.0_8; denpg_ad = 0.0_8; denh_ad = 0.0_8; denph_ad = 0.0_8
  rats_ad = 0.0_8; ratg_ad = 0.0_8; rath_ad = 0.0_8
  as_ad(:) = 0.0_8; ag_ad(:) = 0.0_8; aps_ad(:) = 0.0_8; apg_ad(:) = 0.0_8; ah_ad(:) = 0.0_8; aph_ad(:) = 0.0_8

  qpr_ad = 0.0_8; qps_ad = 0.0_8; qpg_ad = 0.0_8; qms_ad = 0.0_8; qmg_ad = 0.0_8; qmh_ad = 0.0_8; qph_ad = 0.0_8
  ntpr_ad = 0.0_8; ntps_ad = 0.0_8; ntpg_ad = 0.0_8; ntms_ad = 0.0_8; ntmg_ad = 0.0_8; ntph_ad = 0.0_8; ntmh_ad = 0.0_8

  ! total
  if ( present(nth_ad) .and. present(qh_ad) ) then
    call dualpol_op_total_ad( &
      zhr_ad,  zhs_ad,  zhg_ad, & 
      zdrr_ad, zdrs_ad, zdrg_ad, & 
      kdpr_ad, kdps_ad, kdpg_ad, & 
      phvr_ad, phvs_ad, phvg_ad, & 
      zh_ad, zdr_ad, kdp_ad, phv_ad, &
      zhps_ad,  zhpg_ad, &
      zdrps_ad, zdrpg_ad, &
      kdpps_ad, kdppg_ad, &
      phvps_ad, phvpg_ad, &
      zhh_ad=zhh_ad, zdrh_ad=zdrh_ad, &
      kdph_ad=kdph_ad, phvh_ad=phvh_ad, &
      zhph_ad=zhph_ad, zdrph_ad=zdrph_ad, &
      kdpph_ad=kdpph_ad, phvph_ad=phvph_ad )
  else
    call dualpol_op_total_ad( &
      zhr_ad,  zhs_ad,  zhg_ad, &
      zdrr_ad, zdrs_ad, zdrg_ad, &
      kdpr_ad, kdps_ad, kdpg_ad, &
      phvr_ad, phvs_ad, phvg_ad, &
      zh_ad, zdr_ad, kdp_ad, phv_ad, &
      zhps_ad,  zhpg_ad, &
      zdrps_ad, zdrpg_ad, &
      kdpps_ad, kdppg_ad, &
      phvps_ad, phvpg_ad)
  endif
  
  ! 1. pure rain
  call set_traj_precip_type('rain')
  call dualpol_op_rain_ad(rain_coefs, wr_ad, dmr_ad, zhr_ad, zdrr_ad, kdpr_ad, phvr_ad)
  if (present(ntr_ad)) then
     call dm_z_2moment_ad(qpr_ad, ntpr_ad, denr_ad, dmr_ad, zr_ad)
  else
     call dm_z_wsm6_ad('rain', qpr_ad, denr_ad, dmr_ad, zr_ad)
  endif
  call watercontent_ad(qpr_ad, wr_ad)

  ! 2. melting snow
  call set_traj_precip_type('snow')
  call dualpol_op_icephase_ad(zs_ad, ws_ad, dens_ad, dms_ad, as_ad, &
      zhs_ad, zdrs_ad, kdps_ad, phvs_ad)
  do i = 0, 12
      call coef_a_ad(snow_coefs(i,0:3), rats_ad, as_ad(i))
  end do
  if (present(nts_ad)) then
    call dm_z_2moment_ad(qms_ad, ntms_ad, dens_ad, dms_ad, zs_ad)
  else
    call dm_z_wsm6_ad('snow', qms_ad, dens_ad, dms_ad, zs_ad)
  endif
  call weticephase_density_ad('snow', rats_ad, dens_ad)
  call watercontent_ad(qms_ad, ws_ad)

  ! 3. melting graupel
  call set_traj_precip_type('graupel')
  call dualpol_op_icephase_ad(zg_ad, wg_ad, deng_ad, dmg_ad, ag_ad, &
      zhg_ad, zdrg_ad, kdpg_ad, phvg_ad)
  do i = 0, 12
      call coef_a_ad(graupel_coefs(i,0:3), ratg_ad, ag_ad(i))
  end do
  if (present(ntg_ad)) then
    call dm_z_2moment_ad(qmg_ad, ntpg_ad, deng_ad, dmg_ad, zg_ad)
  else
    call dm_z_wsm6_ad('graupel', qmg_ad, deng_ad, dmg_ad, zg_ad)
  endif
  call weticephase_density_ad('graupel', ratg_ad, deng_ad)
  call watercontent_ad(qmg_ad, wg_ad)

  ! 4. pure snow
  call set_traj_precip_type('pure snow')
  call dualpol_op_icephase_ad(zps_ad, wps_ad, denps_ad, dmps_ad, aps_ad, &
                              zhps_ad, zdrps_ad, kdpps_ad, phvps_ad)
  if (present(nts_ad)) then
    call dm_z_2moment_ad(qps_ad, ntps_ad, denps_ad, dmps_ad, zps_ad)
  else
    call dm_z_wsm6_ad('pure snow', qps_ad, denps_ad, dmps_ad, zps_ad)
  endif
  call watercontent_ad(qps_ad, wps_ad)
  
  ! 5. pure graupel
  call set_traj_precip_type('pure graupel')
  call dualpol_op_icephase_ad(zpg_ad, wpg_ad, denpg_ad, dmpg_ad, apg_ad, &
                              zhpg_ad, zdrpg_ad, kdppg_ad, phvpg_ad)
  if (present(ntg_ad)) then
    call dm_z_2moment_ad(qpg_ad, ntpg_ad, denpg_ad, dmpg_ad, zpg_ad)
  else
    call dm_z_wsm6_ad('pure graupel', qpg_ad, denpg_ad, dmpg_ad, zpg_ad)
  endif
  call watercontent_ad(qpg_ad, wpg_ad)

  if (present(qh_ad) .and. present(nth_ad)) then
     ! 3. melting hail
     call set_traj_precip_type('hail')
     call dualpol_op_icephase_ad(zmh_ad, wh_ad, denh_ad, dmh_ad, ah_ad, &
        zhh_ad, zdrh_ad, kdph_ad, phvh_ad)
     do i = 0, 12
        call coef_a_ad(hail_coefs(i,0:3), rath_ad, ah_ad(i))
     end do
     call dm_z_2moment_ad(qmh_ad, ntmh_ad, denh_ad, dmh_ad, zmh_ad)
     call weticephase_density_ad('hail', rath_ad, denh_ad)
     call watercontent_ad(qmh_ad, wh_ad)

     ! 4. pure hail
     call set_traj_precip_type('pure hail')
     call dualpol_op_icephase_ad(zph_ad, wph_ad, denph_ad, dmph_ad, aph_ad, &
                              zhph_ad, zdrph_ad, kdpph_ad, phvph_ad)
     call dm_z_2moment_ad(qph_ad, ntph_ad, denph_ad, dmph_ad, zph_ad)
     call watercontent_ad(qph_ad, wph_ad)
  endif

  ! melting scheme
  if (present(qh_ad) .and. present(ntr_ad) .and. present(nts_ad) .and. present(ntg_ad) .and. present(nth_ad)) then
     call melting_scheme_zhang24_ad(qr_ad, qs_ad, qg_ad, qms_ad, qmg_ad, qpr_ad, qps_ad, qpg_ad, rats_ad, ratg_ad, &
          ntr_ad=ntr_ad, nts_ad=nts_ad, ntg_ad=ntg_ad, &
          ntms_ad=ntms_ad, ntmg_ad=ntmg_ad, ntpr_ad=ntpr_ad, ntps_ad=ntps_ad, ntpg_ad=ntps_ad, &
          qh_ad=qh_ad, nth_ad=nth_ad, rath_ad=rath_ad, &
          qmh_ad=qmh_ad, qph_ad=qph_ad, ntmh_ad=ntmh_ad, ntph_ad=ntph_ad )
  else !! WSM6 and Thompson are the same
     call melting_scheme_zhang24_ad(qr_ad, qs_ad, qg_ad, qms_ad, qmg_ad, qpr_ad, qps_ad, qpg_ad, rats_ad, ratg_ad)            
  endif

end subroutine ufo_PPRO_sim1obs_ad

subroutine store_traj(self, iobs)
  implicit none
  type(ufo_PPRO_tlad), intent(inout)  :: self
  integer, intent(in) :: iobs

  self%qr_tray(iobs)=qr_traj; self%qs_tray(iobs)=qs_traj; self%qg_tray(iobs)=qg_traj; self%qh_tray(iobs)=qh_traj
  self%ntr_tray(iobs)=ntr_traj; self%nts_tray(iobs)=nts_traj; self%ntg_tray(iobs)=ntg_traj; self%nth_tray(iobs)=nth_traj

  self%qpr_tray(iobs)=qpr_traj; self%qms_tray(iobs)=qms_traj; self%qmg_tray(iobs)=qmg_traj; self%qmh_tray(iobs)=qmh_traj
  self%qps_tray(iobs)=qps_traj; self%qpg_tray(iobs)=qpg_traj; self%qph_tray(iobs)=qph_traj
  self%ntpr_tray(iobs)=ntpr_traj; self%ntms_tray(iobs)=ntms_traj; self%ntmg_tray(iobs)=ntmg_traj; self%ntmh_tray(iobs)=ntmh_traj
  self%ntps_tray(iobs)=ntps_traj; self%ntpg_tray(iobs)=ntpg_traj; self%ntph_tray(iobs)=ntph_traj

  self%wr_tray(iobs)=wr_traj; self%ws_tray(iobs)=ws_traj; self%wg_tray(iobs)=wg_traj; self%wh_tray(iobs)=wh_traj
  self%wps_tray(iobs)=wps_traj; self%wpg_tray(iobs)=wpg_traj; self%wph_tray(iobs)=wph_traj
  self%dmr_tray(iobs)=dmr_traj; self%dms_tray(iobs)=dms_traj; self%dmg_tray(iobs)=dmg_traj; self%dmh_tray(iobs)=dmh_traj
  self%dmps_tray(iobs)=dmps_traj; self%dmpg_tray(iobs)=dmpg_traj; self%dmph_tray(iobs)=dmph_traj
  self%zr_tray(iobs)=zr_traj; self%zs_tray(iobs)=zs_traj; self%zg_tray(iobs)=zg_traj; self%zmh_tray(iobs)=zh_traj
  self%zps_tray(iobs)=zps_traj; self%zpg_tray(iobs)=zpg_traj; self%zph_tray(iobs)=zph_traj

  self%dens_tray(iobs)=dens_traj; self%deng_tray(iobs)=deng_traj; self%denh_tray(iobs)=denh_traj
  self%rats_tray(iobs)=rats_traj; self%ratg_tray(iobs)=ratg_traj; self%rath_tray(iobs)=rath_traj

  self%zhr_tray(iobs)=zhr_traj; self%zhs_tray(iobs)=zhs_traj; self%zhg_tray(iobs)=zhg_traj; self%zhh_tray(iobs)=zhh_traj
  self%zhps_tray(iobs)=zhps_traj; self%zhpg_tray(iobs)=zhpg_traj; self%zhph_tray(iobs)=zhph_traj
  self%zdrr_tray(iobs)=zdrr_traj; self%zdrs_tray(iobs)=zdrs_traj; self%zdrg_tray(iobs)=zdrg_traj; self%zdrh_tray(iobs)=zdrh_traj
  self%zdrps_tray(iobs)=zdrps_traj; self%zdrpg_tray(iobs)=zdrpg_traj; self%zdrph_tray(iobs)=zdrph_traj
  self%zbarr_tray(iobs)=zbarr_traj; self%zbars_tray(iobs)=zbars_traj; self%zbarg_tray(iobs)=zbarg_traj; self%zbarh_tray(iobs)=zbarh_traj
  self%zbarps_tray(iobs)=zbarps_traj; self%zbarpg_tray(iobs)=zbarpg_traj; self%zbarph_tray(iobs)=zbarph_traj
  self%zvr_tray(iobs)=zvr_traj; self%zvs_tray(iobs)=zvs_traj; self%zvg_tray(iobs)=zvg_traj; self%zvh_tray(iobs)=zvh_traj
  self%zvps_tray(iobs)=zvps_traj; self%zvpg_tray(iobs)=zvpg_traj; self%zvph_tray(iobs)=zvph_traj
  self%kdpr_tray(iobs)=kdpr_traj; self%kdps_tray(iobs)=kdps_traj; self%kdpg_tray(iobs)=kdpg_traj; self%kdph_tray(iobs)=kdph_traj
  self%kdpps_tray(iobs)=kdpps_traj; self%kdppg_tray(iobs)=kdppg_traj; self%kdpph_tray(iobs)=kdpph_traj
  self%phvr_tray(iobs)=phvr_traj; self%phvs_tray(iobs)=phvs_traj; self%phvg_tray(iobs)=phvg_traj; self%phvh_tray(iobs)=phvh_traj
  self%phvps_tray(iobs)=phvps_traj; self%phvpg_tray(iobs)=phvpg_traj; self%phvph_tray(iobs)=phvph_traj

  self%as_tray(iobs,0:12)=as_traj(0:12); self%ag_tray(iobs,0:12)=ag_traj(0:12); self%ah_tray(iobs,0:12)=ah_traj(0:12)

  self%density_air_tray(iobs)=density_air_traj; self%temp_air_tray(iobs)=temp_air_traj; 
  self%zh_tray(iobs)=zh_traj; self%zdr_tray(iobs)=zdr_traj; self%kdp_tray(iobs)=kdp_traj; self%phv_tray(iobs)=phv_traj

end subroutine store_traj

subroutine restore_traj(self, iobs)
  implicit none
  type(ufo_PPRO_tlad), intent(in)  :: self
  integer, intent(in) :: iobs

  qr_traj=self%qr_tray(iobs); qs_traj=self%qs_tray(iobs); qg_traj=self%qg_tray(iobs); qh_traj=self%qh_tray(iobs)
  ntr_traj=self%ntr_tray(iobs); nts_traj=self%nts_tray(iobs); ntg_traj=self%ntg_tray(iobs); nth_traj=self%nth_tray(iobs)

  qpr_traj=self%qpr_tray(iobs); qms_traj=self%qms_tray(iobs); qmg_traj=self%qmg_tray(iobs); qmh_traj=self%qmh_tray(iobs)
  qps_traj=self%qps_tray(iobs); qpg_traj=self%qpg_tray(iobs); qph_traj=self%qph_tray(iobs)
  ntpr_traj=self%ntpr_tray(iobs); ntms_traj=self%ntms_tray(iobs); ntmg_traj=self%ntmg_tray(iobs); ntmh_traj=self%ntmh_tray(iobs)
  ntps_traj=self%ntps_tray(iobs); ntpg_traj=self%ntpg_tray(iobs); ntph_traj=self%ntph_tray(iobs)

  wr_traj=self%wr_tray(iobs); ws_traj=self%ws_tray(iobs); wg_traj=self%wg_tray(iobs); wh_traj=self%wh_tray(iobs)
  wps_traj=self%wps_tray(iobs); wpg_traj=self%wpg_tray(iobs); wph_traj=self%wph_tray(iobs)
  dmr_traj=self%dmr_tray(iobs); dms_traj=self%dms_tray(iobs); dmg_traj=self%dmg_tray(iobs); dmh_traj=self%dmh_tray(iobs)
  dmps_traj=self%dmps_tray(iobs); dmpg_traj=self%dmpg_tray(iobs); dmph_traj=self%dmph_tray(iobs)
  zr_traj=self%zr_tray(iobs); zs_traj=self%zs_tray(iobs); zg_traj=self%zg_tray(iobs); zh_traj=self%zmh_tray(iobs)
  zps_traj=self%zps_tray(iobs); zpg_traj=self%zpg_tray(iobs); zph_traj=self%zph_tray(iobs)

  dens_traj=self%dens_tray(iobs); deng_traj=self%deng_tray(iobs); denh_traj=self%denh_tray(iobs)
  rats_traj=self%rats_tray(iobs); ratg_traj=self%ratg_tray(iobs); rath_traj=self%rath_tray(iobs)

  zhr_traj=self%zhr_tray(iobs); zhs_traj=self%zhs_tray(iobs); zhg_traj=self%zhg_tray(iobs); zhh_traj=self%zhh_tray(iobs)
  zhps_traj=self%zhps_tray(iobs); zhpg_traj=self%zhpg_tray(iobs); zhph_traj=self%zhph_tray(iobs)
  zdrr_traj=self%zdrr_tray(iobs); zdrs_traj=self%zdrs_tray(iobs); zdrg_traj=self%zdrg_tray(iobs); zdrh_traj=self%zdrh_tray(iobs)
  zdrps_traj=self%zdrps_tray(iobs); zdrpg_traj=self%zdrpg_tray(iobs); zdrph_traj=self%zdrph_tray(iobs)
  zbarr_traj=self%zbarr_tray(iobs); zbars_traj=self%zbars_tray(iobs); zbarg_traj=self%zbarg_tray(iobs); zbarh_traj=self%zbarh_tray(iobs)
  zbarps_traj=self%zbarps_tray(iobs); zbarpg_traj=self%zbarpg_tray(iobs); zbarph_traj=self%zbarph_tray(iobs)
  zvr_traj=self%zvr_tray(iobs); zvs_traj=self%zvs_tray(iobs); zvg_traj=self%zvg_tray(iobs); zvh_traj=self%zvh_tray(iobs)
  zvps_traj=self%zvps_tray(iobs); zvpg_traj=self%zvpg_tray(iobs); zvph_traj=self%zvph_tray(iobs)
  kdpr_traj=self%kdpr_tray(iobs); kdps_traj=self%kdps_tray(iobs); kdpg_traj=self%kdpg_tray(iobs); kdph_traj=self%kdph_tray(iobs)
  kdpps_traj=self%kdpps_tray(iobs); kdppg_traj=self%kdppg_tray(iobs); kdpph_traj=self%kdpph_tray(iobs)
  phvr_traj=self%phvr_tray(iobs); phvs_traj=self%phvs_tray(iobs); phvg_traj=self%phvg_tray(iobs); phvh_traj=self%phvh_tray(iobs)
  phvps_traj=self%phvps_tray(iobs); phvpg_traj=self%phvpg_tray(iobs); phvph_traj=self%phvph_tray(iobs)

  as_traj(0:12)=self%as_tray(iobs,0:12); ag_traj(0:12)=self%ag_tray(iobs,0:12); ah_traj(0:12)=self%ah_tray(iobs,0:12)
  ! restore the pointer
  if (self%iband(iobs)==1) then
    aps_traj => sband_snow_a
    apg_traj => sband_graupel_a
    aph_traj => sband_hail_a
  else
    aps_traj => cband_snow_a
    apg_traj => cband_graupel_a
    aph_traj => cband_hail_a
  endif

  density_air_traj=self%density_air_tray(iobs); temp_air_traj=self%temp_air_tray(iobs)
  zh_traj=self%zh_tray(iobs); zdr_traj=self%zdr_tray(iobs); kdp_traj=self%kdp_tray(iobs); phv_traj=self%phv_tray(iobs)

end subroutine restore_traj

subroutine print_traj
  implicit none

  character(len=512) :: buffer

  write(buffer,*) 'qr_traj, qs_traj, qg_traj = ', qr_traj, qs_traj, qg_traj
  call fckit_log%info(buffer); buffer = ''
  write(buffer,*) 'ntr_traj, nts_traj, ntg_traj = ', ntr_traj, nts_traj, ntg_traj
  call fckit_log%info(buffer); buffer = ''

  write(buffer,*) 'qpr_traj, qms_traj, qmg_traj, qps_traj, qpg_traj = ', qpr_traj, qms_traj, qmg_traj, qps_traj, qpg_traj
  call fckit_log%info(buffer); buffer = ''
  write(buffer,*) 'ntpr_traj, ntms_traj, ntmg_traj, ntps_traj, ntpg_traj = ', ntpr_traj, ntms_traj, ntmg_traj, ntps_traj, ntpg_traj
  call fckit_log%info(buffer); buffer = ''

  write(buffer,*) 'wr_traj, ws_traj, wg_traj, wps_traj, wpg_traj = ', wr_traj, ws_traj, wg_traj, wps_traj, wpg_traj
  call fckit_log%info(buffer); buffer = ''
  write(buffer,*) 'dmr_traj, dms_traj, dmg_traj, dmps_traj, dmpg_traj = ', dmr_traj, dms_traj, dmg_traj, dmps_traj, dmpg_traj
  call fckit_log%info(buffer); buffer = ''
  write(buffer,*) 'zr_traj, zs_traj, zg_traj, zps_traj, zpg_traj = ', zr_traj, zs_traj, zg_traj, zps_traj, zpg_traj
  call fckit_log%info(buffer); buffer = ''
  write(buffer,*) 'dens_traj, deng_traj = ', dens_traj, deng_traj
  call fckit_log%info(buffer); buffer = ''
  write(buffer,*) 'rats_traj, ratg_traj = ', rats_traj, ratg_traj
  call fckit_log%info(buffer); buffer = ''
  
  write(buffer,*) 'zhr_traj, zhs_traj, zhg_traj, zhps_traj, zhpg_traj = ', zhr_traj, zhs_traj, zhg_traj, zhps_traj, zhpg_traj
  call fckit_log%info(buffer); buffer = ''
  write(buffer,*) 'zdrr_traj, zdrs_traj, zdrg_traj, zdrps_traj, zdrpg_traj = ', zdrr_traj, zdrs_traj, zdrg_traj, zdrps_traj, zdrpg_traj
  call fckit_log%info(buffer); buffer = ''
  write(buffer,*) 'zbarr_traj, zbars_traj, zbarg_traj, zbarps_traj, zbarpg_traj = ', &
    zbarr_traj, zbars_traj, zbarg_traj, zbarps_traj, zbarpg_traj
  call fckit_log%info(buffer); buffer = ''
  write(buffer,*) 'zvr_traj, zvs_traj, zvg_traj, zvps_traj, zvpg_traj = ', zvr_traj, zvs_traj, zvg_traj, zvps_traj, zvpg_traj
  call fckit_log%info(buffer); buffer = ''
  write(buffer,*) 'kdpr_traj, kdps_traj, kdpg_traj, kdpps_traj, kdppg_traj = ', kdpr_traj, kdps_traj, kdpg_traj, kdpps_traj, kdppg_traj
  call fckit_log%info(buffer); buffer = ''
  write(buffer,*) 'phvr_traj, phvs_traj, phvg_traj, phvps_traj, phvpg_traj = ', phvr_traj, phvs_traj, phvg_traj, phvps_traj, phvpg_traj
  call fckit_log%info(buffer); buffer = ''

  write(buffer,*) 'as_traj = ', as_traj(0:12)
  call fckit_log%info(buffer); buffer = ''
  write(buffer,*) 'ag_traj = ', ag_traj(0:12)
  call fckit_log%info(buffer); buffer = ''
  
  write(buffer,*) 'density_air_traj = ', density_air_traj
  call fckit_log%info(buffer); buffer = ''
  write(buffer,*) 'temp_air_traj = ', temp_air_traj
  call fckit_log%info(buffer); buffer = ''
  write(buffer,*) 'zh_traj, zdr_traj, kdp_traj, phv_traj = ', zh_traj, zdr_traj, kdp_traj, phv_traj
  call fckit_log%info(buffer); buffer = ''

end subroutine print_traj

subroutine allocate_traj_trays(self)
  implicit none
  type(ufo_PPRO_tlad), intent(inout)  :: self

  integer :: nobs 
  
  nobs = self%nlocs

  allocate(self%qr_tray(nobs), self%qs_tray(nobs), self%qg_tray(nobs), self%qh_tray(nobs))
  allocate(self%ntr_tray(nobs), self%nts_tray(nobs), self%ntg_tray(nobs),  self%nth_tray(nobs))

  allocate(self%qpr_tray(nobs), self%qms_tray(nobs), self%qmg_tray(nobs), self%qps_tray(nobs), self%qpg_tray(nobs), self%qph_tray(nobs)) 
  allocate(self%ntpr_tray(nobs), self%ntms_tray(nobs), self%ntmg_tray(nobs), self%ntps_tray(nobs), self%ntpg_tray(nobs), self%ntmh_tray(nobs), self%ntph_tray(nobs))

  allocate(self%wr_tray(nobs), self%ws_tray(nobs), self%wg_tray(nobs), self%wps_tray(nobs), self%wpg_tray(nobs), self%wh_tray(nobs), self%wph_tray(nobs))
  allocate(self%dmr_tray(nobs), self%dms_tray(nobs), self%dmg_tray(nobs), self%dmps_tray(nobs), self%dmpg_tray(nobs), self%dmh_tray(nobs), self%dmph_tray(nobs)) 
  allocate(self%zr_tray(nobs), self%zs_tray(nobs), self%zg_tray(nobs), self%zps_tray(nobs), self%zpg_tray(nobs), self%zmh_tray(nobs), self%zph_tray(nobs))
  allocate(self%dens_tray(nobs), self%deng_tray(nobs), self%denh_tray(nobs))
  allocate(self%rats_tray(nobs), self%ratg_tray(nobs), self%rath_tray(nobs))

  allocate(self%zhr_tray(nobs), self%zhs_tray(nobs), self%zhg_tray(nobs), self%zhps_tray(nobs), self%zhpg_tray(nobs), self%zhh_tray(nobs), self%zhph_tray(nobs))
  allocate(self%zdrr_tray(nobs), self%zdrs_tray(nobs), self%zdrg_tray(nobs), self%zdrps_tray(nobs), self%zdrpg_tray(nobs), self%zdrh_tray(nobs), self%zdrph_tray(nobs))
  allocate(self%zbarr_tray(nobs), self%zbars_tray(nobs), self%zbarg_tray(nobs), self%zbarps_tray(nobs), self%zbarpg_tray(nobs), self%zbarh_tray(nobs), self%zbarph_tray(nobs))
  allocate(self%zvr_tray(nobs), self%zvs_tray(nobs), self%zvg_tray(nobs), self%zvps_tray(nobs), self%zvpg_tray(nobs), self%zvh_tray(nobs), self%zvph_tray(nobs))
  allocate(self%kdpr_tray(nobs), self%kdps_tray(nobs), self%kdpg_tray(nobs), self%kdpps_tray(nobs), self%kdppg_tray(nobs), self%kdph_tray(nobs), self%kdpph_tray(nobs))
  allocate(self%phvr_tray(nobs), self%phvs_tray(nobs), self%phvg_tray(nobs), self%phvps_tray(nobs), self%phvpg_tray(nobs), self%phvh_tray(nobs), self%phvph_tray(nobs))

  allocate(self%as_tray(nobs,0:12), self%ag_tray(nobs,0:12),  self%ah_tray(nobs,0:12))

  allocate(self%density_air_tray(nobs), self%temp_air_tray(nobs))
  allocate(self%zh_tray(nobs), self%zdr_tray(nobs), self%kdp_tray(nobs), self%phv_tray(nobs))
  allocate(self%badbg_tray(nobs))
  self%badbg_tray(:) = .False.

end subroutine allocate_traj_trays

subroutine deallocate_traj_trays(self)
  implicit none
  type(ufo_PPRO_tlad), intent(inout)  :: self

  if (allocated(self%qr_tray)) then
    deallocate(self%qr_tray, self%qs_tray, self%qg_tray, self%qh_tray)
    deallocate(self%ntr_tray, self%nts_tray, self%ntg_tray, self%nth_tray)

    deallocate(self%qpr_tray, self%qms_tray, self%qmg_tray, self%qps_tray, self%qpg_tray, self%qmh_tray, self%qph_tray) 
    deallocate(self%ntpr_tray, self%ntms_tray, self%ntmg_tray, self%ntps_tray, self%ntpg_tray, self%ntmh_tray, self%ntph_tray)

    deallocate(self%wr_tray, self%ws_tray, self%wg_tray, self%wps_tray, self%wpg_tray, self%wh_tray, self%wph_tray)
    deallocate(self%dmr_tray, self%dms_tray, self%dmg_tray, self%dmps_tray, self%dmpg_tray, self%dmh_tray, self%dmph_tray) 
    deallocate(self%zr_tray, self%zs_tray, self%zg_tray, self%zps_tray, self%zpg_tray, self%zmh_tray, self%zph_tray)
    deallocate(self%dens_tray, self%deng_tray, self%denh_tray)
    deallocate(self%rats_tray, self%ratg_tray, self%rath_tray)

    deallocate(self%zhr_tray, self%zhs_tray, self%zhg_tray, self%zhps_tray, self%zhpg_tray, self%zhh_tray, self%zhph_tray)
    deallocate(self%zdrr_tray, self%zdrs_tray, self%zdrg_tray, self%zdrps_tray, self%zdrpg_tray, self%zdrh_tray, self%zdrph_tray)
    deallocate(self%zbarr_tray, self%zbars_tray, self%zbarg_tray, self%zbarps_tray, self%zbarpg_tray, self%zbarh_tray, self%zbarph_tray)
    deallocate(self%zvr_tray, self%zvs_tray, self%zvg_tray, self%zvps_tray, self%zvpg_tray, self%zvh_tray, self%zvph_tray)
    deallocate(self%kdpr_tray, self%kdps_tray, self%kdpg_tray, self%kdpps_tray, self%kdppg_tray, self%kdph_tray, self%kdpph_tray)
    deallocate(self%phvr_tray, self%phvs_tray, self%phvg_tray, self%phvps_tray, self%phvpg_tray, self%phvh_tray, self%phvph_tray)

    deallocate(self%as_tray, self%ag_tray, self%ah_tray)

    deallocate(self%density_air_tray, self%temp_air_tray)
    deallocate(self%zh_tray, self%zdr_tray, self%kdp_tray, self%phv_tray)
    deallocate(self%badbg_tray)
  endif

end subroutine deallocate_traj_trays

! ------------------------------------------------------------------------------

end module ufo_PPRO_tlad_mod
