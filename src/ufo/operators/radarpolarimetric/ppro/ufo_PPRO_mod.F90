! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for radarradialvelocity observation operator

module ufo_PPRO_mod
!> @brief UFO interface module for PPRO dual-polarization radar operator
!>
!> This module provides the UFO interface layer for the PPRO operator.
!> Core physics are in the external ppro-lib library.
!> 
!> Authors:
!>   Zhiquan (Jake) Liu (NCAR/MMM) - Implementation of the Zhang21 operator
!>   Tzu-Chin Tsai (CWA,Taiwan) - Implementation of the TCWA2 operator
!>   Hejun Xie  - Integration into the JEDI-UFO framework;
!>                development of TL/AD components
!>   Tao Sun - Extension to include hail categories and advanced microphysics support
!>   Rong Kong (NCAR/MMM) - Bug fixes, operator tuning and testing,
!>                          development of a multi-operator architecture, and
!>                          modularization with an external ppro library



 use oops_variables_mod
 use obs_variables_mod
 use ufo_vars_mod
 use dualpol_op_mod, only: ppro_init_coefs, ppro_compute_point, ppro_set_operator
 use missing_values_mod
 use fckit_log_module, only: fckit_log

 implicit none

 integer, parameter, public :: nopvar = 4 ! count available operator variables: zh, zdr, kdp, phv 
 character(len=maxvarlen), dimension(nopvar), public :: opvars_list = (/ &
      'equivalentReflectivityFactor ', & ! zh
      'differentialReflectivity     ', & ! zdr  
      'specificDifferentialPhase    ', & ! kdp
      'copolarCorrelationCoefficient'/)  ! phv

!> Fortran derived type for the observation type
 type, public :: ufo_PPRO
 private
   type(obs_variables), public :: obsvars
   type(oops_variables), public :: geovars
   integer, public :: obsvarindices(nopvar) ! Array that maps indices of opvars to indices of hofx columns
                                       ! -1 means not used
   character(len=MAXVARLEN), public :: v_coord ! GeoVaL to use to interpolate in vertical
   character(len=MAXVARLEN), public :: micro_option    ! Choice (enum) of microphysics option
 contains
   procedure :: setup  => ufo_PPRO_setup
   procedure :: simobs => ufo_PPRO_simobs
 end type ufo_PPRO

 private
 
 ! Public subroutines
 public :: ufo_PPRO_sim1obs  ! For TL/AD trajectory setting

 integer :: n_geovars
 character(len=maxvarlen), dimension(:), allocatable :: geovars_list

contains

! ------------------------------------------------------------------------------
subroutine ufo_PPRO_setup(self, yaml_conf)
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
implicit none
class(ufo_PPRO), intent(inout)     :: self
type(fckit_configuration), intent(in) :: yaml_conf

character(kind=c_char,len=:), allocatable :: coord_name
character(kind=c_char,len=:), allocatable :: micro_option
character(kind=c_char,len=:), allocatable :: ppro_operator
character(kind=c_char,len=:), allocatable :: this_varname
character(len=maxvarlen) :: var_string
character(len=512) :: buffer
integer :: i, j
logical :: found

  ! Initialize all as not used variables
  self%obsvarindices(:) = -1

  DO i = 1, self%obsvars%nvars()
    write(buffer,*) 'PPRO obsvars(', i, ') = ', TRIM(self%obsvars%variable(i))
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

  write(buffer,*) 'PPRO obsvarindices = ', self%obsvarindices(:)
  call fckit_log%info(buffer); buffer = ''

  ! YAML option for polarimetric operator (default: Zhang21)
  if( yaml_conf%has("polarimetric operator") ) then
     call yaml_conf%get_or_die("polarimetric operator", ppro_operator)
     call ppro_set_operator(trim(ppro_operator))
  else
     call ppro_set_operator('Zhang21')  ! Default operator
  endif
  
  ! YAML option for microphysics scheme
  call yaml_conf%get_or_die("microphysics option", micro_option)
  if (trim(micro_option) .eq. "Thompson") then
    self%micro_option = 'THOMPSON'
    n_geovars=8
  else if (trim(micro_option) .eq. "WSM6") then
    self%micro_option = 'WSM6'
    n_geovars=7
  else if (trim(micro_option) .eq. "NSSL") then
    self%micro_option = 'NSSL'
    n_geovars=14
  else if (trim(micro_option) .eq. "TCWA2") then
    self%micro_option = 'TCWA2'
    n_geovars=15
  else
    print*, ' microphysics picked is: ', trim(micro_option)
    call abor1_ftn("microphysics option not set or unsupported, aborting")
  endif

  if ( .not. allocated(geovars_list) ) allocate(geovars_list(n_geovars))
  geovars_list(1) = var_airdens
  geovars_list(2) = var_ts
  geovars_list(3) = var_prs
  geovars_list(4) = var_q
  !! WSM6 scheme
  ! rain water mixing ratio, kg/kg
  var_string="var_rain_mixing_ratio"
  if( yaml_conf%has(trim(var_string)) ) then
     call yaml_conf%get_or_die(trim(var_string), this_varname)
     geovars_list(5) = this_varname
  endif
  ! snow mixing ratio, kg/kg
  var_string="var_snow_mixing_ratio"
  if( yaml_conf%has(trim(var_string)) ) then
     call yaml_conf%get_or_die(trim(var_string), this_varname)
     geovars_list(6) = this_varname
  endif
  ! graupel mxing ratio, kg/kg
  var_string="var_graupel_mixing_ratio"
  if( yaml_conf%has(trim(var_string)) ) then
     call yaml_conf%get_or_die(trim(var_string), this_varname)
     geovars_list(7) = this_varname
  endif
  !! Thompson shceme
  if ( trim(self%micro_option) .eq. "THOMPSON" ) then
     ! number concentration of rain water, #/kg
     var_string="var_rain_number_concentration"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(8) = this_varname
     endif
  endif
  !! NSSL scheme
  if ( trim(self%micro_option) .eq. "NSSL" ) then
     ! hail mxing ratio, kg/kg
     var_string="var_hail_mixing_ratio"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(8) = this_varname
     endif
     ! number concentration of rain water, #/kg
     var_string="var_rain_number_concentration"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(9) = this_varname
     endif
     ! number concentration of snow, #/kg
     var_string="var_snow_number_concentration"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(10) = this_varname
     endif
     ! number concentration of graupel, #/kg
     var_string="var_graupel_number_concentration"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(11) = this_varname
     endif
     ! number concentration of hail, #/kg
     var_string="var_hail_number_concentration"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(12) = this_varname
     endif

     ! volume mixing ratio of graupel, m^3/kg
     var_string="var_graupel_vol_mixing_ratio"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(13) = this_varname
     endif
    ! volume mixing ratio of hail, m^3/kg
     var_string="var_hail_vol_mixing_ratio"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(14) = this_varname
     endif

  endif

  !! TCWA2 scheme
  if ( trim(self%micro_option) .eq. "TCWA2" ) then
     ! number concentration of rain water, #/kg
     var_string="var_rain_number_concentration"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(8) = this_varname
     endif
     ! number concentration of snow, #/kg
     var_string="var_snow_number_concentration"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(9) = this_varname
     endif
     ! number concentration of graupel, #/kg
     var_string="var_graupel_number_concentration"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(10) = this_varname
     endif
     ! mixing ratio of cloud, kg/kg
     var_string="var_cloud_mixing_ratio"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(11) = this_varname
     endif
     ! mixing ratio of ice, kg/kg
     var_string="var_ice_mixing_ratio"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(12) = this_varname
     endif
     ! number concentration of ice, #/kg
     var_string="var_ice_number_concentration"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(13) = this_varname
     endif
     ! melted fraction of snow
     var_string="var_snow_melted_fraction"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(14) = this_varname
     endif
     ! melted fraction of graupel
     var_string="var_graupel_melted_fraction"
     if( yaml_conf%has(trim(var_string)) ) then
        call yaml_conf%get_or_die(trim(var_string), this_varname)
        geovars_list(15) = this_varname
     endif
  endif

  ! YAML option for vertical coordinate name
  call yaml_conf%get_or_die("VertCoord",coord_name)
  self%v_coord = coord_name
  if( trim(self%v_coord) .ne. var_z ) then
      write(*,'(A)') 'ERROR: Unsupported vertical coordinate: ' // trim(self%v_coord)
      write(*,'(A)') 'Supported coordinate: height (var_z)'
      call abor1_ftn("ufo_PPRO: incorrect vertical coordinate specified")
  endif

  call self%geovars%push_back(geovars_list)
  call self%geovars%push_back(self%v_coord)

end subroutine ufo_PPRO_setup

! ------------------------------------------------------------------------------
! Code in this routine is for radar reflectivity only
subroutine ufo_PPRO_simobs(self, geovals, obss, nvars, nlocs, hofx)
  use kinds
  use vert_interp_mod
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use obsspace_mod
  implicit none
  class(ufo_PPRO), intent(in)    :: self
  integer, intent(in)               :: nvars, nlocs
  type(ufo_geovals),  intent(in)    :: geovals
  real(c_double),     intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value, intent(in)    :: obss

  ! Local variables
  logical :: variable_present
  integer :: iobs, ivar, nvars_geovars
  real(kind_real),  dimension(:), allocatable :: obsvcoord
  real(kind_real),  dimension(:), allocatable :: zhobs, zdrobs, kdpobs
  type(ufo_geoval), pointer :: vcoordprofile, profile
  real(kind_real),  allocatable :: wf(:)
  integer,          allocatable :: wi(:)
  integer,          allocatable :: iband(:) ! s-band=1 or c-band=2

  character(len=MAXVARLEN) :: geovar

  real(kind_real), allocatable :: tmp(:)
  real(kind_real) :: tmp2

  real(kind_real), allocatable :: fields(:,:)  ! background fields interplated vertically to obs height

  integer       :: i
  integer       :: obsvarindex
  real (kind=8) :: rhow,rho ! kg/m^3
  real (kind=8) :: t ! temperature, Kelvin
  real (kind=8) :: qr, qs, qg, qh ! mixing ratio, g/kg
  real (kind=8) :: qc, qi             ! cloud water and ice mixing ratio, g/kg (TCWA2)
  real (kind=8) :: nr, ns, ng, nh ! number concentration
  real (kind=8) :: ni                 ! ice number concentration (TCWA2)
  real (kind=8) :: vg, vh         ! volume mixing ratio
  real (kind=8) :: smlf, gmlf         ! melted fractions (TCWA2)
  real (kind=8) :: zh, zdr, kdp, phv
  real (kind=8) :: zh_dBZ, zdr_dB
  real(c_double) :: missing

  real (kind=8) :: p, q
  real(kind_real), parameter :: rd=287.04_kind_real
  real(kind_real), parameter :: one=1._kind_real
  real(kind_real), parameter :: D608=0.608_kind_real

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Get height profiles from geovals
  call ufo_geovals_get_var(geovals, self%v_coord, vcoordprofile)

! Get the observation vertical coordinates
  allocate(obsvcoord(nlocs))
  allocate(iband(nlocs))
  allocate(zhobs(nlocs))
  allocate(zdrobs(nlocs))
  allocate(kdpobs(nlocs))
  !call obsspace_get_db(obss, "MetaData", "geometric_height", obsvcoord)
  call obsspace_get_db(obss, "MetaData", "height", obsvcoord)
 
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  call obsspace_get_db(obss, "ObsValue", "equivalentReflectivityFactor", zhobs)
  call obsspace_get_db(obss, "ObsValue", "differentialReflectivity", zdrobs)
  call obsspace_get_db(obss, "ObsValue", "specificDifferentialPhase", kdpobs)
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  variable_present = obsspace_has(obss, "MetaData", "iband")
  if (variable_present) then
     call obsspace_get_db(obss, "MetaData", "iband", iband)
  else
     iband(:) = 1
  endif

! put observation operator code here
! Allocate arrays for interpolation weights

  allocate(wi(nlocs))
  allocate(wf(nlocs))

! Calculate the interpolation weights

  allocate(tmp(vcoordprofile%nval)) ! one profile with 'nval' levels
  do iobs = 1, nlocs
    tmp = vcoordprofile%vals(:,iobs) ! model vertical coordinate profile at iobs location
    tmp2 = obsvcoord(iobs)           ! obs vertical coordinate (single value) at iobs location
    call vert_interp_weights(vcoordprofile%nval, tmp2, tmp, wi(iobs), wf(iobs))
  enddo

! Number of variables in geovars (without the vertical coordinate)
  nvars_geovars = self%geovars%nvars() - 1
  allocate(fields(nvars_geovars,nlocs))

  do ivar = 1, nvars_geovars
! Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)

! Get profiles for this variable from geovals, profiles at nlocs
    call ufo_geovals_get_var(geovals, geovar, profile)
    
! Interpolate from geovals to observational location into hofx
    do iobs = 1, nlocs
      call vert_interp_apply(profile%nval, profile%vals(:,iobs), &
                             & fields(ivar,iobs), wi(iobs), wf(iobs))
    enddo
  enddo

! Apply parameterized dual-pol radar operator
  missing = missing_value(missing)

  ! Initialize coefficients only once (IMPORTANT for performance!)
  call ppro_init_coefs()

  do iobs=1,nlocs
    !if(iband(iobs) == 2) cycle   !Rong Kong temporarily added for testing
    rhow  = fields(1,iobs) ! moist air density, kg/m^3
    t  = fields(2,iobs) ! temperature (K)
    p  = fields(3,iobs) ! pressure (Pa) 
    q  = fields(4,iobs) ! specific humidity
    q = q/(one-q) !convert from specific humidity to mixing ratio
    rho = p/(rd*t*(one+D608*q))  !calculate dry air density, kg/m^3
    !!! print*,'t=',t,'p=',p,'q=','rho = ',rho, " rhow=", rhow
    qr   = 1000.0*fields(5,iobs) ! kg/kg -> g/kg
    qs   = 1000.0*fields(6,iobs) ! kg/kg -> g/kg
    qg   = 1000.0*fields(7,iobs) ! kg/kg -> g/kg

    if ( trim(self%micro_option) .eq. "THOMPSON" ) then
       nr = fields(8,iobs)*rho
    else if ( trim(self%micro_option) .eq. "NSSL" ) then
       qh = 1000.0*fields(8,iobs)
       nr = fields(9,iobs)*rho  ! #/kg   x  kg/m^3  =  #/m^3
       ns = fields(10,iobs)*rho ! #/kg   x  kg/m^3  =  #/m^3
       ng = fields(11,iobs)*rho ! #/kg   x  kg/m^3  =  #/m^3
       nh = fields(12,iobs)*rho ! #/kg   x  kg/m^3  =  #/m^3
       vg = fields(13,iobs)*rho ! m^3/kg x  kg/m^3  =  m^3/m^3
       vh = fields(14,iobs)*rho ! m^3/kg x  kg/m^3  =  m^3/m^3
    else if ( trim(self%micro_option) .eq. "TCWA2" ) then
       nr = fields(8,iobs)*rho     ! #/kg x kg/m^3 = #/m^3
       ns = fields(9,iobs)*rho     ! #/kg x kg/m^3 = #/m^3
       ng = fields(10,iobs)*rho    ! #/kg x kg/m^3 = #/m^3
       qc = 1000.0*fields(11,iobs) ! kg/kg -> g/kg
       qi = 1000.0*fields(12,iobs) ! kg/kg -> g/kg
       ni = fields(13,iobs)*rho    ! #/kg x kg/m^3 = #/m^3
       smlf = fields(14,iobs)      ! melted fraction of snow
       gmlf = fields(15,iobs)      ! melted fraction of graupel
    endif

    if (rho .LT. 0.0) then
      hofx(:,iobs) = missing
      cycle
    endif

    ! Call ppro-lib core computation
    if ( trim(self%micro_option) .eq. "WSM6" ) then
       call ppro_compute_point(iband(iobs), self%micro_option, rho, t, &
                               qr, qs, qg, zh, zdr, kdp, phv)

    else if ( trim(self%micro_option) .eq. "THOMPSON" ) then
       call ppro_compute_point(iband(iobs), self%micro_option, rho, t, &
                               qr, qs, qg, zh, zdr, kdp, phv, &
                          nr=nr)

    else if ( trim(self%micro_option) .eq. "NSSL" ) then
       call ppro_compute_point(iband(iobs), self%micro_option, rho, t, &
                               qr, qs, qg, zh, zdr, kdp, phv, &
                               qh=qh, nr=nr, ns=ns, ng=ng, nh=nh, vg=vg, vh=vh)

    else if ( trim(self%micro_option) .eq. "TCWA2" ) then
       call ppro_compute_point(iband(iobs), self%micro_option, rho, t, &
                               qr, qs, qg, zh, zdr, kdp, phv, &
                               nr=nr, ns=ns, ng=ng, qi=qi, ni=ni, qc=qc, smlf=smlf, gmlf=gmlf)
    end if

    if (zh < 1.0_kind_real)then
       zh = zh + 1.0_kind_real  !K.R.
    endif 

    !if (zdr < 1.0_kind_real)then
    !   zdr = zdr + 1.0E-2_kind_real  !K.R.???
    !endif

    if (zh > 0.0) then
      zh_dBZ = 10.0*log10(zh)
    else
      zh_dBZ = missing
    endif

    if (zdr > 0.0) then
      zdr_dB = 10.0*log10(zdr)
    else
      zdr_dB = missing
    endif

    do ivar=1, nopvar
      obsvarindex = self%obsvarindices(ivar)
      if (obsvarindex .EQ. -1) CYCLE ! VARIABLE NOT USED
      select case(ivar)
        case (1) ! zh
          hofx(obsvarindex,iobs) = zh_dBZ
        case (2) ! zdr
          hofx(obsvarindex,iobs) = zdr_dB
        case (3) ! kdp
          hofx(obsvarindex,iobs) = kdp
        case (4) ! phv
          hofx(obsvarindex,iobs) = phv
      end select
    enddo

  enddo

! Cleanup memory
  deallocate(obsvcoord)
  deallocate(zhobs, zdrobs, kdpobs)
  deallocate(wi)
  deallocate(wf)

  deallocate(fields)

end subroutine ufo_PPRO_simobs

! ==============================================================================
!> @brief Wrapper for TL/AD trajectory setting
!>
!> This function is kept for backward compatibility with the TL/AD module.
!> It calls ppro_compute_point and then sets the trajectory.
! ==============================================================================
subroutine ufo_PPRO_sim1obs(iband, density_air, temp_air, qr, qs, qg, zh, zdr, kdp, phv, qh, & 
                           nr, ns, ng, nh, vg, vh, set_trajectory, height, zhobs, zdrobs, kdpobs)
  use zhang21_tlad_mod, only: qr_traj, qs_traj, qg_traj, &
                              zh_traj, zdr_traj, kdp_traj, phv_traj, &
                              density_air_traj, temp_air_traj
  implicit none
  
  integer, intent(in) :: iband
  real(kind=8), intent(in) :: density_air, temp_air
  real(kind=8), intent(in) :: qr, qs, qg
  real(kind=8), intent(out) :: zh, zdr, kdp, phv
  real(kind=8), intent(in), optional :: qh, nr, ns, ng, nh, vg, vh
  real(kind=8), intent(in), optional :: height, zhobs, zdrobs, kdpobs 
  logical, optional, intent(in) :: set_trajectory
  
  ! Determine scheme type (for backward compatibility)
  character(len=20) :: scheme_type
  
  if (present(nr) .and. present(ns) .and. present(ng) .and. present(nh) .and. present(qh)) then
    scheme_type = 'NSSL'
  elseif (present(nr) .and. present(ns) .and. present(ng)) then
    scheme_type = 'TCWA2'
  elseif (present(nr)) then
    scheme_type = 'THOMPSON'
  else
    scheme_type = 'WSM6'
  endif
  
  ! Call ppro-lib core function
  call ppro_compute_point(iband, scheme_type, density_air, temp_air, &
                          qr, qs, qg, zh, zdr, kdp, phv, &
                          qh, nr, ns, ng, nh, vg, vh)
  
  ! Set trajectory if requested (for TL/AD)
  if (present(set_trajectory)) then
    if (set_trajectory) then
      ! Note: The trajectory variables are set inside ppro-lib functions
      ! Here we just record the final state
      qr_traj = qr
      qs_traj = qs
      qg_traj = qg
      zh_traj = zh
      zdr_traj = zdr
      kdp_traj = kdp
      phv_traj = phv
      density_air_traj = density_air
      temp_air_traj = temp_air
    endif
  endif
  
end subroutine ufo_PPRO_sim1obs

! ------------------------------------------------------------------------------

end module ufo_PPRO_mod
