module ocn_yakl

   use iso_c_binding
   use mpas_kind_types
   use ocn_constants

   implicit none
   private

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: &
        ocn_init_yakl, &
        ocn_finalize_yakl, & 
        ocn_yakl_fence, &
        ocn_yakl_update_host, &
        ocn_yakl_update_device, &
        ocn_yakl_create, &
        ocn_yakl_delete, &
        ocn_yakl_wrap, &
        ocn_yakl_unwrap, &
        ocn_yakl_type

   type, bind(C) :: ocn_yakl_type
    integer(c_int), dimension(10) :: ashape
    type(c_ptr) :: aloc
   end type

!   type ocn_yakl_type
!      type(c_ptr) :: wptr
!   end type


   interface ocn_yakl_update_host
      module procedure ocn_yakl_update_host_d1d
      module procedure ocn_yakl_update_host_d2d
   end interface

   interface ocn_yakl_update_device
      module procedure ocn_yakl_update_device_d1d
      module procedure ocn_yakl_update_device_d2d
   end interface

   interface ocn_yakl_create
      module procedure ocn_yakl_create_d1d
      module procedure ocn_yakl_create_d2d
   end interface

   interface ocn_yakl_wrap
      module procedure ocn_yakl_wrap_array_d1d
      module procedure ocn_yakl_wrap_array_d2d
      !module procedure ocn_yakl_wray_array_i1d
      !module procedure ocn_yakl_wray_array_i2d
   end interface

   interface ocn_yakl_unwrap
      module procedure ocn_yakl_unwrap_array_d1d
      module procedure ocn_yakl_unwrap_array_d2d
      !module procedure ocn_yakl_wray_array_i1d
      !module procedure ocn_yakl_wray_array_i2d
   end interface

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_yakl
!
!> \brief  initialize yakl library
!> \author Rob Aulwes
!> \date   October 2021
!> \details
!>  This routine initializes the YAKL library
!
!-----------------------------------------------------------------------
   subroutine ocn_init_yakl!{{{
   
      interface
        subroutine init_yakl() bind(C,name="init_yakl")
          use iso_c_binding
        end subroutine init_yakl
      end interface
   
      call init_yakl
     
   end subroutine ocn_init_yakl

!***********************************************************************
!
!  routine ocn_init_yakl
!
!> \brief  initialize yakl library
!> \author Rob Aulwes
!> \date   October 2021
!> \details
!>  This routine initializes the YAKL library
!
!-----------------------------------------------------------------------
   subroutine ocn_finalize_yakl!{{{
   
      interface
        subroutine finalize_yakl() bind(C,name="finalize_yakl")
          use iso_c_binding
        end subroutine finalize_yakl
      end interface
   
      call finalize_yakl
     
   end subroutine ocn_finalize_yakl

!***********************************************************************
!
!  routine ocn_yakl_fence
!
!> \brief  barrier to wait for all YAKL actions to complete
!> \author Rob Aulwes
!> \date   October 2021
!> \details
!>  This routine initializes the YAKL library
!
!-----------------------------------------------------------------------
   subroutine ocn_yakl_fence!{{{
   
      interface
        subroutine yakl_fence() bind(C,name="yakl_fence")
          use iso_c_binding
        end subroutine yakl_fence
      end interface
   
      call yakl_fence
     
   end subroutine ocn_yakl_fence

!***********************************************************************
!
!  routine ocn_yakl_create_d1d
!
!> \brief  initialize yakl library
!> \author Rob Aulwes
!> \date   October 2021
!> \details
!>  This routine initializes the YAKL library
!
!-----------------------------------------------------------------------
   subroutine ocn_yakl_create_d1d(yarr, nn)!{{{
   
      type(c_ptr), intent(out) :: yarr
      integer, intent(in), value :: nn

     interface
       function yakl_create_d1d(nn) bind(C,name="yakl_create_d1d")
          use iso_c_binding
          integer, intent(in), value :: nn
          type(c_ptr) :: yakl_create_d1d
        end function yakl_create_d1d
      end interface
   
      yarr = yakl_create_d1d(nn)
     
   end subroutine ocn_yakl_create_d1d

!***********************************************************************
!
!  routine ocn_yakl_create_d2d
!
!> \brief  initialize yakl library
!> \author Rob Aulwes
!> \date   October 2021
!> \details
!>  This routine initializes the YAKL library
!
!-----------------------------------------------------------------------
   subroutine ocn_yakl_create_d2d(yarr, nn, mm)!{{{
   
      type(c_ptr), intent(out) :: yarr
      integer, intent(in), value :: nn, mm

     interface
       function yakl_create_d2d(nn, mm) bind(C,name="yakl_create_d2d")
          use iso_c_binding
          integer, intent(in), value :: nn, mm
          type(c_ptr) :: yakl_create_d2d
        end function yakl_create_d2d
      end interface
   
      yarr = yakl_create_d2d(nn, mm)
     
   end subroutine ocn_yakl_create_d2d

!***********************************************************************
!
!  routine ocn_yakl_create_d1d
!
!> \brief  initialize yakl library
!> \author Rob Aulwes
!> \date   October 2021
!> \details
!>  This routine initializes the YAKL library
!
!-----------------------------------------------------------------------
   subroutine ocn_yakl_delete(yarr)!{{{
   
      type(c_ptr), intent(in) :: yarr

     interface
       subroutine yakl_delete(yarr) bind(C,name="yakl_delete")
          use iso_c_binding
      	  type(c_ptr), intent(in), value :: yarr	
        end subroutine yakl_delete
      end interface
   
      call yakl_delete(yarr)
     
   end subroutine ocn_yakl_delete

!***********************************************************************
!
!  routine ocn_yakl_wrap_array_d2d
!
!> \brief  initialize yakl library
!> \author Rob Aulwes
!> \date   October 2021
!> \details
!>  This routine initializes the YAKL library
!
!-----------------------------------------------------------------------
   subroutine ocn_yakl_wrap_array_d1d(yarr, arr)!{{{
   
      real (kind=RKIND), dimension(:), intent(in) :: arr
      type(c_ptr), intent(out) :: yarr
      
      integer :: n,m

     interface
       function wrap_double_array_1d(var1, nn) bind(C,name="wrap_array_1d")
          use iso_c_binding
          real(8), dimension(*), intent(in) :: var1
          integer, intent(in), value :: nn
          type(c_ptr) :: wrap_double_array_1d
        end function wrap_double_array_1d
      end interface
   
      n = size(arr,1)
      yarr = wrap_double_array_1d(arr, n)
     
   end subroutine ocn_yakl_wrap_array_d1d

!***********************************************************************
!
!  routine ocn_yakl_wrap_array_d2d
!
!> \brief  initialize yakl library
!> \author Rob Aulwes
!> \date   October 2021
!> \details
!>  This routine initializes the YAKL library
!
!-----------------------------------------------------------------------
   subroutine ocn_yakl_wrap_array_d2d(yarr, arr)!{{{
   
      real (kind=RKIND), dimension(:,:), intent(in) :: arr
      type(c_ptr), intent(out) :: yarr
      
      integer :: n,m

     interface
       function wrap_double_array_2d(var1, nn, mm) bind(C,name="wrap_array_2d")
          use iso_c_binding
          real(8), dimension(*), intent(in) :: var1
          integer, intent(in), value :: nn,mm
          type(c_ptr) :: wrap_double_array_2d
        end function wrap_double_array_2d
      end interface
   
      n = size(arr,1)
      m = size(arr,2)
      yarr = wrap_double_array_2d(arr, n, m)
     
   end subroutine ocn_yakl_wrap_array_d2d

!***********************************************************************
!
!  routine ocn_yakl_wrap_array_d2d
!
!> \brief  initialize yakl library
!> \author Rob Aulwes
!> \date   October 2021
!> \details
!>  This routine initializes the YAKL library
!
!-----------------------------------------------------------------------
   subroutine ocn_yakl_unwrap_array_d1d(yarr, arr)!{{{
   
     type(c_ptr), intent(in) :: yarr
     real (kind=RKIND), dimension(:), intent(in) :: arr

     interface
       subroutine unwrap_1d(var1) bind(C,name="unwrap_1d")
          use iso_c_binding
          type(c_ptr), value :: var1
        end subroutine
      end interface
   
      call unwrap_1d(yarr)
     
   end subroutine ocn_yakl_unwrap_array_d1d

!***********************************************************************
!
!  routine ocn_yakl_wrap_array_d2d
!
!> \brief  initialize yakl library
!> \author Rob Aulwes
!> \date   October 2021
!> \details
!>  This routine initializes the YAKL library
!
!-----------------------------------------------------------------------
   subroutine ocn_yakl_unwrap_array_d2d(yarr, arr)!{{{
   
      type(c_ptr), intent(in) :: yarr
      real (kind=RKIND), dimension(:,:), intent(in) :: arr

     interface
       subroutine unwrap_2d(var1) bind(C,name="unwrap_2d")
          use iso_c_binding
          type(c_ptr), value :: var1
        end subroutine
      end interface
   
      call unwrap_2d(yarr)
     
   end subroutine ocn_yakl_unwrap_array_d2d

!***********************************************************************
!
!  routine ocn_yakl_update_host_d2d
!
!> \brief  initialize yakl library
!> \author Rob Aulwes
!> \date   October 2021
!> \details
!>  This routine initializes the YAKL library
!
!-----------------------------------------------------------------------
   subroutine ocn_yakl_update_host_d2d(yarr, arr)!{{{
   
      real (kind=RKIND), dimension(:,:), intent(inout) :: arr
      type(c_ptr), intent(out) :: yarr
      
     interface
       subroutine update_host_d2d(var1, arr) bind(C,name="update_host_d2d")
          use iso_c_binding
          real(8), dimension(*), intent(in) :: arr
          type(c_ptr), value :: var1
        end subroutine update_host_d2d
      end interface
   
      call update_host_d2d(yarr, arr)
     
   end subroutine ocn_yakl_update_host_d2d

!***********************************************************************
!
!  routine ocn_yakl_update_host_d1d
!
!> \brief  initialize yakl library
!> \author Rob Aulwes
!> \date   October 2021
!> \details
!>  This routine initializes the YAKL library
!
!-----------------------------------------------------------------------
   subroutine ocn_yakl_update_host_d1d(yarr, arr)!{{{
   
      real (kind=RKIND), dimension(:), intent(inout) :: arr
      type(c_ptr), intent(out) :: yarr
      
     interface
       subroutine update_host_d1d(var1, arr) bind(C,name="update_host_d1d")
          use iso_c_binding
          real(8), dimension(*), intent(in) :: arr
          type(c_ptr), value :: var1
        end subroutine update_host_d1d
      end interface
   
      call update_host_d1d(yarr, arr)
     
   end subroutine ocn_yakl_update_host_d1d

!***********************************************************************
!
!  routine ocn_yakl_update_device_d1d
!
!> \brief  initialize yakl library
!> \author Rob Aulwes
!> \date   October 2021
!> \details
!>  This routine initializes the YAKL library
!
!-----------------------------------------------------------------------
   subroutine ocn_yakl_update_device_d1d(yarr, arr)!{{{
   
      real (kind=RKIND), dimension(:), intent(in) :: arr
      type(c_ptr), intent(out) :: yarr
      
     interface
       subroutine update_device_d1d(var1, arr) bind(C,name="update_device_d1d")
          use iso_c_binding
          real(8), dimension(*), intent(in) :: arr
          type(c_ptr), value :: var1
        end subroutine update_device_d1d
      end interface
   
      call update_device_d1d(yarr, arr)
     
   end subroutine ocn_yakl_update_device_d1d

!***********************************************************************
!
!  routine ocn_yakl_update_host_d2d
!
!> \brief  initialize yakl library
!> \author Rob Aulwes
!> \date   October 2021
!> \details
!>  This routine initializes the YAKL library
!
!-----------------------------------------------------------------------
   subroutine ocn_yakl_update_device_d2d(yarr, arr)!{{{
   
      real (kind=RKIND), dimension(:,:), intent(in) :: arr
      type(c_ptr), intent(out) :: yarr
      
     interface
       subroutine update_device_d2d(var1, arr) bind(C,name="update_device_d2d")
          use iso_c_binding
          real(8), dimension(*), intent(in) :: arr
          type(c_ptr), value :: var1
        end subroutine update_device_d2d
      end interface
   
      call update_device_d2d(yarr, arr)
     
   end subroutine ocn_yakl_update_device_d2d

!***********************************************************************

end module ocn_yakl
