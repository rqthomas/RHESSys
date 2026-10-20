! Minimal stand-in for CESM share shr_infnan_mod. FATES only does
!   use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
!   use shr_infnan_mod, only : isnan => shr_infnan_isnan
! so shr_infnan_nan is a quiet-NaN real(r8) constant (intrinsic assignment
! then covers scalars and arrays of any rank) and assignment(=) is a token
! generic that only has to exist to satisfy the import.
module shr_infnan_mod
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use, intrinsic :: iso_fortran_env, only : int64
   use shr_kind_mod, only : r8 => shr_kind_r8
   implicit none
   private
   public :: shr_infnan_nan, shr_infnan_isnan, assignment(=)

   real(r8), parameter :: shr_infnan_nan = &
        transfer(int(z'7FF8000000000000', int64), 1.0_r8)

   type :: shr_infnan_unused_type
      integer :: i = 0
   end type shr_infnan_unused_type

   interface assignment(=)
      module procedure assign_unused
   end interface

   interface shr_infnan_isnan
      module procedure isnan_r8
   end interface
contains
   elemental logical function isnan_r8(x)
      real(r8), intent(in) :: x
      isnan_r8 = ieee_is_nan(x)
   end function isnan_r8
   subroutine assign_unused(a, b)
      type(shr_infnan_unused_type), intent(out) :: a
      integer,                      intent(in)  :: b
      a%i = b
   end subroutine assign_unused
end module shr_infnan_mod
