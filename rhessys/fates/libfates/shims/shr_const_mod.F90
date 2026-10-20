! Minimal stand-in for CESM share shr_const_mod (only what FATES uses).
module shr_const_mod
   use shr_kind_mod, only : r8 => shr_kind_r8
   implicit none
   public
   real(r8), parameter :: shr_const_cday = 86400.0_r8   ! sec in calendar day
end module shr_const_mod
