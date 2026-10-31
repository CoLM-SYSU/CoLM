program diagnostics
   use MOD_Precision
   use MOD_MT_HydroDiagnostics
   implicit none
   type(mt_water_diagnostics) :: a,b
   call a%reset(4,-3)
   if(lbound(a%etroot,1)/=1.or.ubound(a%etroot,1)/=4) error stop 'soil bounds'
   if(lbound(a%snow_qout_layer,1)/=-3.or.ubound(a%snow_qout_layer,1)/=0) error stop 'snow bounds'
   a%etroot=1; a%etroot_actual=2; a%rsub_source_layer=3
   a%wblc_ice_sink=4; a%snow_qout_layer=5; a%permeable_soil=.true.
   a%etroot_aquifer=6; a%etroot_surface=7; a%rsub_source_surface=8; a%rsub_source_aquifer=9
   a%dew_overflow=10; a%frost_displaced=11; a%late_runoff=12
   a%imperv_evap_wdsrf=13; a%imperv_evap_soil=14; a%imperv_subl_soil=15; a%qgtop=16
   b=a
   call a%reset(4,-3)
   if(any(a%etroot/=0).or.any(a%etroot_actual/=0).or.any(a%rsub_source_layer/=0)) error stop 'soil reset'
   if(any(a%wblc_ice_sink/=0).or.any(a%snow_qout_layer/=0).or.any(a%permeable_soil)) error stop 'water reset'
   if(any([a%etroot_aquifer,a%etroot_surface,a%rsub_source_surface,a%rsub_source_aquifer, &
      a%dew_overflow,a%frost_displaced,a%late_runoff,a%imperv_evap_wdsrf, &
      a%imperv_evap_soil,a%imperv_subl_soil,a%qgtop]/=0)) error stop 'scalar reset'
   if(any(b%etroot/=1).or.any(b%wblc_ice_sink/=4).or.any(b%snow_qout_layer/=5).or.b%qgtop/=16) &
      error stop 'reset changed another object'
   call b%reset_soil(2)
   if(size(b%etroot)/=2.or.size(b%etroot_actual)/=2.or.size(b%rsub_source_layer)/=2) error stop 'soil resize'
   if(any(b%etroot/=0).or.any(b%etroot_actual/=0).or.any(b%rsub_source_layer/=0)) error stop 'resized reset'
   if(any(b%wblc_ice_sink/=4).or.any(b%snow_qout_layer/=5).or.b%qgtop/=16) error stop 'soil reset clobber'
   call b%reset(2,0)
   if(lbound(b%snow_qout_layer,1)/=0.or.ubound(b%snow_qout_layer,1)/=0) error stop 'thin snow bounds'
   if(size(b%wblc_ice_sink)/=2.or.size(b%permeable_soil)/=2) error stop 'water resize'
   call b%reset(1,1)
   if(size(b%snow_qout_layer)/=0) error stop 'no snow extent'
   call b%reset(3,-2)
   if(lbound(b%snow_qout_layer,1)/=-2.or.ubound(b%snow_qout_layer,1)/=0) error stop 'snow returned'
   if(any(b%snow_qout_layer/=0).or.any(b%etroot/=0)) error stop 'stale state after resize'
   print *, 'PASS diagnostics: reset, resize, snow bounds, independent objects and partial reset'
end program
