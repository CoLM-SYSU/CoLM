program soil_column
   use MOD_Precision
   use MOD_Namelist
   use MOD_Hydro_SoilWater
   use MOD_Hydro_SoilFunction
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifndef LEGACY_INTERFACE
   use MOD_MT_HydroDiagnostics, only: mt_water_diagnostics
#endif
   implicit none
   integer, parameter :: n=4, np=5
   integer :: case_id, j, step
   real(r8) :: z(n), zi(0:n), por(n), residual(n), psi(n), ks(n), prms(np,n)
   real(r8) :: dt, qtop, etr, roots(n), rootflux(n), runoff, infil, pond, zwt, wa
   real(r8) :: liquid(n), mass(n), smp(n), hk(n), qlayer(0:n), balance
   real(r8) :: root_demand(n), root_actual(n), root_aquifer, root_surface
   real(r8) :: runoff_layer(n), runoff_aquifer, runoff_surface
   logical :: permeable(n)
#ifndef LEGACY_INTERFACE
   type(mt_water_diagnostics) :: diag
#endif
   open(10,file='column.bin',access='stream',form='unformatted',status='replace')
   zi = [0._r8,100._r8,300._r8,700._r8,1500._r8]
   z = (zi(:n-1)+zi(1:))/2._r8
   por=0.45_r8; residual=0.05_r8; psi=-100._r8; ks=0.01_r8
   prms(1,:)=0.002_r8; prms(2,:)=1.5_r8; prms(3,:)=0.5_r8
   do j=1,n
      call get_derived_parameters_vGM(psi(j),prms(1,j),prms(2,j),prms(4,j),prms(5,j))
   enddo
   dt=1800._r8
   do case_id=1,6
      permeable=.true.
      roots=[0.4_r8,0.3_r8,0.2_r8,0.1_r8]
      qtop=0.001_r8; etr=0.0002_r8; runoff=0.0001_r8
      pond=0._r8; zwt=2000._r8
      DEF_USE_PLANTHYDRAULICS=.false.
      select case(case_id)
      case(2)
         zwt=50._r8; pond=3._r8
      case(3)
         permeable(2)=.false.
      case(4)
         permeable=.false.
         qtop=0._r8
      case(5)
         DEF_USE_PLANTHYDRAULICS=.true.
         rootflux=[-0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8]
      case(6)
         zwt=5000._r8; etr=0.02_r8; qtop=0._r8
      end select
      call get_water_equilibrium_state(zwt,n,mass,smp,hk,wa,z,zi,por,residual,psi,ks,np,prms)
      liquid=mass/(zi(1:)-zi(:n-1))
      if (.not.DEF_USE_PLANTHYDRAULICS) rootflux=roots*etr
#ifndef LEGACY_INTERFACE
      call diag%reset(n,-3)
#endif
      do step=1,3
#ifndef LEGACY_INTERFACE
         ! Lower-level soil resets must preserve diagnostics captured by WATER_VSF.
         diag%snow_qout_layer=17._r8
         diag%wblc_ice_sink=23._r8
         diag%qgtop=29._r8
#endif
         call soil_water_vertical_movement(n,dt,z,zi,permeable,por,residual,psi,ks,np,prms,por(n), &
            qtop,etr,roots,rootflux,runoff,infil,pond,zwt,wa,liquid,smp,hk,qlayer, &
#ifdef LEGACY_INTERFACE
            root_demand,root_actual,root_aquifer,root_surface,runoff_layer,runoff_surface,runoff_aquifer, &
#else
            diag, &
#endif
            1.e-3_r8,balance)
#ifndef LEGACY_INTERFACE
         root_demand=diag%etroot; root_actual=diag%etroot_actual
         root_aquifer=diag%etroot_aquifer; root_surface=diag%etroot_surface
         runoff_layer=diag%rsub_source_layer
         runoff_surface=diag%rsub_source_surface; runoff_aquifer=diag%rsub_source_aquifer
         if(any(diag%snow_qout_layer/=17._r8)) error stop 'soil reset lost snow flux'
         if(any(diag%wblc_ice_sink/=23._r8)) error stop 'soil reset lost ice sink'
         if(diag%qgtop/=29._r8) error stop 'soil reset lost top flux'
#endif
         if (.not.all(ieee_is_finite([liquid,smp,hk,qlayer,infil,pond,zwt,wa,balance, &
            root_demand,root_actual,root_aquifer,root_surface,runoff_layer,runoff_surface,runoff_aquifer]))) &
            error stop 'non-finite water state or diagnostic'
         if(abs(balance)>1.e-3_r8) error stop 'column water balance'
         if(case_id==3.and.root_actual(2)/=0._r8) error stop 'root uptake from impermeable layer'
         if(case_id==5.and.root_actual(1)>=0._r8) error stop 'root return was not exercised'
         if(case_id==4.and.root_aquifer<=0._r8) error stop 'aquifer uptake was not exercised'
         if(case_id==2.and.step>1.and.root_surface<=0._r8) error stop 'pond uptake was not exercised'
         if(abs(sum(root_actual)+root_aquifer+root_surface-etr*dt)>1.e-10_r8) &
            error stop 'root pathway budget'
         if(abs(sum(runoff_layer)+runoff_aquifer+runoff_surface-runoff*dt)>1.e-10_r8) &
            error stop 'runoff pathway budget'
         write(10) liquid,smp,hk,qlayer,infil,pond,zwt,wa,balance, &
            root_demand,root_actual,root_aquifer,root_surface,runoff_layer,runoff_surface,runoff_aquifer
      enddo
   enddo
   close(10)
   print *, 'PASS 6 soil columns x 3 steps: finite states, water balance, root pathways'
end program
