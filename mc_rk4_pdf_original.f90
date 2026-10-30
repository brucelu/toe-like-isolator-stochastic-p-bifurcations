program mc_rk4_pdf_original
! Original system MC (Ito Euler–Maruyama) + Poincaré sampling
  ! Goal: statistics consistent with slow-flow / FPK in (a,phi)
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)
  real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp

  !================== MC / grid settings ==================
  integer, parameter :: Nx1=1, Nx2=1                 ! initial-condition grid
  integer, parameter :: Nc = Nx1*Nx2
  integer, parameter :: Nperiod=100000               ! total sampled periods (recommend 1e5~5e5)
  integer, parameter :: Nt=5000                     ! 丢弃的瞬态周期数
  integer, parameter :: Na=200, Nphi=480         ! 幅度a和相位phi的直方图 bin 数



  real(dp) :: xlow, xup, vlow, vup
  real(dp) :: a_min, a_max, da, dphi
  real(dp) :: hx, hv                               !初始条件网格的步长
  real(dp) :: dt,Tper,  total_time, transient_time  
     !时间步长,激励周期,总模拟时间和瞬态时间

  !integer  ::nstep_total !总时间步数
  ! real(dp) :: Tper, dt
  integer  :: nstep_per_period

  !================== system parameters (Eq.7) ==================
  real(dp) :: k1,k2,k3,k4, b0,b1,b2, xi1,xi2, rho, Z, g, S_eta, wn, eps,sigma,Omega

  !================== MC variables ==================
  integer :: dx,dy, ip, istep!初始条件网格索引,时间步索引
  real(dp) :: t, x, v
  real(dp) :: U1,U2, rn
  real(dp) :: driftv, dW, noiseAmp
  real(dp) :: a, phi
  !real(dp) :: a, theta, phi
  integer :: ia, jphi
  integer :: sample_count
  real(dp), allocatable :: H(:,:)
!integer :: sample_stride
!sample_stride = 50
  ! ---- continuous-sample statistics for a (before histogram binning effects) ----
integer :: n_a_stat
real(dp) :: mean_a_mc, M2_a_mc, var_a_mc, delta_a
integer :: low_count, high_count
real(dp) :: a_split
integer :: current_state, previous_state
integer :: transition_count


! ===== Variables for one-period lock-in extraction =====
real(dp) :: Cfund, Sfund
real(dp) :: sum_cos, sum_sin
real(dp) :: x_old, v_old
real(dp) :: x_new, v_new
real(dp) :: x_mid, t_mid

! Declaration
!real(dp) :: phi
real(dp) :: phi_raw
real(dp) :: phi_corrected
real(dp) :: delta_phi_det
  !================== helpers ==================
 !real(dp) :: dphi_align   ! 在变量声明区添加

!dphi_align = 0.0_dp!0.383972435439_dp   ! 约 22 度（你当前拟合得到）
a_split = 2.8_dp
delta_phi_det = 0.0505_dp!0.0435_dp
  call random_seed()

  ! ---- parameters (use your values) ----
  k1=0.0313_dp!0.03124999998_dp!0.0313_dp!0.03124999998_dp
  k2=-0.0609_dp!-0.06089241115_dp!-0.0609_dp
  k3=0.0381_dp!0.03808593748_dp!0.0381_dp
  k4=-0.0075_dp!-0.0166226915726781_dp!-0.01_dp

  b0=0.5_dp
  b1=-0.2887_dp
  b2=0.5417_dp

  xi1=0.005_dp
  xi2=0.005_dp
  rho = xi1 + xi2*b0

  Z=0.001_dp   !更容易对齐的量级是 Z ≈ 1e-3,(1e-2不容易对的上)。
  g=0.001_dp
  S_eta=0.005_dp

  wn = sqrt(k1)
  ! IMPORTANT: use the SAME Omega as in modulation (primary resonance scaling)
  eps=0.1_dp;
  !sigma = -1.0_dp;               
!Omega = wn + eps**2*sigma ;  
  Omega=0.1669_dp!0.1543_dp
  sigma=(Omega-wn)/(eps*eps)
  !Omega = wn + (0.1_dp**2)*(-1.0_dp)   ! eps=0.1, sigma=-1

  ! ---- IC domain (x,v) ----
  xlow=-0.5_dp; xup= 0.5_dp
  vlow=-0.1_dp; vup= 0.1_dp
  hx=(xup-xlow)/real(Nx1,dp)
  hv=(vup-vlow)/real(Nx2,dp)

  ! ---- (a,phi) histogram domain ----
  a_min = 2.7_dp
  a_max = 3._dp
  da    = (a_max-a_min)/real(Na,dp)
  dphi  = (2.0_dp*PI)/real(Nphi,dp)

  allocate(H(Na,Nphi))
  H = 0.0_dp
  sample_count = 0

  !新增的
  n_a_stat  = 0
mean_a_mc = 0.0_dp
M2_a_mc   = 0.0_dp
var_a_mc  = 0.0_dp
  
  ! ---- time step ----
  Tper = 2.0_dp*PI/Omega   !一个激励周期
  dt   = 0.005_dp          !小时间步，确保数值稳定性
  ! total_time = 200.0_dp * Tper ! 总模拟时间（覆盖200周期）
!nstep_total = int(total_time/dt)!总步数≈200*Tper/0.001（很大，计算密集）
   ! transient_time = Nt*Tper
  nstep_per_period = int(Tper/dt)
   if (nstep_per_period < 20) nstep_per_period = 20
  dt = Tper/real(nstep_per_period,dp)   ! adjust dt to fit integer steps

  ! noise amplitude in dv:  g*sqrt(S_eta)*sqrt(dt)*N(0,1)
  noiseAmp = g*sqrt(S_eta)

 do dx=1,Nx1
    do dy=1,Nx2
      low_count = 0
      high_count = 0
      previous_state = -1
      current_state = -1
      transition_count = 0

      t = 0.0_dp
      x = xlow + (real(dx,dp)-0.5_dp)*hx
      v = vlow + (real(dy,dp)-0.5_dp)*hv

      do ip=1,Nperiod

  ! Reset one-period lock-in integrals
  sum_cos = 0.0_dp
  sum_sin = 0.0_dp

  do istep=1,nstep_per_period

    ! ------------------------------------------------------------
    ! Generate one standard Gaussian random number
    ! ------------------------------------------------------------
    call random_number(U1)
    call random_number(U2)

    if (U1 < 1.0e-12_dp) U1 = 1.0e-12_dp

    rn = sqrt(-2.0_dp*log(U1))*cos(2.0_dp*PI*U2)
    dW = sqrt(dt)*rn

    ! Save old state
    x_old = x
    v_old = v

    ! ------------------------------------------------------------
    ! Original-system drift evaluated at the old state
    ! ------------------------------------------------------------
    driftv = -k1*x_old                              &
             -k2*x_old*x_old                       &
             -k3*x_old*x_old*x_old                 &
             -k4*x_old*x_old*x_old*x_old           &
             -2.0_dp*rho*v_old                     &
             -2.0_dp*xi2*                          &
              (b1*x_old+b2*x_old*x_old)*v_old      &
             +Z*cos(Omega*t)

    ! ------------------------------------------------------------
    ! Euler-Maruyama update
    ! ------------------------------------------------------------
    x_new = x_old+v_old*dt

    v_new = v_old+driftv*dt+noiseAmp*dW

    ! ------------------------------------------------------------
    ! Midpoint approximation for the Fourier integrals
    ! ------------------------------------------------------------
    t_mid = t+0.5_dp*dt
    x_mid = 0.5_dp*(x_old+x_new)

    sum_cos = sum_cos +                         &
              x_mid*cos(Omega*t_mid)*dt

    sum_sin = sum_sin +                         &
              x_mid*sin(Omega*t_mid)*dt

    ! Accept updated state
    x = x_new
    v = v_new
    t = t+dt

  end do  ! istep

  ! ==============================================================
  ! Extract the fundamental slow amplitude and phase
  ! ==============================================================
  Cfund = (2.0_dp/Tper)*sum_cos
  Sfund = (2.0_dp/Tper)*sum_sin

  if (ip > Nt) then

    a = sqrt(Cfund*Cfund+Sfund*Sfund)/eps
    
    ! Phase extraction and coordinate correction
phi_raw = modulo(atan2(Sfund,Cfund),2.0_dp*PI)

phi_corrected = modulo(                         &
    phi_raw-delta_phi_det+2.0_dp*PI,            &
    2.0_dp*PI)

phi = phi_corrected
    
    if (a<a_split) then
         current_state = 0
  low_count = low_count+1
    else
        current_state = 1
  high_count = high_count+1
    end if
    
    if (previous_state >= 0) then
    if (current_state /= previous_state) then
        transition_count = transition_count + 1
    end if
end if

previous_state = current_state

    if (a>=a_min .and. a<a_max) then

      ia = int((a-a_min)/da)+1
      jphi = int(phi/dphi)+1

      if (jphi<1) jphi=1
      if (jphi>Nphi) jphi=Nphi

      if (ia>=1 .and. ia<=Na) then

        H(ia,jphi) = H(ia,jphi)+1.0_dp
        sample_count = sample_count+1

        ! Continuous-sample statistics
        n_a_stat = n_a_stat+1
        delta_a = a-mean_a_mc
        mean_a_mc = mean_a_mc +                     &
                    delta_a/real(n_a_stat,dp)
        M2_a_mc = M2_a_mc +                         &
                  delta_a*(a-mean_a_mc)

      end if
    end if
    

  end if

  if (ip==Nt) then
    print*,'Transient finished.'
  end if

  if (mod(ip,10000)==0) then
    print*,'ip=',ip,' sample_count=',sample_count
  end if

      end do  ! ip
      
      if (low_count + high_count > 0) then
        print *, 'trajectory ',dx,dy,                    &
                 ' P_high = ',                           &
                 real(high_count,dp) /                   &
                 real(low_count+high_count,dp),          &
                 ' transitions = ',transition_count
      else
        print *, 'trajectory ',dx,dy,' has no valid samples.'
      end if

    end do
  end do

     
  ! ---- normalize to pdf ----
  if (sample_count <= 0) then
    print*, 'No samples counted. Check a-range or parameters.'
   stop
  end if

  !归一化：先除以总采样数（得到频率），再除以bin面积(da*dphi)，得到密度（PDF）
  H = H / real(sample_count,dp) 
  H = H / (da*dphi)    ! density

  call write_outputs(H, Na, Nphi, a_min, da, dphi)
print*, 'Done. sample_count=', sample_count
!新增的
if (n_a_stat > 1) then
  var_a_mc = M2_a_mc / real(n_a_stat - 1, dp)
else
  var_a_mc = 0.0_dp
end if

print*, 'MC continuous-sample stats (a):'
print*, '  n_a_stat = ', n_a_stat
print*, '  mean(a)  = ', mean_a_mc
print*, '  var(a)   = ', var_a_mc
print*, '  std(a)   = ', sqrt(max(var_a_mc,0.0_dp))
 !print*, 'Continuous-time sampling done. samples=', sample_count
contains

  subroutine write_outputs(H,Na,Nphi,a_min,da,dphi)
  integer, intent(in) :: Na,Nphi
  real(dp), intent(in) :: H(Na,Nphi),a_min,da,dphi

  integer :: i,j
  real(dp) :: ac,phic
  real(dp) :: pa(Na),pphi(Nphi)

  pa   = 0.0_dp
  pphi = 0.0_dp

  ! Amplitude marginal density
  do i=1,Na
    do j=1,Nphi
      pa(i) = pa(i)+H(i,j)*dphi
    end do
  end do

  ! Phase marginal density
  do j=1,Nphi
    do i=1,Na
      pphi(j) = pphi(j)+H(i,j)*da
    end do
  end do

  open(10,file='MC_PDF_aphi.txt',status='replace')
  do i=1,Na
    do j=1,Nphi
      write(10,'(ES24.16)',advance='no') H(i,j)
    end do
    write(10,*)
  end do
  close(10)

  open(11,file='MC_pa.txt',status='replace')
  do i=1,Na
    ac = a_min+(real(i,dp)-0.5_dp)*da
    write(11,'(2ES24.16)') ac,pa(i)
  end do
  close(11)

  open(12,file='MC_pphi.txt',status='replace')
  do j=1,Nphi
    phic = (real(j,dp)-0.5_dp)*dphi
    write(12,'(2ES24.16)') phic,pphi(j)
  end do
  close(12)

  print*,'Integral p(a)   = ',sum(pa)*da
  print*,'Integral p(phi) = ',sum(pphi)*dphi

end subroutine write_outputs

end program mc_rk4_pdf_original