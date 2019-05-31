program main
  implicit none
  integer,parameter    :: kd = kind(1.0)
  real(kd),parameter   :: pi = 4.0_kd*atan(1.0_kd)
  real(kd),allocatable :: R(:)              ! 自相关函数值
  real(kd),allocatable :: S(:)              ! 波能谱值
  real(kd),allocatable :: omg(:)            ! 频率
  real(kd),allocatable :: zeta(:)           ! 波面升高
  integer              :: n                 ! 对时间 N 等分
  integer              :: m                 ! 滞后波浪数量
  real(kd)             :: T                 ! 总时间 s
  real(kd)             :: tp(2)             ! 时间节点
  real(kd)             :: dt                ! 时间间隔

  real(kd)             :: tmp
  logical              :: exists
  integer              :: i,k,ierr 
  character(128)        :: filename
 

    filename = "random_wave.csv" 
    open(10,file=trim(adjustl(filename)))
    n = 0
    do
      read(10,*,iostat = ierr )
      if(ierr .ne. 0 ) exit
      n = n + 1
    end do
    close(10) 
    m = n/37 ! 选择合适的值
    allocate( R(0:m),S(0:m),omg(0:m),zeta(n) )
    R = 0.0_kd; S = 0.0_kd; omg = 0.0_kd; zeta = 0.0_kd


    ! 读取波面升高
    open(10,file=trim(adjustl(filename)))
    read(10,*) tp(1),zeta(1)
    do i = 2, n-1
      read(10,*) tmp,zeta(i)
    end do
    read(10,*) tp(2),zeta(n)
    close(10)
    T  = tp(2) - tp(1)
    dt = T/n 
      
    ! 获取各 k 值处的自相关函数
    do k = 0, m
      do i = 1, n-k
        R(k) = R(k) + zeta(i)*zeta(i+k)
      end do
      R(k) = R(k) / real(n-k,kd)
    end do
    
    ! 利用梯形积分获得谱值
    do i = 0, m
      omg(i) = real(i,kd)*pi/m/dt
      do k = 1, m-1
        S(i) = S(i) + 2.0_kd*R(k)*cos( (real(i,kd)*k)/m*pi )
      end do
      S(i) = S(i) + R(0) + R(m)*cos( real(i,kd)*pi )
    end do
    S = S * dt / pi

    ! Hamming平滑
    S(0) = 0.54_kd*S(0)+0.46_kd*S(1)
    do i = 1, m-1
      S(i) = 0.54_kd*S(i) + 0.23_kd*(S(i-1)+S(i+1))
    end do
    S(m) = 0.54_kd*S(m)+0.46_kd*S(m-1)
   ! 输出谱值
    filename = "BNP01.csv"  
    open(10,file=trim(adjustl(filename)))
    do i = 0, m
      write(10,"(f12.7,',',f12.7)") omg(i),S(i)
    end do
    close(10)
    deallocate( R,S,omg,zeta )
end program


