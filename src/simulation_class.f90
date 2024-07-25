module simulation_class
    implicit none
    private

    public :: part
    public :: simulation

    type :: part
        real(8) :: id         !> Particle id number
        real(8) :: pos(3)     !> Particle position
        real(8) :: vel(3)     !> Modeled velocity
    end type part

    type :: simulation
        integer :: npart                       !> Number of particles
        type(part), allocatable :: ps(:)       !> Vector containing particles
        integer :: numbins, rank
        real(8) :: delta
        real(8) :: length
        real(8) :: t, tf, dt
        integer(4) :: step, stepf
        logical :: isRun=.true.
      contains
        procedure :: init
        procedure :: increment_time        !> Step forward particle SDEs/ODEs in time
        procedure :: check_time
        procedure :: check_and_correct_bounds
        procedure :: write_particle_data   !> Routine to write particle information to text
        procedure :: print
    end type simulation

contains

   subroutine init(self, npart, length, numbins, delta, t, tf, dt, step, stepf, rank) 
      implicit none
      class(simulation), intent(out) :: self
      integer, intent(in) :: npart, numbins, stepf, step, rank
      real(8), intent(in) :: length, tf, dt, delta, t
      integer :: i

      self%rank = rank
      self%npart = npart
      self%numbins = numbins
      self%length = length
      self%delta = delta
      
      call random_initialize

      allocate(self%ps(1:self%npart))
      if (.not. allocated(self%ps)) then
         print *, 'Error: Allocation of self%ps failed.'
         stop
      end if
      
      do i=1,self%npart
         self%ps(i)%pos = get_rand_pos(self%length)
         self%ps(i)%vel = 0.0
         self%ps(i)%id  = real(i + rank*npart, 8)
      end do

      self%t = t
      self%tf = tf
      self%dt = dt
      self%step = step
      self%stepf = stepf
   end subroutine init

   function get_rand_pos(length) result(pos)
      real(8), intent(in) :: length
      real(8), dimension(3) :: pos
      logical :: isInSphere
      integer :: i

      isInSphere=.false.
      pos = 0.0

      do while (.not.isInSphere)      
         do i=1,3
            call random_number(pos(i))
            pos(i) = 2.0 * pos(i) - 1.0
            pos(i) = pos(i) * length
         end do
         pos = pos * length
         if (norm2(pos).le.length) then
            isInSphere = .true.
         end if
      end do
   end function get_rand_pos

   function kernel(pos, delta) result(rho)
      real(8), dimension(3), intent(in) :: pos
      real(8), intent(in) :: delta
      real(8) :: rho
      rho = EXP(-0.5 * norm2(pos)**2 / delta**2)
   end function kernel

   !> Initialization of the RNG: seeding is based on /dev/urandom
   !> or if not available on a naive RNG seeded with time and pid.
   !>
   !> This comes from the GFortran website.
   subroutine random_initialize
      use, intrinsic :: iso_fortran_env, only: INT32,INT64
      implicit none
      integer(kind=INT32), allocatable, dimension(:) :: seed
      integer(kind=INT32) :: i,n,un,istat,dt(8),pid
      integer(kind=INT64) :: t
      ! Get seed size
      call random_seed(size=n)
      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=un,file="/dev/urandom",access="stream",form="unformatted",action="read",status="old",iostat=istat)
      if (istat.eq.0) then
         read(un) seed
         close(un)
      else
         ! Fallback to XOR:ing the current time and pid. The PID is useful in
         ! case one launches multiple instances of the same program in parallel
         call system_clock(t)
         if (t.eq.0) then
            call date_and_time(values=dt)
            t=(dt(1)-1970)*365_INT64*24*60*60*1000 &
            & +dt(2)*31_INT64*24*60*60*1000+dt(3)*24_INT64*60*60*1000 &
            & +dt(5)*60*60*1000+dt(6)*60*1000+dt(7)*1000+dt(8)
         end if
         pid=getpid()
         t=ieor(t,int(pid,kind(t)))
         do i=1,n
            seed(i)=lcg(t)
         end do
      end if
      ! Place the seed
      call random_seed(put=seed)
   contains
      !> Simple RNG, sufficient for seeding a better RNG
      function lcg(s) result(v)
         integer(kind=INT64) :: s
         integer(kind=INT32) :: v
         if (s.eq.0) then
            s=104729
         else
            s=mod(s,4294967296_INT64)
         end if
         s=mod(s*279470273_INT64,4294967291_INT64)
         v=int(mod(s,int(huge(0),INT64)),kind(0))
      end function lcg
   end subroutine random_initialize

   !> Sampling of a WP real from a normal distribution
   !> Adapted from the following Fortran 77 code:
   !> ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM,
   !> IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
   !> VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
   !>
   !> The function random_normal() returns a normally distributed
   !> pseudo-random number with zero mean and unit variance.
   !> The algorithm uses the ratio of uniforms method of A.J. Kinderman
   !> and J.F. Monahan augmented with quadratic bounding curves.
   function random_normal(m,sd) result(res)
      implicit none
      !> Mean and standard deviation of the distribution
      real(8), intent(in), optional :: m,sd
      real(8) :: res
      real(8), parameter :: s =+0.449871
      real(8), parameter :: t =-0.386595
      real(8), parameter :: a =+0.19600
      real(8), parameter :: b =+0.25472
      real(8), parameter :: r1=+0.27597
      real(8), parameter :: r2=+0.27846
      real(8) :: u,v,x,y,q
      ! Generate P=(u,v) uniform in rectangle enclosing acceptance region
      try: do
         call random_number(u)
         call random_number(v)
         v=1.7156*(v-0.5)
         ! Evaluate the quadratic form
         x=u-s
         y=abs(v)-t
         q=x**2+y*(a*y-b*x)
         ! Accept P if inside inner ellipse
         if (q.lt.r1) exit try
         ! Reject P if outside outer ellipse
         if (q.gt.r2) cycle try
         ! Accept P if inside acceptance region
         if (v**2.lt.-4.0*log(u)*u**2) exit try
      end do try
      ! Return ratio of P's coordinates as the normal deviate
      res=v/u
      ! Modify to give correct mean and standard deviation
      if (present(sd)) res=res*sd
      if (present(m))  res=res+m
   end function random_normal

   subroutine increment_time(this, a, b)
      implicit none
      class(simulation), intent(inout) :: this
      type(part) :: p
      real(8), intent(in) :: a, b
      real(8), dimension(3) :: dWi, dWj
      real(8) :: rho, mean
      integer :: i
      mean=0.0

      do i=1,this%npart
         p = this%ps(i)
         rho = kernel(p%pos, this%delta)
         dWi = [random_normal(m=mean,sd=sqrt(this%dt)), random_normal(m=mean,sd=sqrt(this%dt)), random_normal(m=mean,sd=sqrt(this%dt))]
         dWj = [random_normal(m=mean,sd=sqrt(this%dt)), random_normal(m=mean,sd=sqrt(this%dt)), random_normal(m=mean,sd=sqrt(this%dt))]
         p%vel = (1.0 - a*this%dt)*p%vel + b*((1.0 - rho)*dWi + (rho - 1.0)*dWj)/sqrt(1.0 + rho**2)
         p%pos = p%pos + p%vel*this%dt
         this%ps(i) = p
      end do

      this%t = this%t + this%dt
      this%step = this%step + 1
   end subroutine increment_time

   subroutine check_and_correct_bounds(this)
      implicit none
      class(simulation), intent(inout) :: this
      real(8) :: r
      integer :: i

      do i=1,this%npart
         r = norm2(this%ps(i)%pos)
         if (r.gt.this%length) then
            this%ps(i)%pos = get_rand_pos(this%length)
         end if
      end do
   end subroutine check_and_correct_bounds

   subroutine check_time(this)
      implicit none
      class(simulation), intent(inout) :: this
      
      if (this%t.ge.this%tf) this%isRun = .false.
      if (this%step.ge.this%stepf) this%isRun = .false.
   end subroutine check_time

   subroutine write_particle_data(this, filename)
      use mpi
      implicit none
      class(simulation), intent(inout) :: this
      character(len=*), intent(in) :: filename
      integer :: status(MPI_STATUS_SIZE)
      integer :: i,io, ierr, particle_size
      integer(8) :: disp
      integer :: unit
      
      ! Get size of one particle
      particle_size = 7*8
      ! Calc displacement for parallel write
      disp = this%rank*this%npart*particle_size
      
      ! Cautionary
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! Do the I/O
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, io, ierr)
      call MPI_FILE_WRITE_AT_ALL(io, disp, this%ps, this%npart*7, MPI_DOUBLE_PRECISION, status, ierr)
      call MPI_FILE_CLOSE(io, ierr)
  end subroutine write_particle_data

  subroutine print(this)
   implicit none
   class(simulation), intent(in) :: this

   print *, "Time :: ", this%t, "  ", "Step :: ", this%step
  end subroutine print

end module simulation_class
