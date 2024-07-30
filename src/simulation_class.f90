module simulation_class
    implicit none
    private

    public :: part
    public :: simulation

    real(8), public :: PI = 3.1415926535897932385

    type :: part
        real(8) :: id         !> Particle id number
        real(8) :: pos(3)     !> Particle position
        real(8) :: vel(3)     !> Modeled velocity
    end type part

    type :: simulation
        integer :: npart                       !> Number of particles
        type(part), allocatable :: ps(:)       !> Vector containing particles
        integer :: numbins, rank, numproc
        real(8) :: delta
        real(8) :: length
        real(8) :: t, tf, dt
        integer(4) :: step, stepf
        logical :: isRun=.true.
        integer(4) :: particles_removed

        real(8), dimension(:), allocatable :: rdf !> Radial distribution function
        real(8), dimension(:), allocatable :: w1  !> Mean relative velocity conditioned on r
        real(8), dimension(:), allocatable :: w2  !> Second moment of rel. vel. conditioned on r
        integer, dimension(:), allocatable :: c   !> Bin counter
      contains
        procedure :: init
        procedure :: increment_time        !> Step forward particle SDEs/ODEs in time
        procedure :: check_time
        procedure :: check_and_correct_bounds
        procedure :: write_particle_data   !> Routine to write particle information to text
        procedure :: print
        procedure :: compute_rdf
        procedure :: compute_w1
        procedure :: compute_w2
        procedure :: write_rdf
        procedure :: write_w1
        procedure :: write_w2
    end type simulation

contains

   subroutine init(self, npart, length, numbins, delta, t, tf, dt, step, stepf, rank, numproc) 
      implicit none
      class(simulation), intent(out) :: self
      integer, intent(in) :: npart, numbins, stepf, step, rank, numproc
      real(8), intent(in) :: length, tf, dt, delta, t
      integer :: i

      call random_initialize

      ! MPI
      self%numproc = numproc
      self%rank = rank

      ! Domain and stats
      self%numbins = numbins
      self%length = length
      self%delta = delta
      self%particles_removed = 0

      ! Particle storage
      self%npart = npart
      allocate(self%ps(1:self%npart))
      if (.not. allocated(self%ps)) then
         print *, 'Error: Allocation of self%ps failed.'
         stop
      end if

      ! Stat storage
      allocate(self%rdf(1:self%numbins)); self%rdf=0.0
      allocate(self%w1( 1:self%numbins)); self%w1 =0.0
      allocate(self%w2( 1:self%numbins)); self%w2 =0.0
      allocate(self%c(  1:self%numbins)); self%c  =0
      
      ! Particles
      do i=1,self%npart
         self%ps(i)%pos = get_rand_pos(self%length)
         self%ps(i)%vel = 0.0
         self%ps(i)%id  = real(i + rank*npart, 8)
      end do

      ! Time tracking
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
      real(8), dimension(3,3) :: temp, bij, Q
      real(8), dimension(3) :: dWi, dWj, i1, i2, i3, rll, rt1, rt2
      real(8) :: rho_ll, rho_nn, mean, std, rll_dot, rt1_dot, rt2_dot
      integer :: i, ind

      ! Identity
      i1 = [1.0, 0.0, 0.0]
      i2 = [0.0, 1.0, 0.0]
      i3 = [0.0, 0.0, 1.0]

      ! Random normal parameters
      mean=0.0; std = sqrt(this%dt)

      do i=1,this%npart
         ! Reset arrays
         Q    = 0.0
         temp = 0.0
         bij  = 0.0

         ! Grab particle from storage
         p = this%ps(i)

         ! Get filter kernel
         rho_ll = kernel(p%pos,6.2832d0/16.0)
         rho_nn = kernel(p%pos,6.2832d0/16.0)

         ! Get dW for particle i
         dWi = [random_normal(m=mean,sd=std), & 
         &      random_normal(m=mean,sd=std), & 
         &      random_normal(m=mean,sd=std)]
         ! Get dW for particle j
         dWj = [random_normal(m=mean,sd=std), & 
         &      random_normal(m=mean,sd=std), & 
         &      random_normal(m=mean,sd=std)]

         ! Compute the reference axes
         rll = p%pos
         rll_dot = dot_product(rll,rll)+epsilon(1.0)

         ! Generate an orthonormal set based on max slip component
         ind = maxloc(abs(rll),dim=1)
         select case (ind)
         case(1)
            ! Max slip in x-direction
            rt1     = -(rll(2)/rll_dot)*rll
            rt1(2)  = 1.0 + rt1(2)
            rt1_dot = dot_product(rt1,rt1)+epsilon(1.0)
         case(2)
            ! Max slip in y-direction
            rt1     = -(rll(1)/rll_dot)*rll
            rt1(1)  = 1.0 + rt1(1)
            rt1_dot = dot_product(rt1,rt1)+epsilon(1.0)
         case(3)
            ! Max slip in z-direction
            rt1     = -(rll(1)/rll_dot)*rll
            rt1(1)  = 1.0 + rt1(1)
            rt1_dot = dot_product(rt1,rt1)+epsilon(1.0)                
         end select

         ! rt2 right-hand coordinate system (cross product)
         rt2(1) = rll(2)*rt1(3) - rll(3)*rt1(2)
         rt2(2) = rll(3)*rt1(1) - rll(1)*rt1(3)
         rt2(3) = rll(1)*rt1(2) - rll(2)*rt1(1)
         rt2_dot = dot_product(rt2,rt2)+epsilon(1.0)

         ! Normalize basis vectors
         rll = rll/sqrt(rll_dot)
         rt1 = rt1/sqrt(rt1_dot)
         rt2 = rt2/sqrt(rt2_dot)

         ! Construct rotation matrices
         Q(:,1) = rll
         Q(:,2) = rt1
         Q(:,3) = rt2

         ! Multiply Q by b^dag
         temp(1,1) = Q(1,1)*rho_ll 
         temp(2,1) = Q(2,1)*rho_ll
         temp(3,1) = Q(3,1)*rho_ll

         temp(1,2) = Q(1,2)*rho_nn 
         temp(2,2) = Q(2,2)*rho_nn
         temp(3,2) = Q(3,2)*rho_nn

         temp(1,3) = Q(1,3)*rho_nn
         temp(2,3) = Q(2,3)*rho_nn
         temp(3,3) = Q(3,3)*rho_nn

         ! Multiply Q*b^dag (temp) by Q^T  to get R tensor
         bij(1,1) = temp(1,1)*Q(1,1) + temp(1,2)*Q(1,2) + temp(1,3)*Q(1,3)
         bij(1,2) = temp(1,1)*Q(2,1) + temp(1,2)*Q(2,2) + temp(1,3)*Q(2,3)
         bij(1,3) = temp(1,1)*Q(3,1) + temp(1,2)*Q(3,2) + temp(1,3)*Q(3,3)

         bij(2,1) = temp(2,1)*Q(1,1) + temp(2,2)*Q(1,2) + temp(2,3)*Q(1,3)
         bij(2,2) = temp(2,1)*Q(2,1) + temp(2,2)*Q(2,2) + temp(2,3)*Q(2,3)
         bij(2,3) = temp(2,1)*Q(3,1) + temp(2,2)*Q(3,2) + temp(2,3)*Q(3,3)

         bij(3,1) = temp(3,1)*Q(1,1) + temp(3,2)*Q(1,2) + temp(3,3)*Q(1,3)
         bij(3,2) = temp(3,1)*Q(2,1) + temp(3,2)*Q(2,2) + temp(3,3)*Q(2,3)
         bij(3,3) = temp(3,1)*Q(3,1) + temp(3,2)*Q(3,2) + temp(3,3)*Q(3,3)

         ! Update velocity
         ! p%vel = (1.0 - a*this%dt)*p%vel + b*((1.0 - rho)*dWi + (rho - 1.0)*dWj)/sqrt(1.0 + rho**2)
         p%vel(1) = (1.0 - a*this%dt)*p%vel(1) + & 
         &          b*(dot_product(i1-bij(1,:),dWi) + dot_product(bij(1,:)-i1,dWj))/sqrt(1.0 + sum(bij(1,:)**2))
         p%vel(2) = (1.0 - a*this%dt)*p%vel(2) + & 
         &          b*(dot_product(i2-bij(2,:),dWi) + dot_product(bij(2,:)-i2,dWj))/sqrt(1.0 + sum(bij(2,:)**2))
         p%vel(3) = (1.0 - a*this%dt)*p%vel(3) + & 
         &          b*(dot_product(i3-bij(3,:),dWi) + dot_product(bij(3,:)-i3,dWj))/sqrt(1.0 + sum(bij(3,:)**2))
         ! Update position
         p%pos = p%pos + p%vel*this%dt
         ! Put particle back in storage
         this%ps(i) = p
      end do
      
      ! Update time and step counter
      this%t = this%t + this%dt
      this%step = this%step + 1
   end subroutine increment_time

   subroutine check_and_correct_bounds(this)
      use mpi
      implicit none
      class(simulation), intent(inout) :: this
      real(8) :: r
      integer :: i, ierr

      ! Reset counter
      this%particles_removed=0

      do i=1,this%npart
         r = norm2(this%ps(i)%pos)
         if (r.gt.this%length) then
            ! this%ps(i)%pos = get_rand_pos(this%length)
            this%ps(i)%pos = this%ps(i)%pos - 2.0*this%ps(i)%pos
            this%particles_removed = this%particles_removed + 1
         end if
      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE, this%particles_removed, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
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

   write(*,'(A, F8.2, A, I5, A, I5)')  &
      & "Time :: ", this%t, "   Step :: ", this%step, "  Particles removed :: ", this%particles_removed
  end subroutine print


  subroutine compute_rdf(this)
   use mpi
   implicit none
   class(simulation), intent(inout) :: this
   real(8) :: r, density, vol
   integer :: i, ir, ierr
   real(8), dimension(:), allocatable :: tmp_rdf

   allocate(tmp_rdf(1:this%numbins)); tmp_rdf=0.0

   do i=1,this%npart
      r  = norm2(this%ps(i)%pos)
      ir = ceiling(r/this%delta)
      ir = min(this%numbins, ir)
      tmp_rdf(ir) = tmp_rdf(ir) + 1     
   end do
   
   density = this%npart/(4.0/3.0*PI*this%length**3)
   
   do i=1,this%numbins
      vol = 4.0/3.0*PI*((this%delta*i)**3 - (this%delta*(i-1))**3)
      tmp_rdf(i) = tmp_rdf(i) / density / vol
   end do

   call MPI_ALLREDUCE(tmp_rdf, this%rdf, this%numbins, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

   this%rdf=this%rdf/this%numproc
  end subroutine compute_rdf

  subroutine compute_w1(this)
   use mpi
   implicit none
   class(simulation), intent(inout) :: this
   real(8) :: r
   real(8), dimension(3) :: rhat
   integer :: i, ierr, ir
   
   ! Reset
   this%w1=0.0
   this%c=0

   ! Compute and bin relative velocity
   do i=1,this%npart
      r  = norm2(this%ps(i)%pos)
      ir = ceiling(r/this%delta)
      ir = min(this%numbins, ir)
      rhat = this%ps(i)%pos/r
      
      this%w1(ir) = this%w1(ir) + dot_product(this%ps(i)%vel,rhat) 
      this%c(ir) = this%c(ir) + 1
   end do

   ! Normalize by bin counter
   do i=1,this%numbins
      if (this%c(i).gt.0) then
         this%w1(i) = this%w1(i)/this%c(i)
      end if
   end do

   ! Reduce and normalize by number of processes
   call MPI_ALLREDUCE(MPI_IN_PLACE, this%w1, this%numbins, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   this%w1 = this%w1/this%numproc
  end subroutine compute_w1

  subroutine compute_w2(this)
   use mpi
   implicit none
   class(simulation), intent(inout) :: this
   real(8) :: r
   real(8), dimension(3) :: rhat
   integer :: i, ierr, ir

   ! Reset
   this%w2=0.0
   this%c=0

   ! Compute and bin relative velocity
   do i=1,this%npart
      r  = norm2(this%ps(i)%pos)
      ir = ceiling(r/this%delta)
      ir = min(this%numbins, ir)
      rhat = this%ps(i)%pos/r
      
      this%w2(ir) = this%w2(ir) + dot_product(this%ps(i)%vel,rhat)**2
      this%c(ir) = this%c(ir) + 1
   end do

   ! Normalize by bin counter
   do i=1,this%numbins
      if (this%c(i).gt.0) then
         this%w2(i) = this%w2(i)/this%c(i)
      end if
   end do

   ! Reduce and normalize by number of processes
   call MPI_ALLREDUCE(MPI_IN_PLACE, this%w2, this%numbins, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   this%w2 = this%w2/this%numproc
  end subroutine compute_w2

  subroutine write_rdf(this, filename)
   implicit none
   class(simulation), intent(in) :: this
   character(len=*), intent(in) :: filename
   integer :: i
   integer :: unit
   
   ! Open the binary file for writing
   open(newunit=unit, file=filename, status='replace', form='unformatted', access='stream')

   ! Write velocities to the binary file
   do i = 1, this%numbins
      write(unit) this%rdf(i)
   end do

   ! Close the binary file
   close(unit) 
  end subroutine write_rdf

  subroutine write_w1(this, filename)
   implicit none
   class(simulation), intent(in) :: this
   character(len=*), intent(in) :: filename
   integer :: i
   integer :: unit
   
   ! Open the binary file for writing
   open(newunit=unit, file=filename, status='replace', form='unformatted', access='stream')

   ! Write velocities to the binary file
   do i = 1, this%numbins
      write(unit) this%w1(i)
   end do

   ! Close the binary file
   close(unit) 
  end subroutine write_w1

  subroutine write_w2(this, filename)
   implicit none
   class(simulation), intent(in) :: this
   character(len=*), intent(in) :: filename
   integer :: i
   integer :: unit
   
   ! Open the binary file for writing
   open(newunit=unit, file=filename, status='replace', form='unformatted', access='stream')

   ! Write velocities to the binary file
   do i = 1, this%numbins
      write(unit) this%w2(i)
   end do

   ! Close the binary file
   close(unit) 
  end subroutine write_w2


end module simulation_class
