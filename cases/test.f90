program test
   use simulation_class
   use mpi
   implicit none

   type(simulation) :: sim

   integer :: np, nbins, stepf, i
   real(8) :: L, t, tf, dt, delta
   real(8) :: a, b

   ! MPI Stuff
   integer :: ierr, rank, num_procs, fh, status(MPI_STATUS_SIZE)
   integer(kind=MPI_OFFSET_KIND) :: disp

   ! Init MPI and get details on comm
   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

   ! Sim parameters
   L = 6.2832
   nbins = 64
   delta = L/nbins
   np = 250000

   t = 0.0
   tf = 11.0
   dt = 0.001

   a = 15.379520019471148
   b = 8.3739639785560840

   call sim%init(npart=np, length=L, numbins=nbins, delta=delta, t = t, tf = tf, dt = dt, step = 0, stepf=9999, rank=rank)
   call sim%write_particle_data("./outs/init.dat")

   ! Integrate
   do while (sim%isRun)
      call sim%check_and_correct_bounds()
      call sim%increment_time(a, b)
      call sim%check_time()
      if (sim%rank.eq.0) call sim%print()
      if (mod(sim%step, 100).eq.0) then
         call sim%write_particle_data("./outs/test.dat")
      end if
   end do

   call MPI_FINALIZE(ierr)
end program test