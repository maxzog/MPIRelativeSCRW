program test
   use simulation_class
   use mpi
   implicit none

   type(simulation) :: sim

   integer :: np, nbins, stepf, i, period
   real(8) :: L, t, tf, dt, delta
   real(8) :: a, b

   ! MPI Stuff
   integer :: ierr, rank, num_procs, fh, status(MPI_STATUS_SIZE)
   integer(kind=MPI_OFFSET_KIND) :: disp

   character(len=30) :: filename

   ! Init MPI and get details on comm
   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

   ! Sim parameters
   L = 3.5
   nbins = 128
   delta = L/nbins
   np = 1000000


   ! State compute and write period
   period = 250

   t = 0.0
   tf = 50.0
   dt = 0.001

   a = 15.379520019471148
   b = 8.3739639785560840

   call sim%init(npart=np, length=L, numbins=nbins, delta=delta, t = t, &
                & tf = tf, dt = dt, step = 0, stepf=999999, rank=rank, numproc=num_procs)
   call sim%write_particle_data("./outs/init.dat")
   call sim%compute_rdf()
   if (rank.eq.0) call sim%write_rdf("./outs/rdf_init.dat")

   ! Make sure we're all ready to rock
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   ! Integrate
   do while (sim%isRun)
      ! Check for out-of-bound particles and replace them
      call sim%check_and_correct_bounds()
      ! Update the equations
      call sim%increment_time(a, b)
      ! Check how we're doing on time
      call sim%check_time()
      ! Print time and step counter
      if (sim%rank.eq.0) call sim%print()
      if (mod(sim%step, period).eq.0) then
         ! Compute stats
         call sim%compute_rdf()
         call sim%compute_w1()
         call sim%compute_w2()
         ! Write rdf
         write(filename, '("./outs/rdf_", I0, ".dat")') sim%step
         if (sim%rank.eq.0) call sim%write_rdf(filename)
         ! Write first moment of velocity | r
         write(filename, '("./outs/w1_", I0, ".dat")') sim%step
         if (sim%rank.eq.0) call sim%write_w1(filename)
         ! Write second moment of velocity | r
         write(filename, '("./outs/w2_", I0, ".dat")') sim%step
         if (sim%rank.eq.0) call sim%write_w2(filename)
      end if
   end do

   ! Write final particle state to file
   write(filename, '("./outs/part_", I0, ".dat")') sim%step
   call sim%write_particle_data(filename)
   ! Compute stats
   call sim%compute_rdf()
   call sim%compute_w1()
   call sim%compute_w2()
   ! Write rdf
   write(filename, '("./outs/rdf_", I0, ".dat")') sim%step
   if (sim%rank.eq.0) call sim%write_rdf(filename)
   ! Write first moment of velocity | r
   write(filename, '("./outs/w1_", I0, ".dat")') sim%step
   if (sim%rank.eq.0) call sim%write_w1(filename)
   ! Write second moment of velocity | r
   write(filename, '("./outs/w2_", I0, ".dat")') sim%step
   if (sim%rank.eq.0) call sim%write_w2(filename)

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   call MPI_FINALIZE(ierr)
end program test
