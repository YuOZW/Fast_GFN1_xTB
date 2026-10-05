! Concurrent initial policy reads, bounded teams and nested/no-OpenMP operation.
program test_resources
   use xtb_gfn1_fast_policy, only : TGFN1FastPolicy,getGFN1FastPolicy,getGFN1HessianResources
   use iso_fortran_env, only : int64
   implicit none
   type(TGFN1FastPolicy) :: policies(4),policy
   integer :: workers(4),serialWorkers,expected,i
   integer(int64) :: sizes(4),estimated,cap
   character(len=16) :: argument
   logical :: openmp
   call get_command_argument(1,argument);read(argument,*) expected
   openmp=.false.
   !$ openmp=.true.
   ! No prior initialization: all workers request the immutable policy.
   !$omp parallel do num_threads(4) default(none) shared(policies,workers,sizes)
   do i=1,4
      call getGFN1FastPolicy(policies(i))
      call getGFN1HessianResources(113,350,200,.true.,workers(i),sizes(i))
   enddo
   !$omp end parallel do
   if(openmp.and.any(workers/=1)) error stop 'nested response team requested'
   do i=2,4
      if(policies(i)%analyticHessianMaxMiB/=policies(1)%analyticHessianMaxMiB) error stop 'policy publication'
      if(sizes(i)/=sizes(1)) error stop 'concurrent resource estimates'
   enddo
   call getGFN1FastPolicy(policy)
   call getGFN1HessianResources(113,350,200,.true.,serialWorkers,estimated)
   if(serialWorkers/=expected) error stop 'requested/budgeted response workers'
   cap=int(policy%analyticHessianMaxMiB,int64)*1048576_int64
   if(estimated>cap.and.serialWorkers/=1) error stop 'over-budget parallel workspace'
   call getGFN1HessianResources(3,8,6,.false.,serialWorkers,estimated)
   if(serialWorkers/=1) error stop 'small-system serial crossover'
   print *, 'PASS concurrent policy, resource budget, nested and small-system teams; workers=',expected
end program
