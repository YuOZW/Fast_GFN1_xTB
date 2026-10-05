module xtb_mctc_accuracy
  implicit none
  integer, parameter :: wp=kind(1.0d0)
end module


module xtb_mctc_convert
  use xtb_mctc_accuracy, only: wp
  implicit none
  real(wp), parameter :: autoev=27.211386245988_wp
end module

module xtb_mctc_blas
  use xtb_mctc_accuracy, only: wp
  implicit none
contains
  subroutine mctc_gemm(a,b,c,transa,transb,alpha,beta)
    real(wp), intent(in) :: a(:,:), b(:,:)
    real(wp), intent(out) :: c(:,:)
    character(len=1), intent(in), optional :: transa,transb
    real(wp), intent(in), optional :: alpha,beta
    character(len=1) :: ta,tb
    integer :: m,n,k,lda,ldb,ldc
    real(wp) :: aa,bb
    ta='N'; tb='N'; aa=1.0_wp; bb=0.0_wp
    if(present(transa))ta=transa; if(present(transb))tb=transb
    if(present(alpha))aa=alpha; if(present(beta))bb=beta
    if(ta=='N'.or.ta=='n')then; m=size(a,1); k=size(a,2); lda=max(1,size(a,1))
    else; m=size(a,2); k=size(a,1); lda=max(1,size(a,1)); endif
    if(tb=='N'.or.tb=='n')then; n=size(b,2); ldb=max(1,size(b,1))
    else; n=size(b,1); ldb=max(1,size(b,1)); endif
    ldc=max(1,size(c,1))
    call dgemm(ta,tb,m,n,k,aa,a,lda,b,ldb,bb,c,ldc)
  end subroutine
end module

module xtb_mctc_lapack
  use xtb_mctc_accuracy, only: wp
  implicit none
  interface lapack_syevd
    module procedure wrap_syevd
  end interface
contains
  subroutine wrap_syevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,info)
    character(len=1),intent(in)::jobz,uplo
    integer,intent(in)::n,lda,lwork,liwork
    real(wp),intent(inout)::a(lda,*),work(*)
    real(wp),intent(out)::w(*)
    integer,intent(inout)::iwork(*)
    integer,intent(out)::info
    call dsyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,info)
  end subroutine
end module
