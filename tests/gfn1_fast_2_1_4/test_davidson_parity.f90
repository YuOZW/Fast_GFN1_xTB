program test_214
  use xtb_mctc_accuracy, only: wp
  use xtb_gfn1_davidson, only: gfn1BlockDavidsonMF, validateGFN1DavidsonEigensystem
  implicit none
  integer, parameter :: n=160,nr=32
  real(wp) :: S(n,n),H0(n,n),H(n,n),Sfull(n,n),Hfull(n,n)
  real(wp) :: eval0(n),evalf(n),evald(nr),C0(n,n),Cf(n,n),Cd(n,nr)
  real(wp),allocatable :: pairS(:),pairH(:),hc(:,:),sc(:,:),carry(:,:)
  integer,allocatable :: ml(:,:)
  integer :: np,i,j,k,info,lwork,liwork,it,ncarry,nrestart,caseid
  real(wp),allocatable :: work(:)
  integer,allocatable :: iwork(:)
  logical :: ok,vok
  real(wp) :: maxr,maxn,maxo,maxp,ediff,goodr,goodn,goodo,goodp,maxediff

  np=n*(n+1)/2; allocate(pairS(np),pairH(np),ml(2,np),hc(n,nr),sc(n,nr),carry(n,nr+16))
  S=0;H0=0
  do i=1,n
    S(i,i)=1.0_wp + 0.02_wp*sin(real(i,wp))
    H0(i,i)=-4.0_wp + 0.035_wp*real(i-1,wp)
    do j=max(1,i-5),i-1
      S(i,j)=0.012_wp*cos(real(i+2*j,wp)); S(j,i)=S(i,j)
      H0(i,j)=0.025_wp*sin(real(2*i+j,wp)); H0(j,i)=H0(i,j)
    enddo
  enddo
  C0=H0; Sfull=S; call fullsolve(C0,Sfull,eval0,info); if(info/=0)stop 1
  ! C0 now eigenvectors
  maxediff=0.0_wp;goodr=0;goodn=0;goodo=0;goodp=0
  do caseid=1,12
    H=H0
    do i=1,n
      H(i,i)=H(i,i)+1.0e-3_wp*real(caseid,wp)*sin(0.17_wp*real(i,wp))
      do j=max(1,i-4),i-1
        H(i,j)=H(i,j)+2.0e-4_wp*real(caseid,wp)*cos(real(i+j,wp))
        H(j,i)=H(i,j)
      enddo
    enddo
    k=0
    do i=1,n
      do j=1,i
        k=k+1; ml(1,k)=i;ml(2,k)=j;pairS(k)=S(i,j);pairH(k)=H(i,j)
      enddo
    enddo
    Hfull=H; Sfull=S; call fullsolve(Hfull,Sfull,evalf,info); if(info/=0)stop 2
    call gfn1BlockDavidsonMF(pairH,pairS,S,np,ml,C0(:,1:nr+8),nr,evald,Cd,ok,it,maxr,carry,ncarry,nrestart)
    if(.not.ok)then; write(*,*)'Davidson failed case',caseid;stop 3;endif
    call validateGFN1DavidsonEigensystem(pairH,pairS,np,ml,Cd,evald,hc,sc,maxr,maxn,maxo,maxp,vok)
    if(.not.vok)then; write(*,*)'validator failed case',caseid,maxr,maxn,maxo,maxp;stop 4;endif
    ediff=maxval(abs(evald-evalf(1:nr)));maxediff=max(maxediff,ediff)
    goodr=max(goodr,maxr);goodn=max(goodn,maxn);goodo=max(goodo,maxo);goodp=max(goodp,maxp)
    if(ediff>1d-10)then; write(*,*)'eig mismatch',caseid,ediff;stop 5;endif
    C0=Hfull; eval0=evalf
  enddo

  ! Non-adjacent duplicate: old adjacent-only orthogonality check could miss C1=C3.
  Cd=C0(:,1:nr); evald=eval0(1:nr); Cd(:,3)=Cd(:,1)
  call validateGFN1DavidsonEigensystem(pairH,pairS,np,ml,Cd,evald,hc,sc,maxr,maxn,maxo,maxp,vok)
  if(vok)then; write(*,*)'validator accepted non-adjacent duplicate';stop 6;endif
  if(maxo<0.9_wp)then; write(*,*)'orthogonality defect not seen',maxo;stop 7;endif

  write(*,'(a,es12.4)') '12-case full parity eig max = ',maxediff
  write(*,'(a,es12.4)') 'valid max relative residual = ',goodr
  write(*,'(a,es12.4)') 'valid max S-norm error = ',goodn
  write(*,'(a,es12.4)') 'valid max S-orth error = ',goodo
  write(*,'(a,es12.4)') 'valid max projected-H error = ',goodp
  write(*,'(a,es12.4)') 'duplicate full-Gram orth error = ',maxo
  write(*,*) 'PASS'
contains
  subroutine fullsolve(A,B,w,info)
    real(wp),intent(inout)::A(n,n),B(n,n)
    real(wp),intent(out)::w(n)
    integer,intent(out)::info
    real(wp),allocatable::ww(:)
    integer,allocatable::ii(:)
    integer::lw,li
    lw=1+6*n+2*n*n;li=3+5*n
    allocate(ww(lw),ii(li))
    call dsygvd(1,'V','U',n,A,n,B,n,w,ww,lw,ii,li,info)
    deallocate(ww,ii)
  end subroutine
end program
