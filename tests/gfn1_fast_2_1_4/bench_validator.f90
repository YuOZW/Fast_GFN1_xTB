program b
 use xtb_mctc_accuracy,only:wp
 use xtb_gfn1_davidson,only:validateGFN1DavidsonEigensystem
 implicit none
 integer,parameter::n=1024,k=96,np=n,reps=20
 real(wp),allocatable::ph(:),ps(:),c(:,:),e(:),hc(:,:),sc(:,:)
 integer,allocatable::ml(:,:)
 integer::i,r,rate,t0,t1
 real(wp)::mr,mn,mo,mp,dt
 logical::ok
 allocate(ph(np),ps(np),c(n,k),e(k),hc(n,k),sc(n,k),ml(2,np))
 c=0
 do i=1,n
   ph(i)=-5d0+0.01d0*i;ps(i)=1d0;ml(:,i)=i
 enddo
 do i=1,k;c(i,i)=1d0;e(i)=ph(i);enddo
 call system_clock(count_rate=rate);call system_clock(t0)
 do r=1,reps
  call validateGFN1DavidsonEigensystem(ph,ps,np,ml,c,e,hc,sc,mr,mn,mo,mp,ok)
 enddo
 call system_clock(t1)
 dt=real(t1-t0,wp)/rate/reps
 print '(a,f10.6,a,l1)','strict validator seconds/call = ',dt,' ok=',ok
end program
