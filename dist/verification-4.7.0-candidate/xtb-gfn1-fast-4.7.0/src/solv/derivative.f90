! This file is part of xtb.
! SPDX-License-Identifier: LGPL-3.0-or-later
! Small scalar derivative algebra for the three-variable Born kernels and
! six independent inertia entries. No molecular-size automatic arrays.
module xtb_solv_derivative
   use xtb_mctc_accuracy, only : wp
   implicit none
   private
   public :: TDerivative,derivativeVariable,derivativeExp,operator(+),operator(-),operator(*),operator(/),operator(**)
   type :: TDerivative
      real(wp) :: value=0.0_wp
      real(wp) :: gradient(6)=0.0_wp
      real(wp) :: hessian(6,6)=0.0_wp
   end type
   interface operator(+)
      module procedure add_dd,add_sd,add_ds
   end interface
   interface operator(-)
      module procedure sub_dd,sub_sd,sub_ds,neg_d
   end interface
   interface operator(*)
      module procedure mul_dd,mul_sd,mul_ds
   end interface
   interface operator(/)
      module procedure div_dd,div_sd,div_ds
   end interface
   interface operator(**)
      module procedure power_ds,power_di
   end interface
contains
pure function derivativeVariable(value,index) result(d)
   real(wp), intent(in) :: value
   integer, intent(in) :: index
   type(TDerivative) :: d
   if(index<1.or.index>6) error stop 'internal derivative variable index'
   d%value=value;d%gradient(index)=1.0_wp
end function
pure function chain(a,value,first,second) result(d)
   type(TDerivative), intent(in) :: a
   real(wp), intent(in) :: value,first,second
   type(TDerivative) :: d
   d%value=value;d%gradient=first*a%gradient
   d%hessian=first*a%hessian+second*spread(a%gradient,1,6)*spread(a%gradient,2,6)
end function
pure function derivativeExp(a) result(d)
   type(TDerivative), intent(in) :: a
   type(TDerivative) :: d
   real(wp) :: value
   value=exp(a%value);d=chain(a,value,value,value)
end function
pure function power_ds(a,p) result(d)
   type(TDerivative), intent(in) :: a
   real(wp), intent(in) :: p
   type(TDerivative) :: d
   if(p==0.0_wp) then
      d=TDerivative(1.0_wp)
   else if(p==1.0_wp) then
      d=a
   else
      d=chain(a,a%value**p,p*a%value**(p-1.0_wp),p*(p-1.0_wp)*a%value**(p-2.0_wp))
   endif
end function
pure function power_di(a,p) result(d)
   type(TDerivative), intent(in) :: a
   integer, intent(in) :: p
   type(TDerivative) :: d
   d=power_ds(a,real(p,wp))
end function
pure function add_dd(a,b) result(d)
   type(TDerivative), intent(in) :: a,b
   type(TDerivative) :: d
   d%value=a%value+b%value;d%gradient=a%gradient+b%gradient;d%hessian=a%hessian+b%hessian
end function
pure function add_sd(a,b) result(d)
   real(wp), intent(in) :: a
   type(TDerivative), intent(in) :: b
   type(TDerivative) :: d
   d=b;d%value=a+b%value
end function
pure function add_ds(a,b) result(d)
   type(TDerivative), intent(in) :: a
   real(wp), intent(in) :: b
   type(TDerivative) :: d
   d=add_sd(b,a)
end function
pure function sub_dd(a,b) result(d)
   type(TDerivative), intent(in) :: a,b
   type(TDerivative) :: d
   d%value=a%value-b%value;d%gradient=a%gradient-b%gradient;d%hessian=a%hessian-b%hessian
end function
pure function sub_sd(a,b) result(d)
   real(wp), intent(in) :: a
   type(TDerivative), intent(in) :: b
   type(TDerivative) :: d
   d=neg_d(b);d%value=a-b%value
end function
pure function sub_ds(a,b) result(d)
   type(TDerivative), intent(in) :: a
   real(wp), intent(in) :: b
   type(TDerivative) :: d
   d=a;d%value=a%value-b
end function
pure function neg_d(a) result(d)
   type(TDerivative), intent(in) :: a
   type(TDerivative) :: d
   d%value=-a%value;d%gradient=-a%gradient;d%hessian=-a%hessian
end function
pure function mul_dd(a,b) result(d)
   type(TDerivative), intent(in) :: a,b
   type(TDerivative) :: d
   d%value=a%value*b%value;d%gradient=a%gradient*b%value+a%value*b%gradient
   d%hessian=a%hessian*b%value+a%value*b%hessian &
      & +spread(a%gradient,1,6)*spread(b%gradient,2,6)+spread(b%gradient,1,6)*spread(a%gradient,2,6)
end function
pure function mul_sd(a,b) result(d)
   real(wp), intent(in) :: a
   type(TDerivative), intent(in) :: b
   type(TDerivative) :: d
   d%value=a*b%value;d%gradient=a*b%gradient;d%hessian=a*b%hessian
end function
pure function mul_ds(a,b) result(d)
   type(TDerivative), intent(in) :: a
   real(wp), intent(in) :: b
   type(TDerivative) :: d
   d=mul_sd(b,a)
end function
pure function div_dd(a,b) result(d)
   type(TDerivative), intent(in) :: a,b
   type(TDerivative) :: d
   d=a*(b**(-1.0_wp))
end function
pure function div_sd(a,b) result(d)
   real(wp), intent(in) :: a
   type(TDerivative), intent(in) :: b
   type(TDerivative) :: d
   d=a*(b**(-1.0_wp))
end function
pure function div_ds(a,b) result(d)
   type(TDerivative), intent(in) :: a
   real(wp), intent(in) :: b
   type(TDerivative) :: d
   d=(1.0_wp/b)*a
end function
end module
