program p
 use xtb_gfn1_fast_policy
 type(TGFN1FastPolicy)::x
 call getGFN1FastPolicy(x)
 if(.not.x%davidsonAudit) stop 1
 print *,'PASS audit policy'
end program
