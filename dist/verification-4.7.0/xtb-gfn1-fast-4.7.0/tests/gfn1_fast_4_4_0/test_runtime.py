"""Nodal-safe GFN1 derivatives, legacy diagnostic fallback, and Gradient timings."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import re
import statistics

ROOT=Path(__file__).resolve().parents[2]


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


hess=load('hess','tests/gfn1_fast_4_2_0/test_runtime.py')
smoke=load('smoke','tests/gfn1_fast_2_2_6/test_windows_smoke.py')
grad=load('grad','tests/gfn1_fast_2_4_0/test_runtime.py')
NODAL='3\nnodal asymmetric water\nO 0 0 0\nH .758602 0 .504284\nH -.3 .65 .62\n'
LEGACY={'XTB_GFN1_FAST_DISABLE_DIRECT_H0_GRADIENT':'1'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--exe',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--reference',type=Path)
    args=parser.parse_args()
    exe,output=args.exe.resolve(),args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    report={'exe':str(exe),'exe_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'cases':[],'benchmarks':[]}
    numerical=hess.run(exe,output,'nodal_numerical',NODAL)
    analytic=hess.run(exe,output,'nodal_analytic',NODAL,policy=hess.ENABLE)
    if not analytic[0]['analytic_used']:
        raise AssertionError('Nodal analytic Hessian rejected')
    error=hess.compare(numerical,analytic,'nodal full CLI Hessian',2e-6)
    report['cases'].append({'case':'nodal_hessian','max_abs_H_difference_Eh_bohr2':error})
    for name,policy in [('legacy',{**LEGACY,**hess.ENABLE}),
                        ('generic',{**hess.ENABLE,'XTB_GFN1_FAST_DISABLE_GRADIENT_KERNEL':'1'})]:
        candidate=hess.run(exe,output,name+'_fallback',NODAL,policy=policy)
        ref=hess.run(exe,output,name+'_numerical',NODAL,policy={k:v for k,v in policy.items() if k not in hess.ENABLE})
        if candidate[0]['analytic_used'] is not False:
            raise AssertionError('Legacy Gradient requested with incompatible analytic response')
        error=hess.compare(ref,candidate,name+' derivative-policy fallback',1e-10)
        report['cases'].append({'case':name+'_fallback','max_abs_H_difference_Eh_bohr2':error})
    for name,options in [('gas',()),('alpb',('--alpb','water')),('gbsa',('--gbsa','water'))]:
        # The CLI clamps --acc below 1e-4. Strict 1e-7 SCC accuracy is tested
        # through the native API. Keep the actual CLI condition explicit.
        opts=('--acc','0.0001',*options)
        # This GBSA geometry retains the exact O_py/H1_s overlap node while
        # moving away from a 1e-6 SASA quadrature screening boundary. The
        # original geometry's discontinuity is independently reproduced by
        # probe_solvent_gradient and diagnose_energy_fd, including stock xtb.
        xyz=NODAL.replace('-.3 .65 .62','-.3 .67 .62') if name=='gbsa' else NODAL
        center=smoke.run_case(exe,output,'fd_'+name+'_center',xyz,opts)
        errors=[];step=2e-4*.529177210903
        for coordinate in range(9):
            energies=[]
            for sign in (-1,1):
                lines=xyz.splitlines();index=coordinate//3+2;fields=lines[index].split()
                fields[coordinate%3+1]=str(float(fields[coordinate%3+1])+sign*step)
                lines[index]=' '.join(fields)
                result=smoke.run_case(exe,output,f'fd_{name}_{coordinate}_{sign}','\n'.join(lines)+'\n',opts)
                energies.append(result[0]['energy_Eh'])
            fd=(energies[1]-energies[0])/(4e-4)
            errors.append(abs(fd-center[1][coordinate]))
        if max(errors)>1e-7:
            raise AssertionError(f'Nodal {name} Gradient differs from Energy derivatives: {errors}')
        report['cases'].append({'case':'nodal_'+name+'_energy_fd','accuracy':1e-4,
                               'xyz':xyz,'max_abs_G_difference_Eh_bohr':max(errors)})
        print('PASS all-coordinate nodal Energy derivatives',name,max(errors),flush=True)
    if args.reference:
        stock=args.reference.resolve();report['stock_sha256']=hashlib.sha256(stock.read_bytes()).hexdigest()
        taxol=(ROOT/'assets/inputs/xyz/taxol.xyz').read_text()
        for name,xyz,options in [('disilane',grad.DISILANE,()),('taxol_gas',taxol,()),
                                  ('taxol_alpb',taxol,('--alpb','water'))]:
            times={k:[] for k in ('stock','legacy_fast','direct_fast')};trials=[]
            rng=random.Random(20261005)
            for repeat in range(-1,7):
                order=list(times);rng.shuffle(order);trial={}
                for mode in order:
                    policy=LEGACY if mode=='legacy_fast' else {}
                    trial[mode]=smoke.run_case(stock if mode=='stock' else exe,output,
                        f'bench_{name}_{repeat}_{mode}',xyz,options,{**policy,'XTB_GFN1_FAST_PROFILE':'0'})
                    if repeat>=0:
                        times[mode].append(trial[mode][0]['wall_seconds'])
                for mode in ('legacy_fast','direct_fast'):
                    grad.compare(trial['stock'],trial[mode],name+' '+mode)
                if repeat>=0:
                    trials.append({'order':order,'results':{k:v[0] for k,v in trial.items()}})
            audit=smoke.run_case(exe,output,name+'_profile_audit',xyz,options)
            if 'direct H0 prefactor=T' not in audit[2]:
                raise AssertionError('Timed direct prefactor was not used')
            medians={k:statistics.median(v) for k,v in times.items()}
            report['benchmarks'].append({'case':name,'threads':1,'profile':False,'warmups':1,
                'samples_seconds':times,'median_seconds':medians,'trials':trials,
                'speedup_over_stock':medians['stock']/medians['direct_fast'],
                'speedup_over_legacy_fast':medians['legacy_fast']/medians['direct_fast']})
            print('BENCH direct H0 Gradient',name,medians,flush=True)
    filename='benchmark_summary.json' if args.reference else 'summary.json'
    (output/filename).write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('PASS nodal CLI derivative consistency',flush=True)


if __name__=='__main__':
    main()
