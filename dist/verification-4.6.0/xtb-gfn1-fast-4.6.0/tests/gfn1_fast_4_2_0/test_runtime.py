"""Actual --hess dispatch, numerical parity, fallback and stock timings."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import random
import re
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('grad', ROOT/'tests/gfn1_fast_2_4_0/test_runtime.py')
grad = importlib.util.module_from_spec(spec)
spec.loader.exec_module(grad)
WATER = '3\nwater\nO 0 0 0\nH .758602 .13 .504284\nH -.758602 -.18 .62\n'
ENABLE = {'XTB_GFN1_FAST_ENABLE_ANALYTIC_HESSIAN': '1'}
DISABLE = {'XTB_GFN1_FAST_DISABLE_ANALYTIC_HESSIAN': '1'}
NUMBER = r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'


def benzene():
    lines=['12','benzene']
    for element, radius in [('C',1.4),('H',2.48)]:
        for i in range(6):
            a=math.pi*i/3
            lines.append(f'{element} {radius*math.cos(a):.12f} {radius*math.sin(a):.12f} 0')
    return '\n'.join(lines)+'\n'


def run(exe,output,name,xyz,options=(),policy=None,step=.0005,control='',threads=1,profile=True,defaults=False,timeout=180):
    work=output/name
    work.mkdir(parents=True,exist_ok=True)
    (work/'input.xyz').write_text(xyz,encoding='ascii')
    hess_control='' if defaults else f'$hess\nsccacc=0.0000001\nstep={step}\n$end\n'
    (work/'control.inp').write_text(hess_control+control+'\n',encoding='ascii')
    env={k:v for k,v in os.environ.items() if not k.startswith('XTB_GFN1_FAST_')}
    runtime=[str(Path(env[k])/'bin') for k in ('CMPLR_ROOT','MKLROOT') if k in env]
    env['PATH']=os.pathsep.join([*runtime,env.get('PATH','')])
    env.update(OMP_NUM_THREADS=str(threads),MKL_NUM_THREADS='1',OMP_STACKSIZE='64M',KMP_STACKSIZE='64M',
               XTBPATH=str(ROOT),XTB_GFN1_FAST_PROFILE=str(int(profile)))
    # Historical parity tests request numerical explicitly through None;
    # an empty dict deliberately exercises the current production default.
    env.update(DISABLE if policy is None else policy)
    accuracy=[] if defaults else ['--acc','0.0000001']
    command=[str(exe),'input.xyz','--gfn','1','--hess','--norestart',*accuracy,'--input','control.inp',*options]
    started=time.perf_counter()
    process=subprocess.run(command,cwd=work,env=env,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=timeout)
    elapsed=time.perf_counter()-started
    log=process.stdout+'\n'+process.stderr
    (work/'run.log').write_text(log,encoding='utf-8')
    if process.returncode or not (work/'hessian').is_file():
        raise AssertionError(f'{name}: actual --hess failed, exit {process.returncode}: {work / "run.log"}')
    n3=3*int(xyz.splitlines()[0])
    matrix=[float(v.replace('D','E').replace('d','e')) for line in (work/'hessian').read_text().splitlines()
            if not line.startswith('$') for v in re.findall(NUMBER,line)]
    if len(matrix)!=n3*n3 or not all(math.isfinite(v) for v in matrix):
        raise AssertionError(f'{name}: invalid Cartesian Hessian file')
    match=re.search(r'GFN1-fast 4\.\d+\.\d+ analytic Hessian: used=([TF])',log)
    reason=re.search(r'dispatch reason = (.*)',log)
    result={'name':name,'process_wall_seconds':elapsed,'analytic_used': match[1]=='T' if match else None,
            'reason':reason[1].strip() if reason else None,'step_bohr':.005 if defaults else step,
            'n3':n3,'log':str(work/'run.log')}
    print(f'PASS --hess {name}: analytic={result["analytic_used"]}, {elapsed:.4f}s',flush=True)
    return result,matrix,log


def compare(reference,candidate,label,tolerance):
    if len(reference[1])!=len(candidate[1]):
        raise AssertionError('Hessian dimensions differ')
    error=max(abs(a-b) for a,b in zip(reference[1],candidate[1]))
    if error>tolerance:
        raise AssertionError(f'{label}: max Hessian difference {error} > {tolerance}')
    print(f'PASS Hessian parity {label}: {error:.3e} Eh/bohr^2',flush=True)
    return error


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--exe',type=Path,required=True)
    parser.add_argument('--reference',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--small-only',action='store_true')
    parser.add_argument('--repeats',type=int,default=7)
    args=parser.parse_args()
    exe,output=args.exe.resolve(),args.output.resolve()
    output.mkdir(parents=True,exist_ok=True)
    report={'exe':str(exe),'exe_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'cases':[],'benchmarks':[]}
    cases=[('water',WATER,()),('hcl','2\nHCl\nCl 0 0 0\nH 1.3 .1 .05\n',()),('water_zero',WATER,('--etemp','0')),
           ('water_highT',WATER,('--etemp','30000')),('water_open',WATER,('--chrg','1','--uhf','1'))]
    if not args.small_only:
        cases += [('disilane',grad.DISILANE,()),('benzene',benzene(),())]
    references={}
    for name,xyz,options in cases:
        numerical=run(exe,output,name+'_numerical',xyz,options)
        analytic=run(exe,output,name+'_analytic',xyz,options,ENABLE)
        if numerical[0]['analytic_used'] is not False or analytic[0]['analytic_used'] is not True:
            raise AssertionError(f'{name}: numerical/analytic dispatch was not exercised')
        error=compare(numerical,analytic,name,2e-6)
        references[name]=numerical
        report['cases'].append({'case':name,'max_abs_H_difference_Eh_bohr2':error})
    default_cases=[('water_defaults',WATER)]
    if not args.small_only:
        default_cases.append(('disilane_defaults',grad.DISILANE))
    for name,xyz in default_cases:
        numerical=run(exe,output,name+'_numerical',xyz,defaults=True)
        analytic=run(exe,output,name+'_analytic',xyz,policy=ENABLE,defaults=True)
        if not analytic[0]['analytic_used']:
            raise AssertionError(f'{name}: ordinary default-accuracy analytic request was rejected')
        error=compare(numerical,analytic,name,2e-4)
        report['cases'].append({'case':name,'step_bohr':.005,'max_abs_H_difference_Eh_bohr2':error})
    fallback=[('salt',WATER,('--gbsa','water'),ENABLE,'$gbsa\nion_st=0.1\n$end\n'),
              ('disable',WATER,(),{**ENABLE,'XTB_GFN1_FAST_DISABLE_ANALYTIC_HESSIAN':'1'},''),
              ('constraint',WATER,(),ENABLE,'$constrain\ndistance: 1,2,1.0\n$end\n')]
    if not args.small_only:
        fallback.append(('memory',grad.DISILANE,(),{**ENABLE,'XTB_GFN1_FAST_ANALYTIC_HESSIAN_MAX_MIB':'1'},''))
    for name,xyz,options,policy,control in fallback:
        candidate=run(exe,output,name+'_fallback',xyz,options,policy,control=control)
        reference=run(exe,output,name+'_reference',xyz,options,control=control)
        if candidate[0]['analytic_used'] is not False:
            raise AssertionError(f'{name}: unsupported request did not fall back')
        error=compare(reference,candidate,name+'_fallback',1e-10)
        report['cases'].append({'case':name+'_fallback','max_abs_H_difference_Eh_bohr2':error,
                                'dispatch_reason':candidate[0]['reason']})
    if not args.small_only:
        parallel=run(exe,output,'parallel_numerical',grad.DISILANE,threads=4)
        error=compare(references['disilane'],parallel,'OpenMP numerical fallback',1e-9)
        report['cases'].append({'case':'parallel_numerical','max_abs_H_difference_Eh_bohr2':error})
    if args.reference:
        if args.repeats<5:
            parser.error('At least five repeats for performance comparisons')
        stock=args.reference.resolve()
        report['stock_sha256']=hashlib.sha256(stock.read_bytes()).hexdigest()
        report['stock_hessian_compatibility']='Stock numerical formula; OpenMP directives suppressed for Windows ifx descriptor crash; one thread only'
        rng=random.Random(20261005)
        bench=[('disilane',grad.DISILANE),('benzene',benzene())]
        for molecule,xyz in bench:
            times={k:[] for k in ('stock','fast_numerical','fast_analytic')}
            for repeat in range(-1,args.repeats):
                order=list(times);rng.shuffle(order)
                trial={}
                for name in order:
                    trial[name]=run(stock if name=='stock' else exe,output,f'bench_{molecule}_{repeat}_{name}',xyz,
                                    policy=ENABLE if name=='fast_analytic' else DISABLE,step=.0005,profile=False)
                    if repeat>=0:
                        times[name].append(trial[name][0]['process_wall_seconds'])
                compare(trial['stock'],trial['fast_numerical'],'stock numerical '+molecule,1e-9)
                compare(trial['stock'],trial['fast_analytic'],'stock analytic '+molecule,2e-6)
            # Reconfirm the optimized path after timing (profiles are off above).
            diagnostics=run(exe,output,f'bench_{molecule}_profile',xyz,policy=ENABLE)
            if not diagnostics[0]['analytic_used']:
                raise AssertionError('Measured analytic policy silently fell back')
            medians={k:statistics.median(v) for k,v in times.items()}
            report['benchmarks'].append({'molecule':molecule,'profile':False,'warmups':1,'threads':1,
                                         'samples_seconds':times,'median_seconds':medians,
                                         'speedup_over_stock':medians['stock']/medians['fast_analytic'],
                                         'speedup_over_fast_numerical':medians['fast_numerical']/medians['fast_analytic']})
            print('BENCH Hessian',molecule,medians,flush=True)
    filename='benchmark_summary.json' if args.reference else 'summary.json'
    (output/filename).write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('PASS actual CLI Hessian dispatch/regression',flush=True)


if __name__=='__main__':
    main()
