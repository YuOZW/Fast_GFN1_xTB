"""Packed H0/direct-contraction parity, finite differences and timings."""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import statistics

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('smoke', ROOT / 'tests/gfn1_fast_2_2_6/test_windows_smoke.py')
smoke = importlib.util.module_from_spec(spec)
spec.loader.exec_module(smoke)
DISILANE = '''8
asymmetric disilane, s/p/d coverage
Si -1.15 0 0
Si 1.15 0.1 0.05
H -1.65 1.15 0.1
H -1.65 -0.6 1.0
H -1.7 -0.5 -1.05
H 1.7 1.2 0.15
H 1.65 -0.5 1.1
H 1.7 -0.4 -1.0
'''
LEGACY = {'XTB_GFN1_FAST_DISABLE_PACKED_GRADIENT': '1',
          'XTB_GFN1_FAST_DISABLE_DIRECT_H0_GRADIENT': '1',
          'XTB_GFN1_FAST_DISABLE_GRADIENT_CONTRACTION': '1'}


def compare(reference, candidate, label):
    de = abs(reference[0]['energy_Eh']-candidate[0]['energy_Eh'])
    dg = max(abs(a-b) for a,b in zip(reference[1],candidate[1]))
    if de>1e-11 or dg>1e-10:
        raise AssertionError(f'{label}: dE={de}, dG={dg}')
    print(f'PASS gradient parity {label}: dE={de:.3e}, max dG={dg:.3e}', flush=True)
    return {'label': label, 'abs_dE_Eh': de, 'max_dG_Eh_bohr': dg}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--exe',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--repeats',type=int,default=5)
    parser.add_argument('--small-only',action='store_true')
    args=parser.parse_args()
    exe,output=args.exe.resolve(),args.output.resolve()
    summary={'comparisons': [], 'benchmarks': []}
    taxol=(ROOT/'assets/inputs/xyz/taxol.xyz').read_text()
    for name,xyz in [('disilane',DISILANE)] + ([] if args.small_only else [('taxol',taxol)]):
        reference=smoke.run_case(exe,output,name+'_legacy',xyz,overrides=LEGACY)
        for mode,policy in [('both',{}),('packed_only',{'XTB_GFN1_FAST_DISABLE_GRADIENT_CONTRACTION':'1'}),
                            ('contract_only',{'XTB_GFN1_FAST_DISABLE_PACKED_GRADIENT':'1'}),
                            ('generic',{'XTB_GFN1_FAST_DISABLE_GRADIENT_KERNEL':'1'}),
                            ('threads4',{'OMP_NUM_THREADS':'4'})]:
            result=smoke.run_case(exe,output,name+'_'+mode,xyz,overrides=policy)
            summary['comparisons'].append(compare(reference,result,name+'_'+mode))
            if mode=='both' and 'packed H0=T, primitive contraction=T' not in result[2]:
                raise AssertionError('production gradient improvements not exercised')
        for suffix,options in [('alpb',('--alpb','water')),('open',('--chrg','1','--uhf','1'))]:
            ref=smoke.run_case(exe,output,name+'_'+suffix+'_legacy',xyz,options,LEGACY)
            new=smoke.run_case(exe,output,name+'_'+suffix+'_new',xyz,options)
            summary['comparisons'].append(compare(ref,new,name+'_'+suffix))

    # The full Energy at displaced coordinates checks the electronic,
    # repulsion, dispersion, CN and Coulomb gradient together.
    options=('--acc','0.1')
    equilibrium=smoke.run_case(exe,output,'finite_difference_center',DISILANE,options)
    h=1e-4
    for atom,axis in [(0,0),(2,1)]:
        energies=[]
        for sign in (-1,1):
            lines=DISILANE.splitlines(); record=lines[atom+2].split()
            record[axis+1]=str(float(record[axis+1])+sign*h)
            lines[atom+2]=' '.join(record)
            result=smoke.run_case(exe,output,f'fd_{atom}_{axis}_{sign}', '\n'.join(lines)+'\n',options)
            energies.append(result[0]['energy_Eh'])
        fd=(energies[1]-energies[0])/(2*h)*0.529177210903
        actual=equilibrium[1][3*atom+axis]
        if abs(fd-actual)>2e-5:
            raise AssertionError(f'full energy derivative atom {atom} axis {axis}: fd={fd}, gradient={actual}')
        print(f'PASS full Energy finite difference: delta={abs(fd-actual):.3e}',flush=True)
    if not args.small_only:
        for profile in ('0','1'):
            samples={'legacy':[],'new':[]}; kernel={'legacy':[],'new':[]}
            for repeat in range(args.repeats):
                for name,policy in [('legacy',LEGACY),('new',{})]:
                    result,_,log=smoke.run_case(exe,output,f'bench_{profile}_{repeat}_{name}',taxol,
                        overrides={**policy,'XTB_GFN1_FAST_PROFILE':profile})
                    samples[name].append(result['wall_seconds'])
                    match=re.search(r'GFN1 overlap kernel:\s*([0-9.]+)',log)
                    if match: kernel[name].append(float(match.group(1)))
            medians={name:statistics.median(values) for name,values in samples.items()}
            summary['benchmarks'].append({'profile':profile,'samples_seconds':samples,
                'median_seconds':medians,'new_over_legacy':medians['new']/medians['legacy'],
                'median_kernel_seconds':{name:statistics.median(values) for name,values in kernel.items() if values}})
            print(f'BENCH profile={profile}: {medians}',flush=True)
    output.mkdir(parents=True,exist_ok=True)
    (output/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')


if __name__=='__main__':
    main()
