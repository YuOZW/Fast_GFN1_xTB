"""Development probe; prints actual accepted counts, not just parity."""
import importlib.util
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('smoke',ROOT/'tests/gfn1_fast_2_2_6/test_windows_smoke.py')
smoke=importlib.util.module_from_spec(spec); spec.loader.exec_module(smoke)
spec=importlib.util.spec_from_file_location('gradient',ROOT/'tests/gfn1_fast_2_4_0/test_runtime.py')
gradient=importlib.util.module_from_spec(spec); spec.loader.exec_module(gradient)
exe=ROOT/'build-gfn1-fast-3.0.0-windows-ifx-release/xtb.exe'
output=exe.parent/'foe-probes'
force={'XTB_GFN1_FAST_ENABLE_FERMI_OPERATOR':'1',
       'XTB_GFN1_FAST_DISABLE_FERMI_OPERATOR_AUTOTUNE':'1',
       'XTB_GFN1_FAST_FERMI_OPERATOR_MIN_NAO':'1'}
for molecule,xyz in [('disilane',gradient.DISILANE),('taxol',(ROOT/'assets/inputs/xyz/taxol.xyz').read_text())]:
    for temperature in ('0','300','1000','5000'):
        name=molecule+'_'+temperature
        options=('--etemp',temperature)
        reference=smoke.run_case(exe,output,name+'_full',xyz,options)
        candidate=smoke.run_case(exe,output,name+'_foe',xyz,options,force)
        gradient.compare(reference,candidate,name)
        print(candidate[2].split('GFN1-fast 3.0.0 Fermi operator profile')[1].split('#    Occupation')[0].strip(),flush=True)
