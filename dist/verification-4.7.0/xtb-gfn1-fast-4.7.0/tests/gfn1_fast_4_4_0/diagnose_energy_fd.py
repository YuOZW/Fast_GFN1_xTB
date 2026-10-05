"""Record step convergence of Energy derivatives independently of Gradient."""
import argparse
import json
from pathlib import Path

from test_runtime import NODAL, smoke


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--coordinate', type=int, default=5)
    args = parser.parse_args()
    exe, output = args.exe.resolve(), args.output.resolve()
    options = ('--acc', '0.0001', '--gbsa', 'water')
    center = smoke.run_case(exe, output, 'center', NODAL, options)
    rows = []
    for step in (4e-4, 2e-4, 1e-4, 5e-5):
        energies = []
        for sign in (-1, 1):
            lines = NODAL.splitlines()
            index = args.coordinate // 3 + 2
            fields = lines[index].split()
            axis = args.coordinate % 3 + 1
            fields[axis] = str(float(fields[axis]) + sign * step * .529177210903)
            lines[index] = ' '.join(fields)
            result = smoke.run_case(exe, output, f'h{step}_{sign}',
                                    '\n'.join(lines) + '\n', options)
            energies.append(result[0]['energy_Eh'])
        fd = (energies[1] - energies[0]) / (2 * step)
        row = {'step_bohr': step, 'energy_derivative': fd,
               'gradient': center[1][args.coordinate],
               'difference': fd - center[1][args.coordinate]}
        rows.append(row)
        print('ENERGY FD', row, flush=True)
    (output / 'step_convergence.json').write_text(
        json.dumps({'coordinate': args.coordinate, 'rows': rows}, indent=2),
        encoding='utf-8')


if __name__ == '__main__':
    main()
