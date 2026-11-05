import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pathlib import Path
from NCq3solve.families import family, symmetries_for
from NCq3solve.pipeline import Solver

name = sys.argv[1] if len(sys.argv) > 1 else 'KP_real'
htilde = sys.argv[2] if len(sys.argv) > 2 else None
ans, slice_eqs = family('KP_real')
s = Solver(ans, htilde=htilde, symmetries=symmetries_for(name), extra_equations=slice_eqs)
s.generate()
# print(s.check_symmetries(do_log=True))


OUT_DIR = Path(__file__).resolve().parent / 'results'
out_json = OUT_DIR / f'{name}_{htilde}.json'
res = s.solve_all(out_json=str(out_json), timeout=1800)