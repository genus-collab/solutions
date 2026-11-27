# -*- coding: utf-8 -*-
"""rc_rerun_compare.py STORED_ROOT RERUN_ROOT -- critic of the Stage-5 referees: compare the referees' stored logs with the logs of
re-runs from scratch copies.  Timings (decimal seconds, 'd.d s') are replaced by 'T s' on both sides; the header line '# logs/...'
that cyc_formal prepended to its stored stdout captures is dropped.  Prints relative names only; writes logs/rc_rerun_compare.log.
STORED_ROOT = the folder allgenus/stage5/referees; RERUN_ROOT = the scratch folder holding the copies (layout of the critic's run)."""
import os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
STORED, RERUN = sys.argv[1], sys.argv[2]
PAIRS = [('cyc_formal', n, 'cyc_formal/logs/%s.log' % n) for n in
         ('cyc_ref_alg', 'cyc_ref_sol_80', 'cyc_ref_frame_80', 'cyc_ref_sym_30', 'cyc_ref_leib_40', 'cyc_ref_xproj')]
PAIRS += [('cyc_adversary', n, 'HVQm/a/b/c/cyc_adversary/logs/%s.log' % n) for n in
          ('cy_run_64', 'cy_frames_80', 'cy_rows_64_12', 'cy_rows_30_6', 'cy_proj_40', 'cy_extra_40')]
PAIRS += [('pb_steps', n, 'pb_steps/logs/%s.log' % n) for n in ('rps_run_v3_80_20_16_60', 'rps_run_v3_120_28_24_80')]
PAIRS += [('lmi_formal', n, 'lmi_formal/logs/%s.log' % n) for n in
          ('lf_dlmf_60_80', 'lf_quintic_14_60', 'lf_lmi_12_12', 'lf_lmi2_14', 'lf_fit_30')]
PAIRS += [('lmi_adversary', n, 'lmi_adversary/logs/%s.log' % n) for n in
          ('lm_run_eps_10_Q', 'lm_run_eps_12_p1', 'lm_run_eps_14_p2', 'lm_run_q_12_Q', 'lm_run_q_16_p1', 'lm_run_q_16_p2',
           'lm_run_x_12')]
TIME = re.compile(r'\d+\.\d+\s?s\b')


def norm(path):
    with open(path, encoding='utf-8') as fh:
        lines = [ln.rstrip() for ln in fh if not ln.startswith('# logs/')]
    return [TIME.sub('T s', ln) for ln in lines]


out, nbad = [], 0
for ref, name, rel in PAIRS:
    a = norm(os.path.join(STORED, ref, 'logs', name + '.log'))
    b = norm(os.path.join(RERUN, rel))
    verdict = 'ALL CHECKS PASSED' in ' '.join(b)
    same = a == b
    nbad += (not same) or (not verdict)
    diff = [i for i in range(max(len(a), len(b))) if (a[i] if i < len(a) else None) != (b[i] if i < len(b) else None)]
    out.append('%-14s %-28s %4d lines  %s%s' % (ref, name + '.log', len(a), 'IDENTICAL (timings stripped)' if same else
                                                'DIFFERS at lines %s' % diff[:6], '' if verdict else '  [re-run lacks ALL CHECKS PASSED]'))
out.append('%d logs compared: %s' % (len(PAIRS), 'all identical after stripping timings, all re-runs ALL CHECKS PASSED'
                                     if nbad == 0 else '%d PROBLEMS' % nbad))
print('\n'.join(out))
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
with open(os.path.join(HERE, 'logs', 'rc_rerun_compare.log'), 'w', encoding='utf-8') as fh:
    fh.write('# logs/rc_rerun_compare.log (py -I -B rc_rerun_compare.py <referees folder> <scratch re-run folder>)\n')
    fh.write('\n'.join(out) + '\n')
