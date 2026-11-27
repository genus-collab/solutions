# Stage 5: two independent checks (main session, 2026-11-14)

New code, Python stdlib only, sharing nothing with the agents' code (no project, worker, referee or critic module is imported or
read). Both scripts print to stdout only; the logs here were made with `py -I -B script.py N > logs/script_N.log`. Finite checks are
evidence, not proof. Context: `../STAGE5_SYNTHESIS.md` and `../../MECHANISM_STAGE5.md` s.4.2.

| script | what it checks | logs |
|---|---|---|
| `cyc_check.py N` | Theorem CYC (`../critic/CYC_PROOF.md`): Frobenius basis of PF4 at `u = 0` from the recursion; the vanishing period found from the `log u` obstruction (`g = f_1 + (7/10) f_2`, nothing assumed); the factorised operator `(1/y0) D (1/Dt) D (u y0^2 Dt^2) D (1/Dt) D (1/y0)` equals PF4 on `u^0..u^5` and kills `f_0, f_1, f_2` through `u^N`, in four frames `y0 = f0 + a f1 + b f2` and with `g` rescaled. Controls: non-isotropic pairs `t = f1/y0`, `t = (f1 + f2/2)/y0`, and `g + f2/7` | `cyc_check_60.log` (11 s) |
| `lmi_check.py N` | the two-variable chart proof of (Gamma) (`../chart/referee_adversary/REFEREE.md` s.3, Theorem LMI with `s* = 1`, `alpha0 = -1`, `M0 = 0`), exact over Q (symmetric functions of `zeta^1..zeta^4` from `1 + x + .. + x^4`): the chart series equals the one-variable full-Stirling `c(z; L)` at `L = 2, -2/3, 1/2, 3`, and equals `h_(w_0)` at `L = 1`, with `h` computed from its definition by Newton sums, through `z^N`. Controls: the singular factor's Stirling series dropped; twist weight `+5`; `Q e^(2 mu)` for `Q e^mu` | `lmi_check_8.log`, `lmi_check_12.log`, `lmi_check_16.log` (201 s) |
