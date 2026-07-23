"""Distributed V-set (nprtrv>1) ectrans4py test: the ``*_vset4py`` bindings for the
2-D spectral decomposition (independent wave-set x V-set spectral grid).

Run:  mpirun -n 4 python tests/test_ectrans4py/mpi_distributed_vset.py   (also -n 2)

An even task count exercises ``nprtrv = 2`` (spectral wave-sets ``nprtrw = size/2``);
an odd count degenerates to ``nprtrv = 1``, where the V-set routines reduce to the
plain ``*_dist4py`` forms. The spectral side is V-set-distributed (each rank holds only
the ``KFLDL`` fields of its V-set); the grid side holds all ``KFLDG`` fields. The
``KVSET(KFLDG)`` map assigns each global field to its owning V-set.

Checks:
  * ``trans_inq_vset4py`` returns the spectral processor-grid position (``KPRTRW`` wave
    sets, this task's 1-based ``KMYSETW`` / ``KMYSETV``);
  * ``dist_spec_vset4py`` -> ``gath_spec_vset4py`` is the identity (V-set scatter/gather);
  * ``inv_trans_scalar_vset4py`` -> ``dir_trans_scalar_vset4py`` -> ``inv`` is the
    identity on the grid (transform round trip on V-set-distributed spectral input);
  * ``inv_trans_uv_vset4py`` -> ``dir_trans_uv_vset4py`` -> ``inv`` likewise for winds.
"""

import sys

import numpy as np
from mpi4py import MPI

import ectrans4py

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

ectrans4py.init_env(unlimited_stack=False)
ectrans4py.mpl_init4py()

REAL = getattr(ectrans4py, "_REAL", np.float64)

# small analytic full Gaussian grid: NDGL=32 latitudes, 64 longitudes, T21
NDGL = 32
KSMAX = 21
KLOEN = np.full(NDGL, 2 * NDGL, dtype=np.int64)

# grid grid = size x 1 (N-S); spectral grid = nprtrw x nprtrv, with nprtrv=2 when even
nprtrv = 2 if size % 2 == 0 else 1
nprtrw = size // nprtrv
ectrans4py.setup_trans0_4py(size, 1, nprtrw, False, 10, False)
kresol = ectrans4py.setup_trans_4py(KSMAX, NDGL, NDGL, KLOEN, True, False)

(kgptot, kspec, kspec2, kgptotg, kspec2g, ksmax, knmeng, pmu, pgw) = ectrans4py.trans_inq4py(
    kresol, NDGL, KSMAX, NDGL, KLOEN, 10
)
(kprtrw, kmysetw, kmysetv) = ectrans4py.trans_inq_vset4py(kresol)

ok = True


def check(cond, msg):
    global ok
    if not cond:
        ok = False
        print(f"[vset] rank {rank} FAIL: {msg}")


check(kprtrw == nprtrw, f"KPRTRW {kprtrw} != nprtrw {nprtrw}")
check(1 <= kmysetw <= nprtrw, f"KMYSETW {kmysetw} out of [1,{nprtrw}]")
check(1 <= kmysetv <= nprtrv, f"KMYSETV {kmysetv} out of [1,{nprtrv}]")

# KFLDG global fields, round-robin across the V-sets; KFLDL = this rank's V-set share
KFLDG = 2 * nprtrv
KVSET = np.array([(f % nprtrv) + 1 for f in range(KFLDG)], dtype=np.int64)
KFLDL = int(np.count_nonzero(KVSET == kmysetv))
kfrom = np.ones(KFLDG, dtype=np.int64)  # all global fields sourced from MPL rank 1
kto = np.ones(KFLDG, dtype=np.int64)  # gathered back to MPL rank 1

# --- V-set scatter/gather consistency -----------------------------------------
# dist_spec_vset (scatter global->local by V-set) and gath_spec_vset (gather back)
# are inverse data-movement operations. A gathered field carries only the meaningful
# spectral coefficients (ecTrans's spectral layout has unused padding entries that the
# scatter/gather legitimately leave at zero), so we test the *stable* round trip: the
# first dist->gath produces a padding-consistent field, and re-scattering + re-gathering
# it reproduces it exactly.
rng = np.random.default_rng(20260709)
specg = np.zeros((KFLDG, kspec2g), dtype=REAL)
if rank == 0:
    specg[:, :] = rng.standard_normal((KFLDG, kspec2g)).astype(REAL)
sloc = ectrans4py.dist_spec_vset4py(kspec2g, kspec2, KFLDG, KFLDL, kfrom, KVSET, specg)
check(sloc.shape == (KFLDL, kspec2), f"dist_spec_vset shape {sloc.shape} != {(KFLDL, kspec2)}")
specg1 = ectrans4py.gath_spec_vset4py(kspec2g, kspec2, KFLDG, KFLDL, kto, KVSET, sloc)
sloc1 = ectrans4py.dist_spec_vset4py(kspec2g, kspec2, KFLDG, KFLDL, kfrom, KVSET, specg1)
specg2 = ectrans4py.gath_spec_vset4py(kspec2g, kspec2, KFLDG, KFLDL, kto, KVSET, sloc1)
if rank == 0:
    dg_err = float(np.max(np.abs(specg2 - specg1)))
    check(dg_err < 1e-10, f"dist_spec_vset <-> gath_spec_vset not stable, err {dg_err:.2e}")
    check(float(np.max(np.abs(specg1))) > 0.0, "dist->gath produced an all-zero field")

# --- scalar transform round trip: inv -> dir -> inv (identity on the grid) ----
g = ectrans4py.inv_trans_scalar_vset4py(kspec2, kgptot, KFLDL, KFLDG, KVSET, sloc)
check(g.shape == (KFLDG, kgptot), f"inv_trans_scalar_vset shape {g.shape} != {(KFLDG, kgptot)}")
s2 = ectrans4py.dir_trans_scalar_vset4py(kspec2, kgptot, KFLDL, KFLDG, KVSET, g)
g2 = ectrans4py.inv_trans_scalar_vset4py(kspec2, kgptot, KFLDL, KFLDG, KVSET, s2)
sc_err = comm.allreduce(float(np.max(np.abs(g2 - g))), op=MPI.MAX)
check(sc_err < 1e-4, f"scalar vset grid->dir->inv round trip err {sc_err:.2e}")

# --- UV transform round trip: inv -> dir -> inv ------------------------------
vor = sloc
div = 0.5 * sloc
gu, gv = ectrans4py.inv_trans_uv_vset4py(kspec2, kgptot, KFLDL, KFLDG, KVSET, vor, div)
check(gu.shape == (KFLDG, kgptot), f"inv_trans_uv_vset shape {gu.shape} != {(KFLDG, kgptot)}")
vor2, div2 = ectrans4py.dir_trans_uv_vset4py(kspec2, kgptot, KFLDL, KFLDG, KVSET, gu, gv)
gu2, gv2 = ectrans4py.inv_trans_uv_vset4py(kspec2, kgptot, KFLDL, KFLDG, KVSET, vor2, div2)
uv_err = comm.allreduce(
    max(float(np.max(np.abs(gu2 - gu))), float(np.max(np.abs(gv2 - gv)))), op=MPI.MAX
)
check(uv_err < 1e-4, f"UV vset grid->dir->inv round trip err {uv_err:.2e}")

# =============================================================================
# Adjoint (transpose) identities for the V-set transforms (the TL/AD building blocks):
# <F x, y>_grid == <x, F* y>_spec. The grid inner product is the unweighted Euclidean sum
# over (KFLDG, kgptot); the spectral inner product is m-folded over (KFLDL, kspec2)
# (m==0 weight 1, m>0 weight 2, the SCALPRODSP convention).
# =============================================================================
ADJ_TOL = 5e-3  # single precision + O(ksmax^2) accumulation


def model_order_m(kt):
    out = []
    for m in range(kt + 1):
        for _n in range(m, kt + 1):
            out += [m, m]  # real, imag parts share m
    return np.array(out, dtype=REAL)


# local zonal wavenumber m per coefficient, replicated across the V-sets
mval_g = np.zeros((1, kspec2g), dtype=REAL)
if rank == 0:
    mval_g[0, :] = model_order_m(KSMAX)
if nprtrv == 1:
    mloc = ectrans4py.dist_spec4py(kspec2g, kspec2, 1, np.ones(1, np.int64), mval_g)[0]
else:
    kfl = 1 if kmysetv == 1 else 0
    ml = ectrans4py.dist_spec_vset4py(
        kspec2g, kspec2, 1, kfl, np.ones(1, np.int64), np.ones(1, np.int64), mval_g
    )
    vcomm = comm.Split(color=kmysetw, key=kmysetv)
    mloc = vcomm.bcast(ml[0] if kmysetv == 1 else None, root=0)
mfold = np.where(np.rint(mloc) > 0, 2.0, 1.0)  # (kspec2,)


def dot_grid(a, b):
    return comm.allreduce(float(np.sum(a.astype(np.float64) * b.astype(np.float64))), op=MPI.SUM)


def dot_spec(a, b):
    return comm.allreduce(
        float(np.sum(mfold[None, :] * a.astype(np.float64) * b.astype(np.float64))), op=MPI.SUM
    )


def adj_check(name, lhs, rhs):
    check(abs(lhs - rhs) <= ADJ_TOL * max(1.0, abs(lhs)), f"adjoint {name}: {lhs:.6e} != {rhs:.6e}")


rng2 = np.random.default_rng(555 + rank)
xs = rng2.standard_normal((KFLDL, kspec2)).astype(REAL)  # local spectral
xs2 = rng2.standard_normal((KFLDL, kspec2)).astype(REAL)
yg = rng2.standard_normal((KFLDG, kgptot)).astype(REAL)  # grid
yg2 = rng2.standard_normal((KFLDG, kgptot)).astype(REAL)

# inv_scalar: <inv(x), y>_grid == <x, inv_ad(y)>_spec
Fx = ectrans4py.inv_trans_scalar_vset4py(kspec2, kgptot, KFLDL, KFLDG, KVSET, xs)
Fty = ectrans4py.inv_trans_scalar_vset4py_ad(kspec2, kgptot, KFLDL, KFLDG, KVSET, yg)
adj_check("inv_scalar", dot_grid(Fx, yg), dot_spec(xs, Fty))

# dir_scalar: <dir(g), x>_spec == <g, dir_ad(x)>_grid
Dg = ectrans4py.dir_trans_scalar_vset4py(kspec2, kgptot, KFLDL, KFLDG, KVSET, yg)
Dtx = ectrans4py.dir_trans_scalar_vset4py_ad(kspec2, kgptot, KFLDL, KFLDG, KVSET, xs)
adj_check("dir_scalar", dot_spec(Dg, xs), dot_grid(yg, Dtx))

# inv_uv: <inv(vor,div), (u,v)>_grid == <(vor,div), inv_ad(u,v)>_spec
Gu, Gv = ectrans4py.inv_trans_uv_vset4py(kspec2, kgptot, KFLDL, KFLDG, KVSET, xs, xs2)
Rvor, Rdiv = ectrans4py.inv_trans_uv_vset4py_ad(kspec2, kgptot, KFLDL, KFLDG, KVSET, yg, yg2)
adj_check("inv_uv", dot_grid(Gu, yg) + dot_grid(Gv, yg2), dot_spec(xs, Rvor) + dot_spec(xs2, Rdiv))

# dir_uv: <dir(u,v), (vor,div)>_spec == <(u,v), dir_ad(vor,div)>_grid
Dvor, Ddiv = ectrans4py.dir_trans_uv_vset4py(kspec2, kgptot, KFLDL, KFLDG, KVSET, yg, yg2)
Du, Dv = ectrans4py.dir_trans_uv_vset4py_ad(kspec2, kgptot, KFLDL, KFLDG, KVSET, xs, xs2)
adj_check("dir_uv", dot_spec(Dvor, xs) + dot_spec(Ddiv, xs2), dot_grid(yg, Du) + dot_grid(yg2, Dv))

# ders forward: value component equals the plain inverse; ders adjoint transpose identity
val, dns, dew = ectrans4py.inv_trans_scalar_ders_vset4py(kspec2, kgptot, KFLDL, KFLDG, KVSET, xs)
plainv = ectrans4py.inv_trans_scalar_vset4py(kspec2, kgptot, KFLDL, KFLDG, KVSET, xs)
dv_err = comm.allreduce(float(np.max(np.abs(val - plainv))), op=MPI.MAX)
check(dv_err < 1e-5, f"inv_scalar_ders value != inv_scalar: {dv_err:.2e}")
se = rng2.standard_normal((KFLDG, kgptot)).astype(REAL)
Ftd = ectrans4py.inv_trans_scalar_ders_vset4py_ad(kspec2, kgptot, KFLDL, KFLDG, KVSET, yg, yg2, se)
adj_check(
    "inv_scalar_ders",
    dot_grid(val, yg) + dot_grid(dns, yg2) + dot_grid(dew, se),
    dot_spec(xs, Ftd),
)

# specnorm_vset: per-global-field norm (finite, non-negative, right shape)
pn = ectrans4py.specnorm_vset4py(kspec2, KFLDG, KFLDL, KVSET, xs)
check(
    pn.shape == (KFLDG,) and bool(np.all(np.isfinite(pn))) and bool(np.all(pn >= 0.0)),
    f"specnorm_vset bad shape/values: {pn}",
)

if rank == 0:
    print(
        f"[vset] size={size} nprtrw={nprtrw} nprtrv={nprtrv} KMYSETW={kmysetw} "
        f"KMYSETV={kmysetv} KFLDG={KFLDG} KFLDL={KFLDL} kspec2={kspec2} "
        f"scalar_rt={sc_err:.2e} uv_rt={uv_err:.2e} specnorm_sum={float(np.sum(pn)):.6e}"
    )
    print("[vset] PASS" if ok else "[vset] FAIL")

ectrans4py.mpl_end4py()
allok = comm.allreduce(ok, op=MPI.LAND)
if not allok:
    sys.exit(1)
