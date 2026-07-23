"""MPL_INIT_COMM4PY: bind FIAT MPL (and thus all of ecTrans) to a *user*
sub-communicator instead of MPI_COMM_WORLD -- the mechanism that lets ecTrans live on
a compute sub-group when the world is split (e.g. compute vs IO-server ranks).

Run:  mpirun -n 4 python tests/test_ectrans4py/mpi_distributed_comm.py   (even n >= 2)

Splits COMM_WORLD into two colour groups; each group binds MPL to its own
sub-communicator via ``mpl_init_comm4py``, then sets up an independent distributed
transform over just that group and runs a grid -> dir -> inv -> grid round trip.
Confirms MPL reports the *sub-group* size (not the world size) and that the transforms
are collective over the sub-communicator only.
"""

import sys

import numpy as np
from mpi4py import MPI

import ectrans4py

world = MPI.COMM_WORLD
wrank, wsize = world.Get_rank(), world.Get_size()

if wsize % 2 != 0:
    if wrank == 0:
        print(f"[comm] SKIP: need an even task count (got {wsize})")
    sys.exit(0)

# two colour groups; ecTrans will run independently on each
colour = wrank % 2
sub = world.Split(color=colour, key=wrank)
ssize = sub.Get_size()

ectrans4py.init_env(unlimited_stack=False)
# bind MPL to the sub-communicator (mpi4py Fortran handle)
rank1, size1 = ectrans4py.mpl_init_comm4py(sub.py2f())

ok = True


def check(cond, msg):
    global ok
    if not cond:
        ok = False
        print(f"[comm] world-rank {wrank} FAIL: {msg}")


# MPL must see the sub-group, not the world
check(size1 == ssize, f"MPL size {size1} != sub-group size {ssize} (bound to world?)")
check(1 <= rank1 <= ssize, f"MPL rank {rank1} out of [1,{ssize}]")

# independent distributed transform over the sub-group
NDGL, KSMAX = 32, 21
KLOEN = np.full(NDGL, 2 * NDGL, dtype=np.int64)
ectrans4py.setup_trans0_4py(ssize, 1, ssize, False, 10, False)
kresol = ectrans4py.setup_trans_4py(KSMAX, NDGL, NDGL, KLOEN, True, False)
(kgptot, kspec, kspec2, kgptotg, kspec2g, ksmax, knmeng, pmu, pgw) = ectrans4py.trans_inq4py(
    kresol, NDGL, KSMAX, NDGL, KLOEN, 10
)

# local grid counts must sum to the global over the SUB-group (not the world)
check(
    sub.allreduce(int(kgptot), op=MPI.SUM) == kgptotg,
    "local kgptot do not sum to global over the sub-communicator",
)

# grid -> dir -> inv -> grid identity, confined to the sub-group
nfld = 1
rng = np.random.default_rng(100 + colour)
g = rng.standard_normal((nfld, kgptot)).astype(getattr(ectrans4py, "_REAL", np.float64))
s = ectrans4py.dir_trans_scalar_dist4py(kspec2, kgptot, nfld, g)
g2 = ectrans4py.inv_trans_scalar_dist4py(kspec2, kgptot, nfld, s)
s2 = ectrans4py.dir_trans_scalar_dist4py(kspec2, kgptot, nfld, g2)
err = sub.allreduce(float(np.max(np.abs(s2 - s))), op=MPI.MAX)
check(err < 1e-5, f"sub-group dir->inv->dir round trip err {err:.2e}")

if wrank == 0:
    print(f"[comm] world={wsize} groups=2 sub_size={ssize} MPL_size={size1} rt_err={err:.2e}")

ectrans4py.mpl_end4py()
allok = world.allreduce(ok, op=MPI.LAND)
if wrank == 0:
    print("[comm] PASS" if allok else "[comm] FAIL")
if not allok:
    sys.exit(1)
