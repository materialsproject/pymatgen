"""
Calculate spin-orbit spillage given the wavefunctions with and without spin-orbit coupling

Forked from: https://github.com/usnistgov/jarvis
Authors: Kamal Choudhary, Kevin Garrity (NIST)
References: 1) https://doi.org/10.1038/s41598-019-45028-y 
            2) https://doi.org/10.1103/PhysRevB.90.125133
"""

import os, sys
import numpy as np
from pymatgen.analysis.topological.vaspwfc import vaspwfc
import scipy as sp


def isclose(n1, n2, rel_tol=1e-7):
    if abs(n1 - n2) < rel_tol:
        return True
    else:
        return False


def orth(A):
    u, s, vh = np.linalg.svd(A, full_matrices=False)
    M, N = A.shape
    eps = np.finfo(float).eps
    tol = max(M, N) * np.amax(s) * eps
    num = np.sum(s > tol, dtype=int)
    Q = u[:, :num]
    return Q, num


def overlap_so_spinpol(wf_noso, wf_so):
    noso = vaspwfc(fnm=wf_noso, lsorbit=False)
    so = vaspwfc(fnm=wf_so, lsorbit=True)
    noso_k, noso_bands = noso.readWFBand()
    so_k, so_bands = so.readWFBand()

    nelec_list = []
    for nk1 in range(1, noso._nkpts + 1):  # no spin orbit kpoints loop
        knoso = noso._kvecs[nk1 - 1, :]
        for nk2 in range(1, so._nkpts + 1):  # spin orbit
            kso = so._kvecs[nk2 - 1, :]
            if (
                isclose(kso[0], knoso[0])
                and isclose(kso[1], knoso[1])
                and isclose(kso[2], knoso[2])
            ):  # do kpoints match?
                for c, e in enumerate(noso._occs[0, nk1 - 1, :]):
                    if e < 0.5:
                        cup = c
                        break
                for c, e in enumerate(noso._occs[1, nk1 - 1, :]):
                    if e < 0.5:
                        cdn = c
                        break

                nelec_list.append([cup, cdn, cup + cdn])

    n_arr = np.array(nelec_list)

    n_up = int(round(np.mean(n_arr[:, 0])))
    n_dn = int(round(np.mean(n_arr[:, 1])))
    n_tot = int(round(np.mean(n_arr[:, 2])))

    nelec = int(n_tot)

    noso_homo_up = np.max(noso_bands[0, :, n_up - 1])
    noso_lumo_up = np.min(noso_bands[0, :, n_up])

    noso_homo_dn = np.max(noso_bands[1, :, n_dn - 1])
    noso_lumo_dn = np.min(noso_bands[1, :, n_dn])

    so_homo = np.max(so_bands[0, :, nelec - 1])
    so_lumo = np.min(so_bands[0, :, nelec])

    noso_direct_up = np.min(noso_bands[0, :, n_up] - noso_bands[0, :, n_up - 1])
    noso_direct_dn = np.min(noso_bands[1, :, n_dn] - noso_bands[1, :, n_dn - 1])

    so_direct = np.min(so_bands[0, :, nelec] - so_bands[0, :, nelec - 1])

    noso_direct = 1000000.0
    noso_homo = -10000000.0
    noso_lumo = 100000000.0
    for i in range(noso_bands.shape[1]):
        homo_k = max(noso_bands[0, i, n_up - 1], noso_bands[1, i, n_dn - 1])
        lumo_k = min(noso_bands[0, i, n_up], noso_bands[1, i, n_dn])
        noso_direct = min(noso_direct, lumo_k - homo_k)

        noso_homo = max(noso_homo, homo_k)
        noso_lumo = min(noso_lumo, lumo_k)

    gamma_k = []
    kpoints = []
    np.set_printoptions(precision=4)

    x = []
    y = []

    nelec_tot = 0.0
    for nk1 in range(1, noso._nkpts + 1):  # no spin orbit kpoints loop
        knoso = noso._kvecs[nk1 - 1, :]
        for nk2 in range(1, so._nkpts + 1):  # spin orbit
            kso = so._kvecs[nk2 - 1, :]
            if (
                isclose(kso[0], knoso[0])
                and isclose(kso[1], knoso[1])
                and isclose(kso[2], knoso[2])
            ):  # do kpoints match?

                # changes section 2
                nelec_up = n_arr[nk1 - 1, 0]
                nelec_dn = n_arr[nk1 - 1, 1]
                nelec_tot = n_arr[nk1 - 1, 2]

                gamma_k.append(nelec_tot)
                kpoints.append(kso)
                Mmn = 0.0
                vnoso = noso.readBandCoeff(ispin=1, ikpt=nk1, iband=1, norm=False)
                n_noso1 = vnoso.shape[0]
                vnoso = noso.readBandCoeff(ispin=2, ikpt=nk1, iband=1, norm=False)
                n_noso2 = vnoso.shape[0]
                vso = so.readBandCoeff(ispin=1, ikpt=nk2, iband=1, norm=False)
                n_so = vso.shape[0]
                vs = min(n_noso1 * 2, n_so)
                Vnoso = np.zeros((vs, nelec_tot), dtype=complex)
                Vso = np.zeros((vs, nelec_tot), dtype=complex)

                # prepare matricies
                for n1 in range(1, nelec_up + 1):
                    Vnoso[0 : vs // 2, n1 - 1] = noso.readBandCoeff(
                        ispin=1, ikpt=nk1, iband=n1, norm=False
                    )[0 : vs // 2]
                for n1 in range(1, nelec_dn + 1):
                    Vnoso[vs // 2 : vs, n1 - 1 + nelec_up] = noso.readBandCoeff(
                        ispin=2, ikpt=nk1, iband=n1, norm=False
                    )[0 : vs // 2]

                for n1 in range(1, nelec_tot + 1):
                    t = so.readBandCoeff(ispin=1, ikpt=nk2, iband=n1, norm=False)
                    Vso[0 : vs // 2, n1 - 1] = t[0 : vs // 2]
                    Vso[vs // 2 : vs, n1 - 1] = t[n_so // 2 : n_so // 2 + vs // 2]
                Qnoso, num_noso = orth(Vnoso)  # make orthonormal basis?

                Qso, num_so = orth(Vso)

                # evalute main expression

                a = []
                for n1 in range(0, nelec_tot):  # noso occupied bands
                    v1 = Qnoso[:, n1]
                    aa = 0.0
                    for n2 in range(0, nelec_tot):  # so occupied bands
                        v2 = Qso[:, n2]

                        t = np.dot(np.conj(v1), v2)
                        Mmn += t * t.conj()
                        aa += (t * t.conj()).real
                    a.append(aa)
                gamma_k[-1] -= Mmn  # eq 4 in prb 90 125133
                x.append(kso)
                y.append(np.real(gamma_k[-1]))

    gmax = max(np.real(gamma_k))
    nkmax = np.argmax(np.real(gamma_k))
    kmax = kpoints[nkmax]

    print()
    print("                   INDIRECT DIRECT      HOMO/LUMO (eV)")
    print(
        "no spin-orbit gaps",
        "{:+.3f}".format(float(noso_lumo - noso_homo)),
        "{:+.3f}".format(noso_direct),
        "   ",
        [noso_homo, noso_lumo],
    )
    print(
        "spin-orbit gaps   ",
        "{:+.3f}".format(float(so_lumo - so_homo)),
        "{:+.3f}".format(so_direct),
        "   ",
        [so_homo, so_lumo],
    )
    print()

    return np.real(gmax)


if __name__ == "__main__":
    wf_noso = os.path.join(os.path.dirname(__file__), "tests", "WAVECAR-NonSOC")
    wf_so = os.path.join(os.path.dirname(__file__), "tests", "WAVECAR-SOC")
    gamma_max = overlap_so_spinpol(wf_noso, wf_so)
    print("spillage = ", gamma_max)
