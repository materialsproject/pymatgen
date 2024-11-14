"""
Code to calculate spin-orbit spillage.
Modified from JARVIS-Tools
https://www.nature.com/articles/s41598-019-45028-y
https://www.nature.com/articles/s41524-020-0319-4.
"""

from __future__ import annotations

import numpy as np

from pymatgen.io.vasp.outputs import Wavecar


class SOCSpillage:
    """
    Spin-orbit spillage criteria to predict whether a material is topologically non-trivial.
    The spillage criteria physically signifies number of band-inverted electrons.
    A non-zero, high value (generally >0.5) suggests non-trivial behavior.
    """

    def __init__(self, wf_noso="", wf_so=""):
        """
        Requires path to WAVECAR files with and without LSORBIT = .TRUE.

        Args:
            wf_noso : WAVECAR without spin-orbit coupling
            wf_so : WAVECAR with spin-orbit coupling
        """
        self.wf_noso = wf_noso
        self.wf_so = wf_so

    @staticmethod
    def isclose(n1, n2, rel_tol=1e-7):
        """Checking if the numbers are close enough."""
        return abs(n1 - n2) < rel_tol

    @staticmethod
    def orth(A):
        """Helper function to create orthonormal basis."""
        u, s, _vh = np.linalg.svd(A, full_matrices=False)
        n_rows, n_cols = A.shape
        eps = np.finfo(float).eps
        tol = max(n_rows, n_cols) * np.amax(s) * eps
        num = np.sum(s > tol, dtype=np.int64)
        orthonormal_basis = u[:, :num]
        return orthonormal_basis, num

    def overlap_so_spinpol(self):
        """Main function to calculate SOC spillage."""
        no_so = Wavecar(self.wf_noso)
        so = Wavecar(self.wf_so)

        b_cell = np.linalg.inv(no_so.a).T
        tmp = np.linalg.norm(np.dot(np.diff(no_so.kpoints, axis=0), b_cell), axis=1)
        noso_k = np.concatenate(([0], np.cumsum(tmp)))
        noso_bands = np.array(no_so.band_energy)[:, :, :, 0]
        noso_kvecs = np.array(no_so.kpoints)
        noso_occs = np.array(no_so.band_energy)[:, :, :, 2]
        n_kpts_noso = len(noso_k)

        b_cell = np.linalg.inv(so.a).T
        tmp = np.linalg.norm(np.dot(np.diff(so.kpoints, axis=0), b_cell), axis=1)
        so_k = np.concatenate(([0], np.cumsum(tmp)))
        so_bands = np.array([np.array(so.band_energy)[:, :, 0]])
        so_kvecs = np.array(so.kpoints)
        # so_occs = np.array([np.array(so.band_energy)[:, :, 2]])
        so_nkpts = len(so_k)

        n_elec_list: list[list[int]] = []
        cup = cdn = 0
        for nk1 in range(1, n_kpts_noso + 1):  # no spin orbit kpoints loop
            k_no_spin_orbit = noso_kvecs[nk1 - 1, :]
            for nk2 in range(1, so_nkpts + 1):  # spin orbit
                kso = so_kvecs[nk2 - 1, :]
                if (
                    self.isclose(kso[0], k_no_spin_orbit[0])
                    and self.isclose(kso[1], k_no_spin_orbit[1])
                    and self.isclose(kso[2], k_no_spin_orbit[2])
                ):  # do kpoints match?
                    for c, e in enumerate(noso_occs[0, nk1 - 1, :]):
                        if e < 0.5:
                            cup = c
                            break
                    for c, e in enumerate(noso_occs[1, nk1 - 1, :]):
                        if e < 0.5:
                            cdn = c
                            break

                    n_elec_list.append([cup, cdn, cup + cdn])
        n_arr = np.array(n_elec_list)

        n_up = round(np.mean(n_arr[:, 0]))
        n_dn = round(np.mean(n_arr[:, 1]))
        n_tot = round(np.mean(n_arr[:, 2]))

        n_elec = n_tot

        # noso_homo_up = np.max(noso_bands[0, :, n_up - 1])
        # noso_lumo_up = np.min(noso_bands[0, :, n_up])

        # noso_homo_dn = np.max(noso_bands[1, :, n_dn - 1])
        # noso_lumo_dn = np.min(noso_bands[1, :, n_dn])

        so_homo = np.max(so_bands[0, :, n_elec - 1])
        so_lumo = np.min(so_bands[0, :, n_elec])

        # noso_direct_up = np.min(noso_bands[0, :, n_up] - noso_bands[0, :, n_up - 1])
        # noso_direct_dn = np.min(noso_bands[1, :, n_dn] - noso_bands[1, :, n_dn - 1])

        so_direct = np.min(so_bands[0, :, n_elec] - so_bands[0, :, n_elec - 1])

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
        for nk1 in range(1, n_kpts_noso + 1):  # no spin orbit kpoints loop
            k_no_spin_orbit = noso_kvecs[nk1 - 1, :]
            for nk2 in range(1, so_nkpts + 1):  # spin orbit
                kso = so_kvecs[nk2 - 1, :]
                if (
                    self.isclose(kso[0], k_no_spin_orbit[0])
                    and self.isclose(kso[1], k_no_spin_orbit[1])
                    and self.isclose(kso[2], k_no_spin_orbit[2])
                ):  # do kpoints match?
                    # changes section 2
                    nelec_up = n_arr[nk1 - 1, 0]
                    nelec_dn = n_arr[nk1 - 1, 1]
                    nelec_tot = n_arr[nk1 - 1, 2]

                    kpoints.append(kso)
                    Mmn = 0.0
                    vnoso = np.array(
                        no_so.coeffs[0][nk1 - 1][0]
                    )  # noso.readBandCoeff(ispin=1, ikpt=nk1, iband=1, norm=False)
                    n_noso1 = vnoso.shape[0]
                    vnoso = np.array(
                        no_so.coeffs[1][nk1 - 1][0]
                    )  # noso.readBandCoeff(ispin=2, ikpt=nk1, iband=1, norm=False)
                    # n_noso2 = vnoso.shape[0]
                    vso = so.coeffs[nk1 - 1][0].flatten()  # so.readBandCoeff(ispin=1, ikpt=nk2, iband=1, norm=False)
                    n_so = vso.shape[0]

                    vs = min(n_noso1 * 2, n_so)
                    Vnoso = np.zeros((vs, nelec_tot), dtype=complex)
                    Vso = np.zeros((vs, nelec_tot), dtype=complex)

                    if np.array(no_so.coeffs[1][nk1 - 1]).shape[1] == vs // 2:
                        # if nk1==10 and nk2==10:
                        # prepare matrices
                        for n1 in range(1, nelec_up + 1):
                            Vnoso[0 : vs // 2, n1 - 1] = np.array(no_so.coeffs[0][nk1 - 1][n1 - 1])[0 : vs // 2]
                        for n1 in range(1, nelec_dn + 1):
                            Vnoso[vs // 2 : vs, n1 - 1 + nelec_up] = np.array(no_so.coeffs[1][nk1 - 1][n1 - 1])[
                                0 : vs // 2
                            ]
                        for n1 in range(1, nelec_tot + 1):
                            t = so.coeffs[nk2 - 1][n1 - 1].flatten()
                            Vso[0 : vs // 2, n1 - 1] = t[0 : vs // 2]
                            Vso[vs // 2 : vs, n1 - 1] = t[n_so // 2 : n_so // 2 + vs // 2]
                        Qnoso, _num_noso = self.orth(Vnoso)  # make orthonormal basis?

                        Qso, _num_so = self.orth(Vso)

                        gamma_k.append(nelec_tot)
                        a = []
                        for n1 in range(nelec_tot):  # noso occupied bands
                            v1 = Qnoso[:, n1]
                            aa = 0.0
                            for n2 in range(nelec_tot):  # so occupied bands
                                v2 = Qso[:, n2]

                                t = np.dot(np.conj(v1), v2)
                                Mmn += t * t.conj()
                                aa += (t * t.conj()).real
                            a.append(aa)
                        gamma_k[-1] -= Mmn  # eq 4 in prb 90 125133
                        x.append(kso)
                        y.append(np.real(gamma_k[-1]))
                        if gamma_k[-1] > 0.5:
                            print(
                                "nk1 nk2 kpoint gamma_k ",
                                nk1,
                                nk2,
                                kso,
                                k_no_spin_orbit,
                                np.real(gamma_k[-1]),
                                "!!!!!!!!!!",
                            )

        gamma_max = max(np.real(gamma_k))
        n_kmax = np.argmax(np.real(gamma_k))
        k_max = kpoints[n_kmax]

        print("------------------------------------")
        print("\n                   INDIRECT DIRECT      HOMO/LUMO (eV)")
        print(
            "no spin-orbit gaps",
            f"{float(noso_lumo - noso_homo):+.3f}",
            f"{noso_direct:+.3f}",
            "   ",
            [noso_homo, noso_lumo],
        )
        print(
            "spin-orbit gaps   ",
            f"{float(so_lumo - so_homo):+.3f}",
            f"{so_direct:+.3f}",
            "   ",
            [so_homo, so_lumo],
        )
        print("gamma max", np.real(gamma_max), " at k =  ", k_max)
        return gamma_max
