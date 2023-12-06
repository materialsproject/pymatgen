from __future__ import annotations

import json
import unittest
from shutil import which

from monty.serialization import loadfn
from pytest import approx

from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer, BoltztrapRunner
from pymatgen.electronic_structure.core import OrbitalType, Spin
from pymatgen.util.testing import TEST_FILES_DIR

try:
    from ase.io.cube import read_cube
except ImportError:
    read_cube = None

try:
    import fdint
except ImportError:
    fdint = None

x_trans = which("x_trans")


@unittest.skipIf(not x_trans, "No x_trans.")
class TestBoltztrapAnalyzer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bz = BoltztrapAnalyzer.from_files(f"{TEST_FILES_DIR}/boltztrap/transp/")
        cls.bz_bands = BoltztrapAnalyzer.from_files(f"{TEST_FILES_DIR}/boltztrap/bands/")
        cls.bz_up = BoltztrapAnalyzer.from_files(f"{TEST_FILES_DIR}/boltztrap/dos_up/", dos_spin=1)
        cls.bz_dw = BoltztrapAnalyzer.from_files(f"{TEST_FILES_DIR}/boltztrap/dos_dw/", dos_spin=-1)
        cls.bz_fermi = BoltztrapAnalyzer.from_files(f"{TEST_FILES_DIR}/boltztrap/fermi/")

        with open(f"{TEST_FILES_DIR}/Cu2O_361_bandstructure.json") as file:
            d = json.load(file)
            cls.bs = BandStructure.from_dict(d)
            cls.btr = BoltztrapRunner(cls.bs, 1)

    def test_properties(self):
        assert self.bz.gap == approx(1.6644932121620404, abs=1e-4)
        array = self.bz._cond[300][102]
        assert array[0][0] / 1e19 == approx(7.5756518, abs=1e-4)
        assert array[0][2] == approx(-11.14679)
        assert array[1][0] == approx(-88.203286)
        assert array[2][2] == approx(1.7133249e19)
        array = self.bz._seebeck[300][22]
        assert array[0][1] == approx(6.4546074e-22)
        assert array[1][1] == approx(-0.00032073711)
        assert array[1][2] == approx(-2.9868424e-24)
        assert array[2][2] == approx(-0.0003126543)
        array = self.bz._kappa[500][300]
        assert array[0][1] == approx(0.00014524309)
        assert array[1][1] == approx(328834400000000.0)
        assert array[1][2] == approx(3.7758069e-05)
        assert array[2][2] == approx(193943750000000.0)
        assert self.bz._hall[400][800][1][0][0] == approx(9.5623749e-28)
        assert self.bz._hall[400][68][1][2][2] == approx(3.2149951e-26)
        assert self.bz.doping["p"][3] == approx(1e18)
        assert self.bz.mu_doping["p"][300][2] == approx(0.1553770018406)
        assert self.bz.mu_doping["n"][300][-1] == approx(1.6486017632924719, abs=1e-4)
        assert self.bz._cond_doping["n"][800][3][1][1] / 1e16 == approx(1.5564085, abs=1e-4)
        assert self.bz._seebeck_doping["p"][600][2][0][1] / 1e-23 == approx(3.2860613, abs=1e-4)
        assert self.bz._carrier_conc[500][67] == approx(38.22832002)
        assert self.bz.vol == approx(612.97557323964838, abs=1e-4)
        assert self.bz.intrans["scissor"] == approx(0.0, abs=1e-1)
        assert self.bz._hall_doping["n"][700][-1][2][2][2] == approx(5.0136483e-26)
        assert self.bz.dos.efermi == approx(-0.0300005507057)
        assert self.bz.dos.energies[0] == approx(-2.4497049391830448, abs=1e-4)
        assert self.bz.dos.energies[345] == approx(-0.72708823447130944, abs=1e-4)
        assert self.bz.dos.energies[-1] == approx(3.7569398770153524, abs=1e-4)
        assert self.bz.dos.densities[Spin.up][400] == approx(118.70171)
        assert self.bz.dos.densities[Spin.up][200] == approx(179.58562)
        assert self.bz.dos.densities[Spin.up][300] == approx(289.43945)

        assert self.bz_bands._bz_bands.shape == approx((1316, 20))
        assert self.bz_bands._bz_kpoints.shape == approx((1316, 3))
        assert self.bz_up._dos_partial["0"]["pz"][2562] == approx(0.023862958)
        assert self.bz_dw._dos_partial["1"]["px"][3120] == approx(5.0192891)
        assert self.bz_fermi.fermi_surface_data.shape == approx((121, 121, 65))
        assert self.bz_fermi.fermi_surface_data[21][79][19] == approx(-1.8831911809439161, abs=1e-5)

    @unittest.skipIf(not fdint, "No FDINT")
    def test_get_seebeck_eff_mass(self):
        ref = [1.956090529381193, 2.0339311618566343, 1.1529383757896965]
        ref2 = [4258.4072823354145, 4597.0351887125289, 4238.1262696392705]
        sbk_mass_tens_mu = self.bz.get_seebeck_eff_mass(output="tensor", doping_levels=False, temp=300)[3]
        sbk_mass_tens_dop = self.bz.get_seebeck_eff_mass(output="tensor", doping_levels=True, temp=300)["n"][2]
        sbk_mass_avg_mu = self.bz.get_seebeck_eff_mass(output="average", doping_levels=False, temp=300)[3]
        sbk_mass_avg_dop = self.bz.get_seebeck_eff_mass(output="average", doping_levels=True, temp=300)["n"][2]

        for i in range(3):
            assert sbk_mass_tens_mu[i] == approx(ref2[i], abs=1e-1)
            assert sbk_mass_tens_dop[i] == approx(ref[i], abs=1e-4)

        assert sbk_mass_avg_mu == approx(4361.4744008038842, abs=1e-1)
        assert sbk_mass_avg_dop == approx(1.661553842105382, abs=1e-4)

    @unittest.skipIf(not fdint, "No FDINT")
    def test_get_complexity_factor(self):
        ref = [2.7658776815227828, 2.9826088215568403, 0.28881335881640308]
        ref2 = [0.0112022048620205, 0.0036001049607186602, 0.0083028947173193028]
        sbk_mass_tens_mu = self.bz.get_complexity_factor(output="tensor", doping_levels=False, temp=300)[3]
        sbk_mass_tens_dop = self.bz.get_complexity_factor(output="tensor", doping_levels=True, temp=300)["n"][2]
        sbk_mass_avg_mu = self.bz.get_complexity_factor(output="average", doping_levels=False, temp=300)[3]
        sbk_mass_avg_dop = self.bz.get_complexity_factor(output="average", doping_levels=True, temp=300)["n"][2]

        for i in range(3):
            assert sbk_mass_tens_mu[i] == approx(ref2[i], abs=1e-4)
            assert sbk_mass_tens_dop[i] == approx(ref[i], abs=1e-4)

        assert sbk_mass_avg_mu == approx(0.00628677029221, abs=1e-4)
        assert sbk_mass_avg_dop == approx(1.12322832119, abs=1e-4)

    def test_get_seebeck(self):
        ref = [-768.99078999999995, -724.43919999999991, -686.84682999999973]
        for i in range(3):
            assert self.bz.get_seebeck()["n"][800][3][i] == approx(ref[i])
        assert self.bz.get_seebeck(output="average")["p"][800][3] == approx(697.608936667)
        assert self.bz.get_seebeck(output="average", doping_levels=False)[500][520] == approx(1266.7056)
        assert self.bz.get_seebeck(output="average", doping_levels=False)[300][65] == approx(
            -36.2459389333
        )  # TODO: this was originally "eigs"

    def test_get_conductivity(self):
        ref = [5.9043185000000022, 17.855599000000002, 26.462935000000002]
        for i in range(3):
            assert self.bz.get_conductivity()["p"][600][2][i] == approx(ref[i])
        assert self.bz.get_conductivity(output="average")["n"][700][1] == approx(1.58736609667)
        assert self.bz.get_conductivity(output="average", doping_levels=False)[300][457] == approx(2.87163566667)
        assert self.bz.get_conductivity(
            output="average",
            doping_levels=False,
            # TODO: this was originally "eigs"
            relaxation_time=1e-15,
        )[200][63] == approx(16573.0536667)

    def test_get_power_factor(self):
        ref = [6.2736602345523362, 17.900184232304138, 26.158282220458144]
        for i in range(3):
            assert self.bz.get_power_factor()["p"][200][2][i] == approx(ref[i])
        assert self.bz.get_power_factor(output="average")["n"][600][4] == approx(411.230962976)
        assert self.bz.get_power_factor(output="average", doping_levels=False, relaxation_time=1e-15)[500][
            459
        ] == approx(6.59277148467)
        assert self.bz.get_power_factor(output="average", doping_levels=False)[800][61] == approx(
            2022.67064134
        )  # TODO: this was originally "eigs"

    def test_get_thermal_conductivity(self):
        ref = [2.7719565628862623e-05, 0.00010048046886793946, 0.00015874549392499391]
        for i in range(3):
            assert self.bz.get_thermal_conductivity()["p"][300][2][i] == approx(ref[i])
        assert self.bz.get_thermal_conductivity(output="average", relaxation_time=1e-15)["n"][500][0] == approx(
            1.74466575612e-07
        )
        assert self.bz.get_thermal_conductivity(output="average", doping_levels=False)[800][874] == approx(
            8.08066254813
        )
        assert self.bz.get_thermal_conductivity(output="average", doping_levels=False)[200][32] == approx(
            0.0738961845832
        )
        assert self.bz.get_thermal_conductivity(k_el=False, output="average", doping_levels=False)[200][32] == approx(
            0.19429052
        )

    def test_get_zt(self):
        ref = [0.097408810215, 0.29335112354, 0.614673998089]
        for i in range(3):
            assert self.bz.get_zt()["n"][800][4][i] == approx(ref[i])
        assert self.bz.get_zt(output="average", k_l=0.5)["p"][700][2] == approx(0.0170001879916)
        assert self.bz.get_zt(output="average", doping_levels=False, relaxation_time=1e-15)[300][240] == approx(
            0.0041923533238348342
        )

        eigs = self.bz.get_zt(output="eigs", doping_levels=False)[700][65]
        ref_eigs = [0.082420053399668847, 0.29408035502671648, 0.40822061215079392]
        for idx, val in enumerate(ref_eigs):
            assert eigs[idx] == approx(val, abs=1e-5)

    def test_get_average_eff_mass(self):
        ref = [0.76045816788363574, 0.96181142990667101, 2.9428428773308628]
        for i in range(3):
            assert self.bz.get_average_eff_mass()["p"][300][2][i] == approx(ref[i])
        ref = [1.1295783824744523, 1.3898454041924351, 5.2459984671977935]
        ref2 = [6.6648842712692078, 31.492540105738343, 37.986369302138954]
        for i in range(3):
            assert self.bz.get_average_eff_mass()["n"][600][1][i] == approx(ref[i])
            assert self.bz.get_average_eff_mass(doping_levels=False)[300][200][i] == approx(ref2[i])
        ref = [
            [9.61811430e-01, -8.25159596e-19, -4.70319444e-19],
            [-8.25159596e-19, 2.94284288e00, 3.00368916e-18],
            [-4.70319444e-19, 3.00368916e-18, 7.60458168e-01],
        ]
        ref2 = [
            [27.97604444269153, -2.39347589e-17, -1.36897140e-17],
            [-2.39347589e-17, 8.55969097e01, 8.74169648e-17],
            [-1.36897140e-17, 8.74169648e-17, 2.21151980e01],
        ]

        for i in range(3):
            for j in range(3):
                assert self.bz.get_average_eff_mass(output="tensor")["p"][300][2][i][j] == approx(ref[i][j], abs=1e-4)
                assert self.bz.get_average_eff_mass(output="tensor", doping_levels=False)[300][500][i][j] == approx(
                    ref2[i][j], 4
                )
        assert self.bz.get_average_eff_mass(output="average")["n"][300][2] == approx(1.53769093989, abs=1e-4)

    def test_get_carrier_concentration(self):
        assert self.bz.get_carrier_concentration()[300][39] / 1e22 == approx(6.4805156617179151, abs=1e-4)
        assert self.bz.get_carrier_concentration()[300][693] / 1e15 == approx(-6.590800965604750, abs=1e-4)

    def test_get_hall_carrier_concentration(self):
        assert self.bz.get_hall_carrier_concentration()[600][120] / 1e21 == approx(6.773394626767555, abs=1e-4)
        assert self.bz.get_hall_carrier_concentration()[500][892] / 1e21 == approx(-9.136803845741777, abs=1e-4)

    def test_get_symm_bands(self):
        structure = loadfn(f"{TEST_FILES_DIR}/boltztrap/structure_mp-12103.json")
        sbs = loadfn(f"{TEST_FILES_DIR}/boltztrap/dft_bs_sym_line.json")
        kpoints = [kp.frac_coords for kp in sbs.kpoints]
        labels_dict = {k: sbs.labels_dict[k].frac_coords for k in sbs.labels_dict}
        for kpt_line, label_dict in zip([None, sbs.kpoints, kpoints], [None, sbs.labels_dict, labels_dict]):
            sbs_bzt = self.bz_bands.get_symm_bands(structure, -5.25204548, kpt_line=kpt_line, labels_dict=label_dict)
            assert len(sbs_bzt.bands[Spin.up]) == approx(20)
            assert len(sbs_bzt.bands[Spin.up][1]) == approx(143)

    # def test_check_acc_bzt_bands(self):
    #     structure = loadfn(f"{TEST_FILES_DIR}/boltztrap/structure_mp-12103.json")
    #     sbs = loadfn(f"{TEST_FILES_DIR}/boltztrap/dft_bs_sym_line.json")
    #     sbs_bzt = self.bz_bands.get_symm_bands(structure, -5.25204548)
    #     corr, werr_vbm, werr_cbm, warn = BoltztrapAnalyzer.check_acc_bzt_bands(sbs_bzt, sbs)
    #     assert corr[2] == 9.16851750e-05
    #     assert werr_vbm["K-H"] == 0.18260273521047862
    #     assert werr_cbm["M-K"] == 0.071552669981356981
    #     assert not warn

    def test_get_complete_dos(self):
        structure = loadfn(f"{TEST_FILES_DIR}/boltztrap/structure_mp-12103.json")
        cdos = self.bz_up.get_complete_dos(structure, self.bz_dw)
        spins = list(cdos.densities)
        assert Spin.down in spins
        assert Spin.up in spins
        assert cdos.get_spd_dos()[OrbitalType.p].densities[Spin.up][3134] == approx(43.839230100999991)
        assert cdos.get_spd_dos()[OrbitalType.s].densities[Spin.down][716] == approx(6.5383268000000001)

    def test_extreme(self):
        extreme = self.bz.get_extreme("seebeck")
        assert extreme["best"]["carrier_type"] == "n"
        assert extreme["p"]["value"] == approx(1255.365, abs=1e-2)
        assert extreme["n"]["isotropic"]
        assert extreme["n"]["temperature"] == 600

        extreme = self.bz.get_extreme("kappa", maximize=False, min_temp=400, min_doping=1e20)
        assert extreme["best"]["value"] == approx(0.105, abs=1e-2)
        assert extreme["n"]["value"] == approx(0.139, abs=1e-2)
        assert extreme["p"]["temperature"] == 400
        assert extreme["n"]["isotropic"] is False

    def test_as_from_dict(self):
        btr_dict = self.btr.as_dict()
        json_str = json.dumps(btr_dict)
        assert json_str is not None
        assert btr_dict["bs"] is not None
