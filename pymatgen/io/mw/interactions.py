from __future__ import annotations

from monty.json import MSONable


class Interactions(MSONable):
    """Interactions class representing the interactions block in the simulation.

    Attributes:
        coulomb (CoulombInteraction): Coulomb interaction configuration.
        lennard_jones (LennardJonesInteraction): Lennard-Jones interaction configuration.
        fumi_tosi (FumiTosiInteraction): Fumi-Tosi interaction configuration.
        xft (XFTInteraction): XFT interaction configuration.
        daim (DaimInteraction): Daim interaction configuration.
        damping (DampingInteraction): Damping interaction configuration.
        steele (SteeleInteraction): Steele interaction configuration.
    """

    def __init__(
        self,
        coulomb: Coulomb | None = None,
        lennard_jones: LennardJones | None = None,
        fumi_tosi: TosiFumi | None = None,
        xft: XFT | None = None,
        daim: DAIM | None = None,
        damping: Damping | None = None,
        steele: Steele | None = None,
    ):
        self.coulomb = coulomb
        self.lennard_jones = lennard_jones
        self.fumi_tosi = fumi_tosi
        self.xft = xft
        self.daim = daim
        self.damping = damping
        self.steele = steele

    def as_dict(self):
        return {
            "coulomb": self.coulomb.as_dict(),
            "lennard_jones": self.lennard_jones.as_dict(),
            "fumi_tosi": self.fumi_tosi.as_dict(),
            "xft": self.xft.as_dict(),
            "daim": self.daim.as_dict(),
            "damping": self.damping.as_dict(),
            "steele": self.steele.as_dict(),
        }

    @classmethod
    def from_dict(cls, d):
        coulomb = Coulomb.from_dict(d["coulomb"])
        lennard_jones = LennardJones.from_dict(d["lennard_jones"])
        fumi_tosi = TosiFumi.from_dict(d["fumi_tosi"])
        xft = XFT.from_dict(d["xft"])
        daim = DAIM.from_dict(d["daim"])
        damping = Damping.from_dict(d["damping"])
        steele = Steele.from_dict(d["steele"])

        return cls(
            coulomb,
            lennard_jones,
            fumi_tosi,
            xft,
            daim,
            damping,
            steele,
        )


class Coulomb(MSONable):
    """Coulomb class representing the coulomb block in the simulation.

    Attributes:
        coulomb_rcut (float): Cutoff distance for real space interactions. Default is 0.0.
        coulomb_rtol (float): Magnitude below which real-space term is not included in the summation.
                              Default is 1.0e-15.
        coulomb_ktol (float): Magnitude below which reciprocal-space term is not included in the
                              summation. Default is 1.0e-15.
    """

    def __init__(
        self,
        coulomb_rcut: float = 0.0,
        coulomb_rtol: float = 1.0e-15,
        coulomb_ktol: float = 1.0e-15,
    ):
        self.coulomb_rcut = coulomb_rcut
        self.coulomb_rtol = coulomb_rtol
        self.coulomb_ktol = coulomb_ktol

    def as_dict(self):
        return {
            "coulomb_rcut": self.coulomb_rcut,
            "coulomb_rtol": self.coulomb_rtol,
            "coulomb_ktol": self.coulomb_ktol,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            coulomb_rcut=d.get("coulomb_rcut", 0.0),
            coulomb_rtol=d.get("coulomb_rtol", 1.0e-15),
            coulomb_ktol=d.get("coulomb_ktol", 1.0e-15),
        )


class LennardJones(MSONable):
    """LennardJones class representing the lennard-jones block in the simulation.

    Attributes:
        lj_rcut (float): Cutoff distance for all Lennard-Jones interactions in atomic units (a0).
        lj_3D_tail_correction (bool): Optional flag to enable long-range tail correction in 3D fluid system.
                                     Default is False.
        lj_pairs (list[tuple[str, str, float, float]]): List of Lennard-Jones pair parameters.
        lj_rule (str | None): Optional mixing rule to calculate cross Lennard-Jones parameters.
                              Default is None.

    Note: The lj_pairs is a list of tuples, where each tuple represents a Lennard-Jones interaction pair.
          The tuple contains the following elements: (species_name1, species_name2, epsilon, sigma).
          Epsilon (ϵ) is in kJ/mol, and Sigma (o) is in angstroms, not in atomic units.
          If lj_rule is provided, cross Lennard-Jones parameters will be calculated using mixing rules.

    """

    def __init__(
        self,
        lj_rcut: float,
        lj_3D_tail_correction: bool = False,
        lj_pairs: list[tuple[str, str, float, float]] | None = None,
        lj_rule: str | None = None,
    ):
        self.lj_rcut = lj_rcut
        self.lj_3D_tail_correction = lj_3D_tail_correction
        self.lj_pairs = lj_pairs or []
        self.lj_rule = lj_rule

    def as_dict(self):
        return {
            "lj_rcut": self.lj_rcut,
            "lj_3D_tail_correction": self.lj_3D_tail_correction,
            "lj_pairs": self.lj_pairs,
            "lj_rule": self.lj_rule,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            lj_rcut=d["lj_rcut"],
            lj_3D_tail_correction=d.get("lj_3D_tail_correction", False),
            lj_pairs=d.get("lj_pairs", []),
            lj_rule=d.get("lj_rule"),
        )


class TosiFumi(MSONable):
    r"""TosiFumi class representing the fumi-tosi block in the simulation.

    Attributes:
        ft_rcut (float): Cutoff distance for all Fumi-Tosi interactions in atomic units (\(a_0\)).
        ft_3D_pressure_tail_correction (bool): Optional flag to enable long-range tail correction in 3D fluid system.
                                               Default is False.
        ft_pairs (list[tuple[str, str, float, float, float, float, float, float, float, float]]): List of Fumi-Tosi
                                                                                                  pair parameters.
        ft_rule (str | None): Optional mixing rule to calculate cross Fumi-Tosi parameters.
                              Default is None.

    Note:
        The `ft_pairs` is a list of tuples, where each tuple represents a Fumi-Tosi interaction pair.
        The tuple contains the following elements: (species_name1, species_name2, \(\eta_{ij}\), \(B_{ij}\),
                                                   \(C_{ij}\), \(D_{ij}\), \(d_{dd}^{ij}\), \(d_{dq}^{ij}\)).
        All parameters are in atomic units (\(a_0\)).
        If `ft_rule` is provided, cross Fumi-Tosi parameters will be calculated using mixing rules.

    """

    def __init__(
        self,
        ft_rcut: float,
        ft_3D_pressure_tail_correction: bool = False,
        ft_pairs: list[tuple[str, str, float, float, float, float, float, float, float, float]] | None = None,
        ft_rule: str | None = None,
    ):
        self.ft_rcut = ft_rcut
        self.ft_3D_pressure_tail_correction = ft_3D_pressure_tail_correction
        self.ft_pairs = ft_pairs or []
        self.ft_rule = ft_rule

    def as_dict(self):
        return {
            "ft_rcut": self.ft_rcut,
            "ft_3D_pressure_tail_correction": self.ft_3D_pressure_tail_correction,
            "ft_pairs": self.ft_pairs,
            "ft_rule": self.ft_rule,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            ft_rcut=d["ft_rcut"],
            ft_3D_pressure_tail_correction=d.get("ft_3D_pressure_tail_correction", False),
            ft_pairs=d.get("ft_pairs", []),
            ft_rule=d.get("ft_rule"),
        )


class XFT(MSONable):
    r"""XFT class representing the xft block in the simulation.

    Attributes:
        xft_rcut (float): Cutoff distance for all XFT interactions in atomic units (\(a_0\)).
        xft_3D_pressure_tail_correction (bool): Optional flag to enable long-range tail correction in 3D fluid system.
                                                Default is False.
        xft_pairs (list[tuple[str, str, float, float, int, float, float, float, float, float, float]]):
                                                List of XFT pair parameters.
        xft_rule (str | None): Optional mixing rule to calculate cross XFT parameters.
                               Default is None.

    Note:
        The `xft_pairs` is a list of tuples, where each tuple represents an XFT interaction pair.
        The tuple contains the following elements: (species_name1, species_name2, \(\eta_{ij}\), \(B_{ij}\), \(n\),
            \(\eta'_{ij}\), \(B'_{ij}\), \(C_{ij}\), \(D_{ij}\), \(d_{dd}^{ij}\), \(d_{dq}^{ij}\)).
        All parameters are in atomic units (\(a_0\)).
        If `xft_rule` is provided, cross XFT parameters will be calculated using mixing rules.

    """

    def __init__(
        self,
        xft_rcut: float,
        xft_3D_pressure_tail_correction: bool = False,
        xft_pairs: list[tuple[str, str, float, float, int, float, float, float, float, float, float]] | None = None,
        xft_rule: str | None = None,
    ):
        self.xft_rcut = xft_rcut
        self.xft_3D_pressure_tail_correction = xft_3D_pressure_tail_correction
        self.xft_pairs = xft_pairs or []
        self.xft_rule = xft_rule

    def as_dict(self):
        return {
            "xft_rcut": self.xft_rcut,
            "xft_3D_pressure_tail_correction": self.xft_3D_pressure_tail_correction,
            "xft_pairs": self.xft_pairs,
            "xft_rule": self.xft_rule,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            xft_rcut=d["xft_rcut"],
            xft_3D_pressure_tail_correction=d.get("xft_3D_pressure_tail_correction", False),
            xft_pairs=d.get("xft_pairs", []),
            xft_rule=d.get("xft_rule"),
        )


class DAIM(MSONable):
    r"""DAIM class representing the daim block in the simulation.

    Attributes:
        daim_rcut (float): Cutoff distance for all DAIM interactions in atomic units (\(a_0\)).
        daim_3D_pressure_tail_correction (bool): Optional flag to enable long-range tail correction in 3D fluid system.
                                                 Default is False.
        daim_pairs (list[tuple[str, str, float, float, float, float, float, float, float, float, float, float, float]]):
            List of DAIM pair parameters.

    Note:
        The `daim_pairs` is a list of tuples, where each tuple represents a DAIM interaction pair.
        The tuple contains the following elements: (species_name1, species_name2, \(\eta^1_{ij}\), \(\eta^2_{ij}\),
            \(\eta^3_{ij}\), \(B^1_{ij}\), \(B^2_{ij}\), \(B^3_{ij}\), \(C_{ij}\), \(D_{ij}\), \(d_{dd}^{ij}\),
            \(d_{dq}^{ij}\)).
        All parameters are in atomic units (\(a_0\)).

    """

    def __init__(
        self,
        daim_rcut: float,
        daim_3D_pressure_tail_correction: bool = False,
        daim_pairs: list[tuple[str, str, float, float, float, float, float, float, float, float, float, float, float]]
        | None = None,
    ):
        self.daim_rcut = daim_rcut
        self.daim_3D_pressure_tail_correction = daim_3D_pressure_tail_correction
        self.daim_pairs = daim_pairs or []

    def as_dict(self):
        return {
            "daim_rcut": self.daim_rcut,
            "daim_3D_pressure_tail_correction": self.daim_3D_pressure_tail_correction,
            "daim_pairs": self.daim_pairs,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            daim_rcut=d["daim_rcut"],
            daim_3D_pressure_tail_correction=d.get("daim_3D_pressure_tail_correction", False),
            daim_pairs=d.get("daim_pairs", []),
        )


class Damping(MSONable):
    r"""Damping class representing the damping block in the simulation.

    Attributes:
        tt_pairs (list[tuple[str, str, float, int, float]]): List of Tang-Toennies damping function parameters.

    Note:
        The `tt_pairs` is a list of tuples, where each tuple represents a pair of species for the damping function.
        The tuple contains the following elements: (species_name_charge, species_name_dipole, bn, n, cn).

    The Tang-Toennies damping function is defined as:

    .. math::

        f_n(r_{ij}) = 1 - c_n \\exp(-b_n r_{ij}) \\sum_{k=0}^{n} \frac{(b_n r_{ij})^k}{k!}

    where:
        - :math:`f_n(r_{ij})`: Damping function for charge-dipole interactions.
        - :math:`r_{ij}`: Distance between the interacting species.
        - :math:`b_n`: Parameter controlling the damping range.
        - :math:`n`: Exponent of the damping function.
        - :math:`c_n`: Coefficient of the damping function.

    """

    def __init__(
        self,
        tt_pairs: list[tuple[str, str, float, int, float]] | None = None,
    ):
        self.tt_pairs = tt_pairs or []

    def as_dict(self):
        return {
            "tt_pairs": self.tt_pairs,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            tt_pairs=d.get("tt_pairs", []),
        )


class Steele(MSONable):
    r"""Steele class representing the Steele potential for non-structured walls.

    Attributes:
        num_walls (int): The number of different walls (typically 1 or 2).
        steele_rcut (float): The cut-off distance for all Steele interactions in atomic units (a0).

        walls (list[tuple[str, str, float, float, float, float]]): List of wall parameters.
            Each tuple contains the following elements: (electrode_name, species_name,
            rho_w, epsilon_fw, sigma_fw, delta).

    The Steele potential is a one-dimensional potential used to model non-structured walls. The potential is defined as:

    .. math::

        V(|z_i - z_0|) = 2 \\pi \rho_w \\epsilon_{fw} \\sigma_{fw} \\Delta \\left[
            \frac{2}{5} \\left( \frac{\\sigma_{fw}}{z} \right)^{10} - \\left( \frac{\\sigma_{fw}}{z} \right)^4 -
            \frac{\\sigma_{fw}^4}{3 \\Delta (z + 0.61 \\Delta)^3} \right]

    where:
        - :math:`V(|z_i - z_0|)`: Steele potential for non-structured walls.
        - :math:`z_i`: Distance to the wall.
        - :math:`z_0`: Wall position.
        - :math:`\rho_w`: Wall density in angstroms^-3.
        - :math:`\\epsilon_{fw}`: Energy parameter in kJ/mol.
        - :math:`\\sigma_{fw}`: Length parameter in angstroms.
        - :math:`\\Delta`: Parameter in angstroms.

    """

    def __init__(
        self,
        num_walls: int,
        steele_rcut: float,
        walls: list[tuple[str, str, float, float, float, float]] | None = None,
    ):
        self.num_walls = num_walls
        self.steele_rcut = steele_rcut
        self.walls = walls or []

    def as_dict(self):
        return {
            "num_walls": self.num_walls,
            "steele_rcut": self.steele_rcut,
            "walls": self.walls,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            num_walls=d.get("num_walls", 0),
            steele_rcut=d.get("steele_rcut", 0.0),
            walls=d.get("walls", []),
        )
