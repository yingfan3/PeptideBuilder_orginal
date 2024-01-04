"""This module is part of the PeptideBuilder library,
written by Matthew Z. Tien, Dariya K. Sydykova,
Austin G. Meyer, and Claus O. Wilke.

The PeptideBuilder module contains code to generate 3D
structures of peptides. It requires the Geometry module
(also part of the PeptideBuilder library), which contains
default bond lengths and angles for all amino acids.

This module also requires the Bio.PDB module from
Biopython, for structure manipulation.

This file is provided to you under the MIT License."""

import math, warnings
from typing import List, Optional, Union

from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.vectors import Vector, rotaxis, calc_dihedral, calc_angle
import numpy as np

from .Geometry import (
    AlaGeo,
    ArgGeo,
    AsnGeo,
    AspGeo,
    CysGeo,
    GlnGeo,
    GluGeo,
    GlyGeo,
    HisGeo,
    IleGeo,
    LeuGeo,
    LysGeo,
    MetGeo,
    PheGeo,
    ProGeo,
    SerGeo,
    ThrGeo,
    TrpGeo,
    TyrGeo,
    ValGeo,
    geometry,
    Geo,
)


def calculateCoordinates(
        refA: Residue, refB: Residue, refC: Residue, L: float, ang: float, di: float
) -> np.ndarray:
    AV = refA.get_vector()
    BV = refB.get_vector()
    CV = refC.get_vector()

    CA = AV - CV
    CB = BV - CV

    ##CA vector
    AX = CA[0]
    AY = CA[1]
    AZ = CA[2]

    ##CB vector
    BX = CB[0]
    BY = CB[1]
    BZ = CB[2]

    ##Plane Parameters
    A = (AY * BZ) - (AZ * BY)
    B = (AZ * BX) - (AX * BZ)
    G = (AX * BY) - (AY * BX)

    ##Dot Product Constant
    F = math.sqrt(BX * BX + BY * BY + BZ * BZ) * L * math.cos(ang * (math.pi / 180.0))

    ##Constants
    const = math.sqrt(
        math.pow((B * BZ - BY * G), 2)
        * (
                -(F * F) * (A * A + B * B + G * G)
                + (
                        B * B * (BX * BX + BZ * BZ)
                        + A * A * (BY * BY + BZ * BZ)
                        - (2 * A * BX * BZ * G)
                        + (BX * BX + BY * BY) * G * G
                        - (2 * B * BY) * (A * BX + BZ * G)
                )
                * L
                * L
        )
    )
    denom = (
            (B * B) * (BX * BX + BZ * BZ)
            + (A * A) * (BY * BY + BZ * BZ)
            - (2 * A * BX * BZ * G)
            + (BX * BX + BY * BY) * (G * G)
            - (2 * B * BY) * (A * BX + BZ * G)
    )

    X = (
                (B * B * BX * F) - (A * B * BY * F) + (F * G) * (-A * BZ + BX * G) + const
        ) / denom

    if (B == 0 or BZ == 0) and (BY == 0 or G == 0):
        const1 = math.sqrt(
            G * G * (-A * A * X * X + (B * B + G * G) * (L - X) * (L + X))
        )
        Y = ((-A * B * X) + const1) / (B * B + G * G)
        Z = -(A * G * G * X + B * const1) / (G * (B * B + G * G))
    else:
        Y = (
                    (A * A * BY * F) * (B * BZ - BY * G)
                    + G * (-F * math.pow(B * BZ - BY * G, 2) + BX * const)
                    - A * (B * B * BX * BZ * F - B * BX * BY * F * G + BZ * const)
            ) / ((B * BZ - BY * G) * denom)
        Z = (
                    (A * A * BZ * F) * (B * BZ - BY * G)
                    + (B * F) * math.pow(B * BZ - BY * G, 2)
                    + (A * BX * F * G) * (-B * BZ + BY * G)
                    - B * BX * const
                    + A * BY * const
            ) / ((B * BZ - BY * G) * denom)

    # Get the new Vector from the origin
    D = Vector(X, Y, Z) + CV
    with warnings.catch_warnings():
        # ignore inconsequential warning
        warnings.simplefilter("ignore")
        temp = calc_dihedral(AV, BV, CV, D) * (180.0 / math.pi)

    di = di - temp
    rot = rotaxis(math.pi * (di / 180.0), CV - BV)
    D = (D - BV).left_multiply(rot) + BV

    return D.get_array()


def Gly(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: GlyGeo):
    SG = []
    CD1 = []

    return SG, CD1


def Ala(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: AlaGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")

    SG = [CB]
    CD1 = [CI]

    return SG, CD1


def Ser(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: SerGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    oxygen_g = calculateCoordinates(
        NB, SG, CB, geo.CB_OG_length, geo.SG_CB_OG_angle, geo.NB_SG_CB_OG_diangle
    )
    OG = Atom("OG", oxygen_g, 0.0, 1.0, " ", " OG", 0, "O")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    oxygen_j = calculateCoordinates(
        N, CD1, CI, geo.CI_OJ_length, geo.CD1_CI_OJ_angle, geo.N_CD1_CI_OJ_diangle
    )
    OJ = Atom("OJ", oxygen_j, 0.0, 1.0, " ", " OJ", 0, "O")

    SG = [CB, OG]
    CD1 = [CI, OJ]

    return SG, CD1


def Cys(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: CysGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    sulfur_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_SG2_length, geo.SG_CB_SG2_angle, geo.NB_SG_CB_SG2_diangle
    )
    SG2 = Atom("SG2", sulfur_g2, 0.0, 1.0, " ", " SG2", 0, "S")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    sulfur_g1 = calculateCoordinates(
        N, CD1, CI, geo.CI_SG1_length, geo.CD1_CI_SG1_angle, geo.N_CD1_CI_SG1_diangle
    )
    SG1 = Atom("SG1", sulfur_g1, 0.0, 1.0, " ", " SG1", 0, "S")

    SG = [CB, SG2]
    CD1 = [CI, SG1]

    return SG, CD1


def Val(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: ValGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g1 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG1_length, geo.SG_CB_CG1_angle, geo.NB_SG_CB_CG1_diangle
    )
    CG1 = Atom("CG1", carbon_g1, 0.0, 1.0, " ", " CG1", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j1 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ1_length, geo.CD1_CI_CJ1_angle, geo.N_CD1_CI_CJ1_diangle
    )
    CJ1 = Atom("CJ1", carbon_j1, 0.0, 1.0, " ", " CJ1", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")

    SG = [CB, CG1, CG2]
    CD1 = [CI, CJ1, CJ2]

    return SG, CD1


def Ile(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: IleGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g1 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG1_length, geo.SG_CB_CG1_angle, geo.NB_SG_CB_CG1_diangle
    )
    CG1 = Atom("CG1", carbon_g1, 0.0, 1.0, " ", " CG1", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    carbon_d2 = calculateCoordinates(
        SG, CB, CG1, geo.CG1_CD2_length, geo.CB_CG1_CD2_angle, geo.SG_CB_CG1_CD2_diangle
    )
    CD2 = Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j1 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ1_length, geo.CD1_CI_CJ1_angle, geo.N_CD1_CI_CJ1_diangle
    )
    CJ1 = Atom("CJ1", carbon_j1, 0.0, 1.0, " ", " CJ1", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")
    carbon_k2 = calculateCoordinates(
        CD1, CI, CJ1, geo.CJ1_CK2_length, geo.CI_CJ1_CK2_angle, geo.CD1_CI_CJ1_CK2_diangle
    )
    CK2 = Atom("CK2", carbon_k2, 0.0, 1.0, " ", " CK2", 0, "C")

    SG = [CB, CG1, CG2, CD2]
    CD1 = [CI, CJ1, CJ2, CK2]

    return SG, CD1


def Leu(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: LeuGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    carbon_d2 = calculateCoordinates(
        SG, CB, CG2, geo.CG2_CD2_length, geo.CB_CG2_CD2_angle, geo.SG_CB_CG2_CD2_diangle
    )
    CD2 = Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    carbon_d3 = calculateCoordinates(
        SG, CB, CG2, geo.CG2_CD3_length, geo.CB_CG2_CD3_angle, geo.SG_CB_CG2_CD3_diangle
    )
    CD3 = Atom("CD3", carbon_d3, 0.0, 1.0, " ", " CD3", 0, "C")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")
    carbon_k2 = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_CK2_length, geo.CI_CJ2_CK2_angle, geo.CD1_CI_CJ2_CK2_diangle
    )
    CK2 = Atom("CK2", carbon_k2, 0.0, 1.0, " ", " CK2", 0, "C")
    carbon_k3 = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_CK3_length, geo.CI_CJ2_CK3_angle, geo.CD1_CI_CJ2_CK3_diangle
    )
    CK3 = Atom("CK3", carbon_k3, 0.0, 1.0, " ", " CK3", 0, "C")

    SG = [CB, CG2, CD2, CD3]
    CD1 = [CI, CJ2, CK2, CK3]

    return SG, CD1


def Thr(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: ThrGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    oxygen_g1 = calculateCoordinates(
        NB, SG, CB, geo.CB_OG1_length, geo.SG_CB_OG1_angle, geo.NB_SG_CB_OG1_diangle
    )
    OG1 = Atom("OG1", oxygen_g1, 0.0, 1.0, " ", " OG1", 0, "O")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    oxygen_j1 = calculateCoordinates(
        N, CD1, CI, geo.CI_OJ1_length, geo.CD1_CI_OJ1_angle, geo.N_CD1_CI_OJ1_diangle
    )
    OJ1 = Atom("OJ1", oxygen_j1, 0.0, 1.0, " ", " OJ1", 0, "O")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")

    SG = [CB, OG1, CG2]
    CD1 = [CI, OJ1, CJ2]

    return SG, CD1


def Arg(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: ArgGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    carbon_d2 = calculateCoordinates(
        SG, CB, CG2, geo.CG2_CD2_length, geo.CB_CG2_CD2_angle, geo.SG_CB_CG2_CD2_diangle
    )
    CD2 = Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    nitrogen_e = calculateCoordinates(
        CB, CG2, CD2, geo.CD2_NE_length, geo.CG2_CD2_NE_angle, geo.CB_CG2_CD2_NE_diangle
    )
    NE = Atom("NE", nitrogen_e, 0.0, 1.0, " ", " NE", 0, "N")
    carbon_z = calculateCoordinates(
        CG2, CD2, NE, geo.NE_CZ_length, geo.CD2_NE_CZ_angle, geo.CG2_CD2_NE_CZ_diangle
    )
    CZ = Atom("CZ", carbon_z, 0.0, 1.0, " ", " CZ", 0, "C")
    nitrogen_h1 = calculateCoordinates(
        CD2, NE, CZ, geo.CZ_NH1_length, geo.NE_CZ_NH1_angle, geo.CD2_NE_CZ_NH1_diangle
    )
    NH1 = Atom("NH1", nitrogen_h1, 0.0, 1.0, " ", " NH1", 0, "N")
    nitrogen_h2 = calculateCoordinates(
        CD2, NE, CZ, geo.CZ_NH2_length, geo.NE_CZ_NH2_angle, geo.CD2_NE_CZ_NH2_diangle
    )
    NH2 = Atom("NH2", nitrogen_h2, 0.0, 1.0, " ", " NH2", 0, "N")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")
    carbon_k2 = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_CK2_length, geo.CI_CJ2_CK2_angle, geo.CD1_CI_CJ2_CK2_diangle
    )
    CK2 = Atom("CK2", carbon_k2, 0.0, 1.0, " ", " CK2", 0, "C")
    nitrogen_m = calculateCoordinates(
        CI, CJ2, CK2, geo.CK2_NM_length, geo.CJ2_CK2_NM_angle, geo.CI_CJ2_CK2_NM_diangle
    )
    NM = Atom("NM", nitrogen_m, 0.0, 1.0, " ", " NM", 0, "N")
    carbon_n = calculateCoordinates(
        CJ2, CK2, NM, geo.NM_CN_length, geo.CK2_NM_CN_angle, geo.CJ2_CK2_NM_CN_diangle
    )
    CN = Atom("CN", carbon_n, 0.0, 1.0, " ", " CN", 0, "C")
    nitrogen_q1 = calculateCoordinates(
        CK2, NM, CN, geo.CN_NQ1_length, geo.NM_CN_NQ1_angle, geo.CK2_NM_CN_NQ1_diangle
    )
    NQ1 = Atom("NQ1", nitrogen_q1, 0.0, 1.0, " ", " NQ1", 0, "N")
    nitrogen_q2 = calculateCoordinates(
        CK2, NM, CN, geo.CN_NQ2_length, geo.NM_CN_NQ2_angle, geo.CK2_NM_CN_NQ2_diangle
    )
    NQ2 = Atom("NQ2", nitrogen_q2, 0.0, 1.0, " ", " NQ2", 0, "N")

    SG = [CB, CG2, CD2, NE, CZ, NH1, NH2]
    CD1 = [CI, CJ2, CK2, NM, CN, NQ1, NQ2]

    return SG, CD1


def Lys(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: LysGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    carbon_d = calculateCoordinates(
        SG, CB, CG2, geo.CG2_CD_length, geo.CB_CG2_CD_angle, geo.SG_CB_CG2_CD_diangle
    )
    CD = Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    carbon_e = calculateCoordinates(
        CB, CG2, CD, geo.CD_CE_length, geo.CG2_CD_CE_angle, geo.CB_CG2_CD_CE_diangle
    )
    CE = Atom("CE", carbon_e, 0.0, 1.0, " ", " CE", 0, "C")
    nitrogen_z = calculateCoordinates(
        CG2, CD, CE, geo.CE_NZ_length, geo.CD_CE_NZ_angle, geo.CG2_CD_CE_NZ_diangle
    )
    NZ = Atom("NZ", nitrogen_z, 0.0, 1.0, " ", " NZ", 0, "N")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")
    carbon_k = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_CK_length, geo.CI_CJ2_CK_angle, geo.CD1_CI_CJ2_CK_diangle
    )
    CK = Atom("CK", carbon_k, 0.0, 1.0, " ", " CK", 0, "C")
    carbon_m = calculateCoordinates(
        CI, CJ2, CK, geo.CK_CM_length, geo.CJ2_CK_CM_angle, geo.CI_CJ2_CK_CM_diangle
    )
    CM = Atom("CM", carbon_m, 0.0, 1.0, " ", " CM", 0, "C")
    nitrogen_z2 = calculateCoordinates(
        CJ2, CK, CM, geo.CM_NZ2_length, geo.CK_CM_NZ2_angle, geo.CJ2_CK_CM_NZ2_diangle
    )
    NZ2 = Atom("NZ2", nitrogen_z2, 0.0, 1.0, " ", " NZ2", 0, "N")

    SG = [CB, CG2, CD, CE, NZ]
    CD1 = [CI, CJ2, CK, CM, NZ2]

    return SG, CD1


def Asp(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: AspGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    oxygen_d3 = calculateCoordinates(
        SG, CB, CG2, geo.CG2_OD3_length, geo.CB_CG2_OD3_angle, geo.SG_CB_CG2_OD3_diangle
    )
    OD3 = Atom("OD3", oxygen_d3, 0.0, 1.0, " ", " OD3", 0, "O")
    oxygen_d4 = calculateCoordinates(
        SG, CB, CG2, geo.CG2_OD4_length, geo.CB_CG2_OD4_angle, geo.SG_CB_CG2_OD4_diangle
    )
    OD4 = Atom("OD4", oxygen_d4, 0.0, 1.0, " ", " OD4", 0, "O")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")
    oxygen_k3 = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_OK3_length, geo.CI_CJ2_OK3_angle, geo.CD1_CI_CJ2_OK3_diangle
    )
    OK3 = Atom("OK3", oxygen_k3, 0.0, 1.0, " ", " OK3", 0, "O")
    oxygen_k4 = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_OK4_length, geo.CI_CJ2_OK4_angle, geo.CD1_CI_CJ2_OK4_diangle
    )
    OK4 = Atom("OK4", oxygen_k4, 0.0, 1.0, " ", " OK4", 0, "O")

    SG = [CB, CG2, OD3, OD4]
    CD1 = [CI, CJ2, OK3, OK4]

    return SG, CD1


def Asn(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: AsnGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    oxygen_d3 = calculateCoordinates(
        SG, CB, CG2, geo.CG2_OD3_length, geo.CB_CG2_OD3_angle, geo.SG_CB_CG2_OD3_diangle
    )
    OD3 = Atom("OD3", oxygen_d3, 0.0, 1.0, " ", " OD3", 0, "O")
    nitrogen_d2 = calculateCoordinates(
        SG, CB, CG2, geo.CG2_ND2_length, geo.CB_CG2_ND2_angle, geo.SG_CB_CG2_ND2_diangle
    )
    ND2 = Atom("ND2", nitrogen_d2, 0.0, 1.0, " ", " ND2", 0, "N")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")
    oxygen_k3 = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_OK3_length, geo.CI_CJ2_OK3_angle, geo.CD1_CI_CJ2_OK3_diangle
    )
    OK3 = Atom("OK3", oxygen_k3, 0.0, 1.0, " ", " OK3", 0, "O")
    nitrogen_k2 = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_NK2_length, geo.CI_CJ2_NK2_angle, geo.CD1_CI_CJ2_NK2_diangle
    )
    NK2 = Atom("NK2", nitrogen_k2, 0.0, 1.0, " ", " NK2", 0, "N")
    SG = [CB, CG2, OD3, ND2]
    CD1 = [CI, CJ2, OK3, NK2]

    return SG, CD1


def Glu(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: GluGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    carbon_d2 = calculateCoordinates(
        SG, CB, CG2, geo.CG2_CD2_length, geo.CB_CG2_CD2_angle, geo.SG_CB_CG2_CD2_diangle
    )
    CD2 = Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    oxygen_e1 = calculateCoordinates(
        CB, CG2, CD2, geo.CD2_OE1_length, geo.CG2_CD2_OE1_angle, geo.CB_CG2_CD2_OE1_diangle
    )
    OE1 = Atom("OE1", oxygen_e1, 0.0, 1.0, " ", " OE1", 0, "O")
    oxygen_e2 = calculateCoordinates(
        CB, CG2, CD2, geo.CD2_OE2_length, geo.CG2_CD2_OE2_angle, geo.CB_CG2_CD2_OE2_diangle
    )
    OE2 = Atom("OE2", oxygen_e2, 0.0, 1.0, " ", " OE2", 0, "O")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")
    carbon_k2 = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_CK2_length, geo.CI_CJ2_CK2_angle, geo.CD1_CI_CJ2_CK2_diangle
    )
    CK2 = Atom("CK2", carbon_k2, 0.0, 1.0, " ", " CK2", 0, "C")
    oxygen_m1 = calculateCoordinates(
        CI, CJ2, CK2, geo.CK2_OM1_length, geo.CJ2_CK2_OM1_angle, geo.CI_CJ2_CK2_OM1_diangle
    )
    OM1 = Atom("OM1", oxygen_m1, 0.0, 1.0, " ", " OM1", 0, "O")
    oxygen_m2 = calculateCoordinates(
        CI, CJ2, CK2, geo.CK2_OM2_length, geo.CJ2_CK2_OM2_angle, geo.CI_CJ2_CK2_OM2_diangle
    )
    OM2 = Atom("OM2", oxygen_m2, 0.0, 1.0, " ", " OM2", 0, "O")
    SG = [CB, CG2, CD2, OE1, OE2]
    CD1 = [CI, CJ2, CK2, OM1, OM2]

    return SG, CD1


def Gln(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: GlnGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    carbon_d = calculateCoordinates(
        SG, CB, CG2, geo.CG2_CD_length, geo.CB_CG2_CD_angle, geo.SG_CB_CG2_CD_diangle
    )
    CD = Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    oxygen_e1 = calculateCoordinates(
        CB, CG2, CD, geo.CD_OE1_length, geo.CG2_CD_OE1_angle, geo.CB_CG2_CD_OE1_diangle
    )
    OE1 = Atom("OE1", oxygen_e1, 0.0, 1.0, " ", " OE1", 0, "O")
    nitrogen_e2 = calculateCoordinates(
        CB, CG2, CD, geo.CD_NE2_length, geo.CG2_CD_NE2_angle, geo.CB_CG2_CD_NE2_diangle
    )
    NE2 = Atom("NE2", nitrogen_e2, 0.0, 1.0, " ", " NE2", 0, "N")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")
    carbon_k = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_CK_length, geo.CI_CJ2_CK_angle, geo.CD1_CI_CJ2_CK_diangle
    )
    CK = Atom("CK", carbon_k, 0.0, 1.0, " ", " CK", 0, "C")
    oxygen_m1 = calculateCoordinates(
        CI, CJ2, CK, geo.CK_OM1_length, geo.CJ2_CK_OM1_angle, geo.CI_CJ2_CK_OM1_diangle
    )
    OM1 = Atom("OM1", oxygen_m1, 0.0, 1.0, " ", " OM1", 0, "O")
    nitrogen_m2 = calculateCoordinates(
        CI, CJ2, CK, geo.CK_NM2_length, geo.CJ2_CK_NM2_angle, geo.CI_CJ2_CK_NM2_diangle
    )
    NM2 = Atom("NM2", nitrogen_m2, 0.0, 1.0, " ", " NM2", 0, "N")
    SG = [CB, CG2, CD, OE1, NE2]
    CD1 = [CI, CJ2, CK, OM1, NM2]

    return SG, CD1


def Met(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: MetGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    sulfur_d = calculateCoordinates(
        SG, CB, CG2, geo.CG2_SD_length, geo.CB_CG2_SD_angle, geo.SG_CB_CG2_SD_diangle
    )
    SD = Atom("SD", sulfur_d, 0.0, 1.0, " ", " SD", 0, "S")
    carbon_e = calculateCoordinates(
        CB, CG2, SD, geo.SD_CE_length, geo.CG2_SD_CE_angle, geo.CB_CG2_SD_CE_diangle
    )
    CE = Atom("CE", carbon_e, 0.0, 1.0, " ", " CE", 0, "C")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")
    sulfur_k = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_SK_length, geo.CI_CJ2_SK_angle, geo.CD1_CI_CJ2_SK_diangle
    )
    SK = Atom("SK", sulfur_k, 0.0, 1.0, " ", " SK", 0, "S")
    carbon_m = calculateCoordinates(
        CI, CJ2, SK, geo.SK_CM_length, geo.CJ2_SK_CM_angle, geo.CI_CJ2_SK_CM_diangle
    )
    CM = Atom("CM", carbon_m, 0.0, 1.0, " ", " CM", 0, "C")
    SG = [CB, CG2, SD, CE]
    CD1 = [CI, CJ2, SK, CM]

    return SG, CD1


def His(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: HisGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    nitrogen_d1 = calculateCoordinates(
        SG, CB, CG2, geo.CG2_ND1_length, geo.CB_CG2_ND1_angle, geo.SG_CB_CG2_ND1_diangle
    )
    ND1 = Atom("ND1", nitrogen_d1, 0.0, 1.0, " ", " ND1", 0, "N")
    carbon_d2 = calculateCoordinates(
        SG, CB, CG2, geo.CG2_CD2_length, geo.CB_CG2_CD2_angle, geo.SG_CB_CG2_CD2_diangle
    )
    CD2 = Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    carbon_e2 = calculateCoordinates(
        CB, CG2, ND1, geo.ND1_CE2_length, geo.CG2_ND1_CE2_angle, geo.CB_CG2_ND1_CE2_diangle
    )
    CE2 = Atom("CE2", carbon_e2, 0.0, 1.0, " ", " CE2", 0, "C")
    nitrogen_e2 = calculateCoordinates(
        CB, CG2, CD2, geo.CD2_NE2_length, geo.CG2_CD2_NE2_angle, geo.CB_CG2_CD2_NE2_diangle
    )
    NE2 = Atom("NE2", nitrogen_e2, 0.0, 1.0, " ", " NE2", 0, "N")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")
    nitrogen_k1 = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_NK1_length, geo.CI_CJ2_NK1_angle, geo.CD1_CI_CJ2_NK1_diangle
    )
    NK1 = Atom("NK1", nitrogen_k1, 0.0, 1.0, " ", " NK1", 0, "N")
    carbon_k2 = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_CK2_length, geo.CI_CJ2_CK2_angle, geo.CD1_CI_CJ2_CK2_diangle
    )
    CK2 = Atom("CK2", carbon_k2, 0.0, 1.0, " ", " CK2", 0, "C")
    carbon_m2 = calculateCoordinates(
        CI, CJ2, NK1, geo.NK1_CM2_length, geo.CJ2_NK1_CM2_angle, geo.CI_CJ2_NK1_CM2_diangle
    )
    CM2 = Atom("CM2", carbon_m2, 0.0, 1.0, " ", " CM2", 0, "C")
    nitrogen_m2 = calculateCoordinates(
        CI, CJ2, CK2, geo.CK2_NM2_length, geo.CJ2_CK2_NM2_angle, geo.CI_CJ2_CK2_NM2_diangle
    )
    NM2 = Atom("NM2", nitrogen_m2, 0.0, 1.0, " ", " NM2", 0, "N")
    SG = [CB, CG2, ND1, CD2, CE2, NE2]
    CD1 = [CI, CJ2, NK1, CK2, CM2, NM2]

    return SG, CD1


def Pro(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: ProGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    carbon_d = calculateCoordinates(
        SG, CB, CG2, geo.CG2_CD_length, geo.CB_CG2_CD_angle, geo.SG_CB_CG2_CD_diangle
    )
    CD = Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")
    carbon_k = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_CK_length, geo.CI_CJ2_CK_angle, geo.CD1_CI_CJ2_CK_diangle
    )
    CK = Atom("CK", carbon_k, 0.0, 1.0, " ", " CK", 0, "C")
    SG = [CB, CG2, CD]
    CD1 = [CI, CJ2, CK]

    return SG, CD1


def Phe(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: PheGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    carbon_d2 = calculateCoordinates(
        SG, CB, CG2, geo.CG2_CD2_length, geo.CB_CG2_CD2_angle, geo.SG_CB_CG2_CD2_diangle
    )
    CD2 = Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    carbon_d3 = calculateCoordinates(
        SG, CB, CG2, geo.CG2_CD3_length, geo.CB_CG2_CD3_angle, geo.SG_CB_CG2_CD3_diangle
    )
    CD3 = Atom("CD3", carbon_d3, 0.0, 1.0, " ", " CD3", 0, "C")
    carbon_e2 = calculateCoordinates(
        CB, CG2, CD2, geo.CD2_CE2_length, geo.CG2_CD2_CE2_angle, geo.CB_CG2_CD2_CE2_diangle
    )
    CE2 = Atom("CE2", carbon_e2, 0.0, 1.0, " ", " CE2", 0, "C")
    carbon_e3 = calculateCoordinates(
        CB, CG2, CD3, geo.CD3_CE3_length, geo.CG2_CD3_CE3_angle, geo.CB_CG2_CD3_CE3_diangle
    )
    CE3 = Atom("CE3", carbon_e3, 0.0, 1.0, " ", " CE3", 0, "C")
    carbon_z = calculateCoordinates(
        CG2, CD2, CE2, geo.CE2_CZ_length, geo.CD2_CE2_CZ_angle, geo.CG2_CD2_CE2_CZ_diangle
    )
    CZ = Atom("CZ", carbon_z, 0.0, 1.0, " ", " CZ", 0, "C")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")
    carbon_k2 = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_CK2_length, geo.CI_CJ2_CK2_angle, geo.CD1_CI_CJ2_CK2_diangle
    )
    CK2 = Atom("CK2", carbon_k2, 0.0, 1.0, " ", " CK2", 0, "C")
    carbon_k3 = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_CK3_length, geo.CI_CJ2_CK3_angle, geo.CD1_CI_CJ2_CK3_diangle
    )
    CK3 = Atom("CK3", carbon_k3, 0.0, 1.0, " ", " CK3", 0, "C")
    carbon_m2 = calculateCoordinates(
        CI, CJ2, CK2, geo.CK2_CM2_length, geo.CJ2_CK2_CM2_angle, geo.CI_CJ2_CK2_CM2_diangle
    )
    CM2 = Atom("CM2", carbon_m2, 0.0, 1.0, " ", " CM2", 0, "C")
    carbon_m3 = calculateCoordinates(
        CI, CJ2, CK3, geo.CK3_CM3_length, geo.CJ2_CK3_CM3_angle, geo.CI_CJ2_CK3_CM3_diangle
    )
    CM3 = Atom("CM3", carbon_m3, 0.0, 1.0, " ", " CM3", 0, "C")
    carbon_n = calculateCoordinates(
        CJ2, CK2, CM2, geo.CM2_CN_length, geo.CK2_CM2_CN_angle, geo.CJ2_CK2_CM2_CN_diangle
    )
    CN = Atom("CN", carbon_n, 0.0, 1.0, " ", " CN", 0, "C")
    SG = [CB, CG2, CD2, CD3, CE2, CE3, CZ]
    CD1 = [CI, CJ2, CK2, CK3, CM2, CM3, CN]

    return SG, CD1


def Tyr(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1: TyrGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo1.SG_CB_length, geo1.NB_SG_CB_angle, geo1.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo1.CB_CG2_length, geo1.SG_CB_CG2_angle, geo1.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    carbon_d2 = calculateCoordinates(
        SG, CB, CG2, geo1.CG2_CD2_length, geo1.CB_CG2_CD2_angle, geo1.SG_CB_CG2_CD2_diangle
    )
    CD2 = Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    carbon_d3 = calculateCoordinates(
        SG, CB, CG2, geo1.CG2_CD3_length, geo1.CB_CG2_CD3_angle, geo1.SG_CB_CG2_CD3_diangle
    )
    CD3 = Atom("CD3", carbon_d3, 0.0, 1.0, " ", " CD3", 0, "C")
    carbon_e2 = calculateCoordinates(
        CB, CG2, CD2, geo1.CD2_CE2_length, geo1.CG2_CD2_CE2_angle, geo1.CB_CG2_CD2_CE2_diangle
    )
    CE2 = Atom("CE2", carbon_e2, 0.0, 1.0, " ", " CE2", 0, "C")
    carbon_e3 = calculateCoordinates(
        CB, CG2, CD3, geo1.CD3_CE3_length, geo1.CG2_CD3_CE3_angle, geo1.CB_CG2_CD3_CE3_diangle
    )
    CE3 = Atom("CE3", carbon_e3, 0.0, 1.0, " ", " CE3", 0, "C")
    carbon_z = calculateCoordinates(
        CG2, CD2, CE2, geo1.CE2_CZ_length, geo1.CD2_CE2_CZ_angle, geo1.CG2_CD2_CE2_CZ_diangle
    )
    CZ = Atom("CZ", carbon_z, 0.0, 1.0, " ", " CZ", 0, "C")
    oxygen_h = calculateCoordinates(
        CD2, CE2, CZ, geo1.CZ_OH_length, geo1.CE2_CZ_OH_angle, geo1.CD2_CE2_CZ_OH_diangle
    )
    OH = Atom("OH", oxygen_h, 0.0, 1.0, " ", " OH", 0, "O")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo1.CI_CD1_length, geo1.CI_CD1_CG_angle, geo1.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo1.CI_CJ2_length, geo1.CD1_CI_CJ2_angle, geo1.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")
    carbon_k2 = calculateCoordinates(
        CD1, CI, CJ2, geo1.CJ2_CK2_length, geo1.CI_CJ2_CK2_angle, geo1.CD1_CI_CJ2_CK2_diangle
    )
    CK2 = Atom("CK2", carbon_k2, 0.0, 1.0, " ", " CK2", 0, "C")
    carbon_k3 = calculateCoordinates(
        CD1, CI, CJ2, geo1.CJ2_CK3_length, geo1.CI_CJ2_CK3_angle, geo1.CD1_CI_CJ2_CK3_diangle
    )
    CK3 = Atom("CK3", carbon_k3, 0.0, 1.0, " ", " CK3", 0, "C")
    carbon_m2 = calculateCoordinates(
        CI, CJ2, CK2, geo1.CK2_CM2_length, geo1.CJ2_CK2_CM2_angle, geo1.CI_CJ2_CK2_CM2_diangle
    )
    CM2 = Atom("CM2", carbon_m2, 0.0, 1.0, " ", " CM2", 0, "C")
    carbon_m3 = calculateCoordinates(
        CI, CJ2, CK3, geo1.CK3_CM3_length, geo1.CJ2_CK3_CM3_angle, geo1.CI_CJ2_CK3_CM3_diangle
    )
    CM3 = Atom("CM3", carbon_m3, 0.0, 1.0, " ", " CM3", 0, "C")
    carbon_n = calculateCoordinates(
        CJ2, CK2, CM2, geo1.CM2_CN_length, geo1.CK2_CM2_CN_angle, geo1.CJ2_CK2_CM2_CN_diangle
    )
    CN = Atom("CN", carbon_n, 0.0, 1.0, " ", " CN", 0, "C")
    oxygen_q = calculateCoordinates(
        CK2, CM2, CN, geo1.CN_OQ_length, geo1.CM2_CN_OQ_angle, geo1.CK2_CM2_CN_OQ_diangle
    )
    OQ = Atom("OQ", oxygen_q, 0.0, 1.0, " ", " OQ", 0, "O")
    SG = [CB, CG2, CD2, CD3, CE2, CE3, CZ, OH]
    CD1 = [CI, CJ2, CK2, CK3, CM2, CM3, CN, OQ]

    return SG, CD1


def Trp(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: TrpGeo):
    carbon_b = calculateCoordinates(
        CG, NB, SG, geo.SG_CB_length, geo.NB_SG_CB_angle, geo.CG_NB_SG_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g2 = calculateCoordinates(
        NB, SG, CB, geo.CB_CG2_length, geo.SG_CB_CG2_angle, geo.NB_SG_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    carbon_d2 = calculateCoordinates(
        SG, CB, CG2, geo.CG2_CD2_length, geo.CB_CG2_CD2_angle, geo.SG_CB_CG2_CD2_diangle
    )
    CD2 = Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    carbon_d3 = calculateCoordinates(
        SG, CB, CG2, geo.CG2_CD3_length, geo.CB_CG2_CD3_angle, geo.SG_CB_CG2_CD3_diangle
    )
    CD3 = Atom("CD3", carbon_d3, 0.0, 1.0, " ", " CD3", 0, "C")
    nitrogen_e1 = calculateCoordinates(
        CB, CG2, CD2, geo.CD2_NE1_length, geo.CG2_CD2_NE1_angle, geo.CB_CG2_CD2_NE1_diangle
    )
    NE1 = Atom("NE1", nitrogen_e1, 0.0, 1.0, " ", " NE1", 0, "N")
    carbon_e2 = calculateCoordinates(
        CB, CG2, CD3, geo.CD3_CE2_length, geo.CG2_CD3_CE2_angle, geo.CB_CG2_CD3_CE2_diangle
    )
    CE2 = Atom("CE2", carbon_e2, 0.0, 1.0, " ", " CE2", 0, "C")
    carbon_e3 = calculateCoordinates(
        CB, CG2, CD3, geo.CD3_CE3_length, geo.CG2_CD3_CE3_angle, geo.CB_CG2_CD3_CE3_diangle
    )
    CE3 = Atom("CE3", carbon_e3, 0.0, 1.0, " ", " CE3", 0, "C")
    carbon_z2 = calculateCoordinates(
        CG2, CD3, CE2, geo.CE2_CZ2_length, geo.CD3_CE2_CZ2_angle, geo.CG2_CD3_CE2_CZ2_diangle
    )
    CZ2 = Atom("CZ2", carbon_z2, 0.0, 1.0, " ", " CZ2", 0, "C")
    carbon_z3 = calculateCoordinates(
        CG2, CD3, CE3, geo.CE3_CZ3_length, geo.CD3_CE3_CZ3_angle, geo.CG2_CD3_CE3_CZ3_diangle
    )
    CZ3 = Atom("CZ3", carbon_z3, 0.0, 1.0, " ", " CZ3", 0, "C")
    carbon_h2 = calculateCoordinates(
        CD3, CE2, CZ2, geo.CZ2_CH2_length, geo.CE2_CZ2_CH2_angle, geo.CD3_CE2_CZ2_CH2_diangle
    )
    CH2 = Atom("CH2", carbon_h2, 0.0, 1.0, " ", " CH2", 0, "C")

    carbon_i = calculateCoordinates(
        NB, CG, CD1, geo.CI_CD1_length, geo.CI_CD1_CG_angle, geo.CI_CD1_CG_NB_diangle
    )
    CI = Atom("CI", carbon_i, 0.0, 1.0, " ", " CI", 0, "C")
    carbon_j2 = calculateCoordinates(
        N, CD1, CI, geo.CI_CJ2_length, geo.CD1_CI_CJ2_angle, geo.N_CD1_CI_CJ2_diangle
    )
    CJ2 = Atom("CJ2", carbon_j2, 0.0, 1.0, " ", " CJ2", 0, "C")
    carbon_k2 = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_CK2_length, geo.CI_CJ2_CK2_angle, geo.CD1_CI_CJ2_CK2_diangle
    )
    CK2 = Atom("CK2", carbon_k2, 0.0, 1.0, " ", " CK2", 0, "C")
    carbon_k3 = calculateCoordinates(
        CD1, CI, CJ2, geo.CJ2_CK3_length, geo.CI_CJ2_CK3_angle, geo.CD1_CI_CJ2_CK3_diangle
    )
    CK3 = Atom("CK3", carbon_k3, 0.0, 1.0, " ", " CK3", 0, "C")
    nitrogen_m1 = calculateCoordinates(
        CI, CJ2, CK2, geo.CK2_NM1_length, geo.CJ2_CK2_NM1_angle, geo.CI_CJ2_CK2_NM1_diangle
    )
    NM1 = Atom("NM1", nitrogen_m1, 0.0, 1.0, " ", " NM1", 0, "N")
    carbon_m2 = calculateCoordinates(
        CI, CJ2, CK3, geo.CK3_CM2_length, geo.CJ2_CK3_CM2_angle, geo.CI_CJ2_CK3_CM2_diangle
    )
    CM2 = Atom("CM2", carbon_m2, 0.0, 1.0, " ", " CM2", 0, "C")
    carbon_m3 = calculateCoordinates(
        CI, CJ2, CK3, geo.CK3_CM3_length, geo.CJ2_CK3_CM3_angle, geo.CI_CJ2_CK3_CM3_diangle
    )
    CM3 = Atom("CM3", carbon_m3, 0.0, 1.0, " ", " CM3", 0, "C")
    carbon_n2 = calculateCoordinates(
        CJ2, CK3, CM2, geo.CM2_CN2_length, geo.CK3_CM2_CN2_angle, geo.CJ2_CK3_CM2_CN2_diangle
    )
    CN2 = Atom("CN2", carbon_n2, 0.0, 1.0, " ", " CN2", 0, "C")
    carbon_n3 = calculateCoordinates(
        CJ2, CK3, CM3, geo.CM3_CN3_length, geo.CK3_CM3_CN3_angle, geo.CJ2_CK3_CM3_CN3_diangle
    )
    CN3 = Atom("CN3", carbon_n3, 0.0, 1.0, " ", " CN3", 0, "C")
    carbon_q2 = calculateCoordinates(
        CK3, CM2, CN2, geo.CN2_CQ2_length, geo.CM2_CN2_CQ2_angle, geo.CK3_CM2_CN2_CQ2_diangle
    )
    CQ2 = Atom("CQ2", carbon_q2, 0.0, 1.0, " ", " CQ2", 0, "C")

    SG = [CB, CG2, CD2, CD3, NE1, CE2, CE3, CZ2, CZ3, CH2]
    CD1 = [CI, CJ2, CK2, CK3, NM1, CM2, CM3, CN2, CN3, CQ2]
    return SG, CD1


def make_res_of_type_aa(segID: int, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo: Geo,
                        geo1: Geo, geo2: Geo) -> Residue:
    res = Residue((" ", segID, " "), "LAA", "    ")
    res.add(N)
    res.add(CD1)
    res.add(CG)
    res.add(NB)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(SG)
    res.add(OD1)
    res.add(OD2)

    if isinstance(geo1, GlyGeo):
        SG_chain = Gly(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, AlaGeo):
        SG_chain = Ala(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, SerGeo):
        SG_chain = Ser(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, CysGeo):
        SG_chain = Cys(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, ValGeo):
        SG_chain = Val(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, IleGeo):
        SG_chain = Ile(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, LeuGeo):
        SG_chain = Leu(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, ThrGeo):
        SG_chain = Thr(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, ArgGeo):
        SG_chain = Arg(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, LysGeo):
        SG_chain = Lys(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, AspGeo):
        SG_chain = Asp(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, GluGeo):
        SG_chain = Glu(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, AsnGeo):
        SG_chain = Asn(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, GlnGeo):
        SG_chain = Gln(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, MetGeo):
        SG_chain = Met(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, HisGeo):
        SG_chain = His(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, ProGeo):
        SG_chain = Pro(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, PheGeo):
        SG_chain = Phe(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    elif isinstance(geo1, TyrGeo):
        SG_chain = Tyr(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]
    else:
        SG_chain = Trp(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo1)[0]

    if isinstance(geo2, GlyGeo):
        CD1_chain = Gly(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, AlaGeo):
        CD1_chain = Ala(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, SerGeo):
        CD1_chain = Ser(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, CysGeo):
        CD1_chain = Cys(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, ValGeo):
        CD1_chain = Val(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, IleGeo):
        CD1_chain = Ile(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, LeuGeo):
        CD1_chain = Leu(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, ThrGeo):
        CD1_chain = Thr(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, ArgGeo):
        CD1_chain = Arg(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, LysGeo):
        CD1_chain = Lys(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, AspGeo):
        CD1_chain = Asp(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, GluGeo):
        CD1_chain = Glu(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, AsnGeo):
        CD1_chain = Asn(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, GlnGeo):
        CD1_chain = Gln(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, MetGeo):
        CD1_chain = Met(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, HisGeo):
        CD1_chain = His(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, ProGeo):
        CD1_chain = Pro(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, PheGeo):
        CD1_chain = Phe(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    elif isinstance(geo2, TyrGeo):
        CD1_chain = Tyr(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]
    else:
        CD1_chain = Trp(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo2)[1]

    for i in range(len(CD1_chain)):
        res.add(CD1_chain[i])
    for i in range(len(SG_chain)):
        res.add(SG_chain[i])

    return res


# ,sidechain1: Union[Geo, str],sidechain2: Union[Geo, str]
def initialize_res(residue: Union[Geo, str],
                   sidechain1: Union[Geo, str], sidechain2: Union[Geo, str]) -> Structure:
    """Creates a new structure containing a single amino acid. The type and
    geometry of the amino acid are determined by the argument, which has to be
    either a geometry object or a single-letter amino acid code.
    The amino acid will be placed into chain A of model 0."""

    if isinstance(residue, Geo):
        geo = residue
    elif isinstance(residue, str):
        geo = geometry(residue)
    else:
        raise ValueError("Invalid residue argument:", residue)

    if isinstance(sidechain1, Geo):
        geo1 = sidechain1
    elif isinstance(sidechain1, str):
        geo1 = geometry(sidechain1)
    else:
        raise ValueError("Invalid residue argument:", sidechain1)

    if isinstance(sidechain2, Geo):
        geo2 = sidechain2
    elif isinstance(sidechain2, str):
        geo2 = geometry(sidechain2)
    else:
        raise ValueError("Invalid residue argument:", sidechain2)

    segID = 1
    AA = geo.residue_name
    CD1_N_length = geo.CD1_N_length
    CD1_CG_length = geo.CD1_CG_length
    CG_CD1_N_angle = geo.CG_CD1_N_angle

    CD1_coord = np.array([0.0, 0.0, 0.0])
    CG_coord = np.array([CD1_CG_length, 0, 0])
    N_coord = np.array(
        [
            CD1_N_length * math.cos(CG_CD1_N_angle * (math.pi / 180.0)),
            CD1_N_length * math.sin(CG_CD1_N_angle * (math.pi / 180.0)),
            0,
        ]
    )

    N = Atom("N", N_coord, 0.0, 1.0, " ", " N", 0, "N")
    CG = Atom("CG", CG_coord, 0.0, 1.0, " ", " CG", 0, "C")
    CD1 = Atom("CD1", CD1_coord, 0.0, 1.0, " ", " CD1", 0, "C")

    NB_CG_CD1_N_diangle = geo.NB_CG_CD1_N_diangle

    NB_CG_length = geo.NB_CG_length
    CA_NB_CG_angle = geo.CA_NB_CG_angle
    C_CA_NB_CG_diangle = geo.C_CA_NB_CG_diangle

    CA_NB_length = geo.CA_NB_length
    C_CA_length = geo.C_CA_length
    C_CA_NB_angle = geo.C_CA_NB_angle
    NB_CG_CD1_angle = geo.NB_CG_CD1_angle
    CA_NB_CG_CD1_diangle = geo.CA_NB_CG_CD1_diangle

    NB = calculateCoordinates(
        N, CD1, CG, NB_CG_length, NB_CG_CD1_angle, NB_CG_CD1_N_diangle
    )
    NB = Atom("NB", NB, 0.0, 1.0, " ", " NB", 0, "N")
    carbon_a = calculateCoordinates(
        CD1, CG, NB, CA_NB_length, CA_NB_CG_angle, CA_NB_CG_CD1_diangle
    )
    CA = Atom("CA", carbon_a, 0.0, 1.0, " ", " CA", 0, "C")
    carbon = calculateCoordinates(
        CG, NB, CA, C_CA_length, C_CA_NB_angle, C_CA_NB_CG_diangle
    )
    C = Atom("C", carbon, 0.0, 1.0, " ", " C", 0, "C")

    ##Create Carbonyl atom (to be moved later)
    C_O_length = geo.C_O_length
    CA_C_O_angle = geo.CA_C_O_angle
    NB_CA_C_O_diangle = geo.NB_CA_C_O_diangle

    carbonyl = calculateCoordinates(
        NB, CA, C, C_O_length, CA_C_O_angle, NB_CA_C_O_diangle
    )
    O = Atom("O", carbonyl, 0.0, 1.0, " ", " O", 0, "O")
    sulfur_g = calculateCoordinates(
        CD1, CG, NB, geo.NB_SG_length, geo.CG_NB_SG_angle, geo.CD1_CG_NB_SG_diangle
    )
    SG = Atom("SG", sulfur_g, 0.0, 1.0, " ", " SG", 0, "S")
    oxygen_d2 = calculateCoordinates(
        CA, NB, SG, geo.OD2_SG_length, geo.OD2_SG_NB_angle, geo.CA_NB_SG_OD2_diangle
    )
    OD2 = Atom("OD2", oxygen_d2, 0.0, 1.0, " ", " OD2", 0, "O")
    oxygen_d1 = calculateCoordinates(
        CG, NB, SG, geo.OD1_SG_length, geo.OD1_SG_NB_angle, geo.CG_NB_SG_OD1_diangle
    )
    OD1 = Atom("OD1", oxygen_d1, 0.0, 1.0, " ", " OD1", 0, "O")

    res = make_res_of_type_aa(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo, geo1, geo2)

    cha = Chain("A")
    cha.add(res)

    mod = Model(0)
    mod.add(cha)

    struc = Structure("X")
    struc.add(mod)
    return struc


def getReferenceResidue(structure: Structure) -> Residue:
    """Returns the last residue of chain A model 0 of the given structure.

    This function is a helper function that should not normally be called
    directly."""

    # If the following line doesn't work we're in trouble.
    # Likely initialize_res() wasn't called.
    resRef = structure[0]["A"].child_list[-1]

    # If the residue is not an amino acid we're in trouble.
    # Likely somebody is trying to append residues to an existing
    # structure that has non-amino-acid molecules in the chain.
    assert is_aa(resRef)

    return resRef


def add_residue_from_geo(structure: Structure, geo: Geo, geo1: Geo, geo2: Geo) -> Structure:
    """Adds a residue to chain A model 0 of the given structure, and
    returns the new structure. The residue to be added is determined by
    the geometry object given as second argument.

    This function is a helper function and should not normally be called
    directly. Call add_residue() instead."""
    resRef = getReferenceResidue(structure)
    AA = geo.residue_name
    segID = resRef.get_id()[1]
    segID += 1

    N_coord = calculateCoordinates(
        resRef["NB"], resRef["CA"], resRef["C"], geo.peptide_bond, geo.CA_C_N_angle, geo.N_C_CA_NB_diangle
    )
    N = Atom("N", N_coord, 0.0, 1.0, " ", " N", 0, "N")

    CD1_coord = calculateCoordinates(
        resRef["CA"], resRef["C"], N, geo.CD1_N_length, geo.CD1_N_C_angle, geo.CD1_N_C_CA_diangle
    )
    CD1 = Atom("CD1", CD1_coord, 0.0, 1.0, " ", " CD1", 0, "C")

    CG_coord = calculateCoordinates(resRef["C"], N, CD1, geo.CD1_CG_length, geo.CG_CD1_N_angle, geo.CG_CD1_N_C_diangle)
    CG = Atom("CG", CG_coord, 0.0, 1.0, " ", " CG", 0, "C")

    NB = calculateCoordinates(
        N, CD1, CG, geo.NB_CG_length, geo.NB_CG_CD1_angle, geo.NB_CG_CD1_N_diangle
    )
    NB = Atom("NB", NB, 0.0, 1.0, " ", " NB", 0, "N")
    carbon_a = calculateCoordinates(
        CD1, CG, NB, geo.CA_NB_length, geo.CA_NB_CG_angle, geo.CA_NB_CG_CD1_diangle
    )
    CA = Atom("CA", carbon_a, 0.0, 1.0, " ", " CA", 0, "C")
    carbon = calculateCoordinates(
        CG, NB, CA, geo.C_CA_length, geo.C_CA_NB_angle, geo.C_CA_NB_CG_diangle
    )
    C = Atom("C", carbon, 0.0, 1.0, " ", " C", 0, "C")

    carbonyl = calculateCoordinates(
        NB, CA, C, geo.C_O_length, geo.CA_C_O_angle, geo.NB_CA_C_O_diangle
    )
    O = Atom("O", carbonyl, 0.0, 1.0, " ", " O", 0, "O")
    sulfur_g = calculateCoordinates(
        CD1, CG, NB, geo.NB_SG_length, geo.CG_NB_SG_angle, geo.CD1_CG_NB_SG_diangle
    )
    SG = Atom("SG", sulfur_g, 0.0, 1.0, " ", " SG", 0, "S")
    oxygen_d2 = calculateCoordinates(
        CA, NB, SG, geo.OD2_SG_length, geo.OD2_SG_NB_angle, geo.CA_NB_SG_OD2_diangle
    )
    OD2 = Atom("OD2", oxygen_d2, 0.0, 1.0, " ", " OD2", 0, "O")
    oxygen_d1 = calculateCoordinates(
        CG, NB, SG, geo.OD1_SG_length, geo.OD1_SG_NB_angle, geo.CG_NB_SG_OD1_diangle
    )
    OD1 = Atom("OD1", oxygen_d1, 0.0, 1.0, " ", " OD1", 0, "O")

    res = make_res_of_type_aa(segID, N, CD1, CG, NB, CA, C, O, SG, OD1, OD2, geo, geo1, geo2)

    structure[0]["A"].add(res)
    return structure


def add_residue(
        structure: Structure, residue: Union[Geo, str], sidechain1: Union[Geo, str], sidechain2: Union[Geo, str],
) -> Structure:
    """Adds a residue to chain A model 0 of the given structure, and
    returns the new structure. The residue to be added can be specified
    in two ways: either as a geometry object (in which case
    the remaining arguments phi, psi_im1, and omega are ignored) or as a
    single-letter amino-acid code. In the latter case, the optional
    arguments phi, psi_im1, and omega specify the corresponding backbone
    angles.

    When omega is specified, it needs to be a value greater than or equal
    to -360. Values below -360 are ignored."""

    if isinstance(residue, Geo):
        geo = residue
    elif isinstance(residue, str):
        geo = geometry(residue)
    else:
        raise ValueError("Invalid residue argument:", residue)

    if isinstance(sidechain1, Geo):
        geo1 = sidechain1
    elif isinstance(sidechain1, str):
        geo1 = geometry(sidechain1)
    else:
        raise ValueError("Invalid residue argument:", sidechain1)

    if isinstance(sidechain2, Geo):
        geo2 = sidechain2
    elif isinstance(sidechain2, str):
        geo2 = geometry(sidechain2)
    else:
        raise ValueError("Invalid residue argument:", sidechain2)

    return add_residue_from_geo(structure, geo, geo1, geo2)


def add_terminal_OXT(structure: Structure, C_OXT_length: float = 1.23) -> Structure:
    """Adds a terminal oxygen atom ('OXT') to the last residue of chain A model 0 of the given structure, and returns the new structure. The OXT atom object will be contained in the last residue object of the structure.

This function should be used only when the structure object is completed and no further residues need to be appended."""

    rad = 180.0 / math.pi

    # obtain last residue infomation
    resRef = getReferenceResidue(structure)
    N_resRef = resRef["N"]
    CA_resRef = resRef["CA"]
    C_resRef = resRef["C"]
    O_resRef = resRef["O"]

    n_vec = N_resRef.get_vector()
    ca_vec = CA_resRef.get_vector()
    c_vec = C_resRef.get_vector()
    o_vec = O_resRef.get_vector()

    # geometry to bring together residue
    CA_C_OXT_angle = calc_angle(ca_vec, c_vec, o_vec) * rad
    N_CA_C_O_diangle = calc_dihedral(n_vec, ca_vec, c_vec, o_vec) * rad
    N_CA_C_OXT_diangle = N_CA_C_O_diangle - 180.0
    if N_CA_C_O_diangle < 0:
        N_CA_C_OXT_diangle = N_CA_C_O_diangle + 180.0

    # OXT atom creation
    OXT_coord = calculateCoordinates(
        N_resRef, CA_resRef, C_resRef, C_OXT_length, CA_C_OXT_angle, N_CA_C_OXT_diangle
    )
    OXT = Atom("OXT", OXT_coord, 0.0, 1.0, " ", "OXT", 0, "O")

    # modify last residue of the structure to contain the OXT atom
    resRef.add(OXT)
    return structure
