"""This module is part of the PeptideBuilder library,
written by Matthew Z. Tien, Dariya K. Sydykova,
Austin G. Meyer, and Claus O. Wilke.

The Geometry module contains the default geometries of
all 20 amino acids. The main function to be used is the
geometry() function, which returns the default geometry
for the requested amino acid.

This file is provided to you under the MIT License."""

import random
from typing import List


class Geo:
    """Geometry base class"""

    residue_name: str

    # Geometry to bring together residue
    peptide_bond: float
    CA_C_N_angle: float
    C_N_CA_angle: float

    # Backbone coordinates
    N_CA_C_angle: float
    CA_N_length: float
    CA_C_length: float
    phi: float
    psi_im1: float
    omega: float

    # Carbonyl atom
    C_O_length: float
    CA_C_O_angle: float
    N_CA_C_O_diangle: float

    def __repr__(self) -> str:
        repr = ""
        for var in self.__dict__:
            repr += "%s = %s\n" % (var, self.__dict__[var])
        return repr


class GlyGeo(Geo):
    """Geometry of Glycine"""

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.residue_name = "G"


class AlaGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.residue_name = "A"


class SerGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_OJ_length = 1.417
        self.CD1_CI_OJ_angle = 110.773
        self.N_CD1_CI_OJ_diangle = -63.3

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_OG_length = 1.417
        self.SG_CB_OG_angle = 110.773
        self.NB_SG_CB_OG_diangle = -63.3

        self.residue_name = "S"


class CysGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_SG1_length = 1.808
        self.CD1_CI_SG1_angle = 113.8169
        self.N_CD1_CI_SG1_diangle = -62.2

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_SG2_length = 1.808
        self.SG_CB_SG2_angle = 113.8169
        self.NB_SG_CB_SG2_diangle = -62.2

        self.residue_name = "C"


class ValGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ1_length = 1.527
        self.CD1_CI_CJ1_angle = 110.7
        self.N_CD1_CI_CJ1_diangle = 177.2

        self.CI_CJ2_length = 1.527
        self.CD1_CI_CJ2_angle = 110.4
        self.N_CD1_CI_CJ2_diangle = -63.3

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG1_length = 1.527
        self.SG_CB_CG1_angle = 110.7
        self.NB_SG_CB_CG1_diangle = 177.2

        self.CB_CG2_length = 1.527
        self.SG_CB_CG2_angle = 110.4
        self.NB_SG_CB_CG2_diangle = -63.3

        self.residue_name = "V"


class IleGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ1_length = 1.527
        self.CD1_CI_CJ1_angle = 110.7
        self.N_CD1_CI_CJ1_diangle = 59.7

        self.CI_CJ2_length = 1.527
        self.CD1_CI_CJ2_angle = 110.4
        self.N_CD1_CI_CJ2_diangle = -61.6

        self.CJ1_CK2_length = 1.52
        self.CI_CJ1_CK2_angle = 113.97
        self.CD1_CI_CJ1_CK2_diangle = 169.8

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG1_length = 1.527
        self.SG_CB_CG1_angle = 110.7
        self.NB_SG_CB_CG1_diangle = 59.7

        self.CB_CG2_length = 1.527
        self.SG_CB_CG2_angle = 110.4
        self.NB_SG_CB_CG2_diangle = -61.6

        self.CG1_CD2_length = 1.52
        self.CB_CG1_CD2_angle = 113.97
        self.SG_CB_CG1_CD2_diangle = 169.8

        self.residue_name = "I"



class LeuGeo(Geo):
   

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915
        
        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ2_length = 1.53
        self.CD1_CI_CJ2_angle = 116.10
        self.N_CD1_CI_CJ2_diangle = -60.1

        self.CJ2_CK2_length = 1.524
        self.CI_CJ2_CK2_angle = 110.27
        self.CD1_CI_CJ2_CK2_diangle = 174.9

        self.CJ2_CK3_length = 1.525
        self.CI_CJ2_CK3_angle = 110.58
        self.CD1_CI_CJ2_CK3_diangle = 66.7
        
        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG2_length = 1.53
        self.SG_CB_CG2_angle = 116.10
        self.NB_SG_CB_CG2_diangle = -60.1

        self.CG2_CD2_length = 1.524
        self.CB_CG2_CD2_angle = 110.27
        self.SG_CB_CG2_CD2_diangle = 174.9

        self.CG2_CD3_length = 1.525
        self.CB_CG2_CD3_angle = 110.58
        self.SG_CB_CG2_CD3_diangle = 66.7

        self.residue_name = "L"


class ThrGeo(Geo):
   

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_OJ1_length = 1.43
        self.CD1_CI_OJ1_angle = 109.18
        self.N_CD1_CI_OJ1_diangle = 60.0

        self.CI_CJ2_length = 1.53
        self.CD1_CI_CJ2_angle = 111.13
        self.N_CD1_CI_CJ2_diangle = -60.3

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_OG1_length = 1.43
        self.SG_CB_OG1_angle = 109.18
        self.NB_SG_CB_OG1_diangle = 60.0

        self.CB_CG2_length = 1.53
        self.SG_CB_CG2_angle = 111.13
        self.NB_SG_CB_CG2_diangle = -60.3

        self.residue_name = "T"


class ArgGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ2_length = 1.52
        self.CD1_CI_CJ2_angle = 113.83
        self.N_CD1_CI_CJ2_diangle = -65.2

        self.CJ2_CK2_length = 1.52
        self.CI_CJ2_CK2_angle = 111.79
        self.CD1_CI_CJ2_CK2_diangle = -179.2

        self.CK2_NM_length = 1.46
        self.CJ2_CK2_NM_angle = 111.68
        self.CI_CJ2_CK2_NM_diangle = -179.3

        self.NM_CN_length = 1.33
        self.CK2_NM_CN_angle = 124.79
        self.CJ2_CK2_NM_CN_diangle = -178.7

        self.CN_NQ1_length = 1.33
        self.NM_CN_NQ1_angle = 120.64
        self.CK2_NM_CN_NQ1_diangle = 0.0

        self.CN_NQ2_length = 1.33
        self.NM_CN_NQ2_angle = 119.63
        self.CK2_NM_CN_NQ2_diangle = 180.0

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG2_length = 1.52
        self.SG_CB_CG2_angle = 113.83
        self.NB_SG_CB_CG2_diangle = -65.2

        self.CG2_CD2_length = 1.52
        self.CB_CG2_CD2_angle = 111.79
        self.SG_CB_CG2_CD2_diangle = -179.2

        self.CD2_NE_length = 1.46
        self.CG2_CD2_NE_angle = 111.68
        self.CB_CG2_CD2_NE_diangle = -179.3

        self.NE_CZ_length = 1.33
        self.CD2_NE_CZ_angle = 124.79
        self.CG2_CD2_NE_CZ_diangle = -178.7

        self.CZ_NH1_length = 1.33
        self.NE_CZ_NH1_angle = 120.64
        self.CD2_NE_CZ_NH1_diangle = 0.0

        self.CZ_NH2_length = 1.33
        self.NE_CZ_NH2_angle = 119.63
        self.CD2_NE_CZ_NH2_diangle = 180.0

        self.residue_name = "R"

    
class LysGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ2_length = 1.52
        self.CD1_CI_CJ2_angle = 113.83
        self.N_CD1_CI_CJ2_diangle = -64.5

        self.CJ2_CK_length = 1.52
        self.CI_CJ2_CK_angle = 111.79
        self.CD1_CI_CJ2_CK_diangle = -178.1

        self.CK_CM_length = 1.46
        self.CJ2_CK_CM_angle = 111.68
        self.CI_CJ2_CK_CM_diangle = -179.6

        self.CM_NZ2_length = 1.33
        self.CK_CM_NZ2_angle = 124.79
        self.CJ2_CK_CM_NZ2_diangle = 179.6
        
        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG2_length = 1.52
        self.SG_CB_CG2_angle = 113.83
        self.NB_SG_CB_CG2_diangle = -64.5

        self.CG2_CD_length = 1.52
        self.CB_CG2_CD_angle = 111.79
        self.SG_CB_CG2_CD_diangle = -178.1

        self.CD_CE_length = 1.46
        self.CG2_CD_CE_angle = 111.68
        self.CB_CG2_CD_CE_diangle = -179.6

        self.CE_NZ_length = 1.33
        self.CD_CE_NZ_angle = 124.79
        self.CG2_CD_CE_NZ_diangle = 179.6

        self.residue_name = "K"


class AspGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ2_length = 1.52
        self.CD1_CI_CJ2_angle = 113.06
        self.N_CD1_CI_CJ2_diangle = -66.4

        self.CJ2_OK3_length = 1.25
        self.CI_CJ2_OK3_angle = 119.22
        self.CD1_CI_CJ2_OK3_diangle = -46.7

        self.CJ2_OK4_length = 1.25
        self.CI_CJ2_OK4_angle = 118.218
        self.CD1_CI_CJ2_OK4_diangle = 180 + self.CD1_CI_CJ2_OK3_diangle

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG2_length = 1.52
        self.SG_CB_CG2_angle = 113.06
        self.NB_SG_CB_CG2_diangle = -66.4

        self.CG2_OD3_length = 1.25
        self.CB_CG2_OD3_angle = 119.22
        self.SG_CB_CG2_OD3_diangle = -46.7

        self.CG2_OD4_length = 1.25
        self.CB_CG2_OD4_angle = 118.218
        self.SG_CB_CG2_OD4_diangle = 180 + self.SG_CB_CG2_OD3_diangle

        self.residue_name = "D"


class AsnGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ2_length = 1.52
        self.CD1_CI_CJ2_angle = 112.62
        self.N_CD1_CI_CJ2_diangle = -65.5

        self.CJ2_OK3_length = 1.23
        self.CI_CJ2_OK3_angle = 120.85
        self.CD1_CI_CJ2_OK3_diangle = -58.3

        self.CJ2_NK2_length = 1.33
        self.CI_CJ2_NK2_angle = 116.48
        self.CD1_CI_CJ2_NK2_diangle = 180.0 + self.CD1_CI_CJ2_OK3_diangle

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG2_length = 1.52
        self.SG_CB_CG2_angle = 112.62
        self.NB_SG_CB_CG2_diangle = -65.5

        self.CG2_OD3_length = 1.23
        self.CB_CG2_OD3_angle = 120.85
        self.SG_CB_CG2_OD3_diangle = -58.3

        self.CG2_ND2_length = 1.33
        self.CB_CG2_ND2_angle = 116.48
        self.SG_CB_CG2_ND2_diangle = 180.0 + self.SG_CB_CG2_OD3_diangle

        self.residue_name = "N"


class GluGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ2_length = 1.52
        self.CD1_CI_CJ2_angle = 113.82
        self.N_CD1_CI_CJ2_diangle = -63.8

        self.CJ2_CK2_length = 1.52
        self.CI_CJ2_CK2_angle = 113.31
        self.CD1_CI_CJ2_CK2_diangle = -179.8

        self.CK2_OM1_length = 1.25
        self.CJ2_CK2_OM1_angle = 119.02
        self.CI_CJ2_CK2_OM1_diangle = -6.2

        self.CK2_OM2_length = 1.25
        self.CJ2_CK2_OM2_angle = 118.08
        self.CI_CJ2_CK2_OM2_diangle = 180.0 + self.CI_CJ2_CK2_OM1_diangle

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG2_length = 1.52
        self.SG_CB_CG2_angle = 113.82
        self.NB_SG_CB_CG2_diangle = -63.8

        self.CG2_CD2_length = 1.52
        self.CB_CG2_CD2_angle = 113.31
        self.SG_CB_CG2_CD2_diangle = -179.8

        self.CD2_OE1_length = 1.25
        self.CG2_CD2_OE1_angle = 119.02
        self.CB_CG2_CD2_OE1_diangle = -6.2

        self.CD2_OE2_length = 1.25
        self.CG2_CD2_OE2_angle = 118.08
        self.CB_CG2_CD2_OE2_diangle = 180.0 + self.CB_CG2_CD2_OE1_diangle

        self.residue_name = "E"


class GlnGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ2_length = 1.52
        self.CD1_CI_CJ2_angle = 113.75
        self.N_CD1_CI_CJ2_diangle = -60.2

        self.CJ2_CK_length = 1.52
        self.CI_CJ2_CK_angle = 112.78
        self.CD1_CI_CJ2_CK_diangle = -69.6

        self.CK_OM1_length = 1.24
        self.CJ2_CK_OM1_angle = 120.86
        self.CI_CJ2_CK_OM1_diangle = -50.5

        self.CK_NM2_length = 1.33
        self.CJ2_CK_NM2_angle = 116.50
        self.CI_CJ2_CK_NM2_diangle = 180 + self.CI_CJ2_CK_OM1_diangle

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG2_length = 1.52
        self.SG_CB_CG2_angle = 113.75
        self.NB_SG_CB_CG2_diangle = -60.2

        self.CG2_CD_length = 1.52
        self.CB_CG2_CD_angle = 112.78
        self.SG_CB_CG2_CD_diangle = -69.6

        self.CD_OE1_length = 1.24
        self.CG2_CD_OE1_angle = 120.86
        self.CB_CG2_CD_OE1_diangle = -50.5

        self.CD_NE2_length = 1.33
        self.CG2_CD_NE2_angle = 116.50
        self.CB_CG2_CD_NE2_diangle = 180 + self.CB_CG2_CD_OE1_diangle

        self.residue_name = "Q"


class MetGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ2_length = 1.52
        self.CD1_CI_CJ2_angle = 113.68
        self.N_CD1_CI_CJ2_diangle = -64.4

        self.CJ2_SK_length = 1.81
        self.CI_CJ2_SK_angle = 112.69
        self.CD1_CI_CJ2_SK_diangle = -179.6

        self.SK_CM_length = 1.79
        self.CJ2_SK_CM_angle = 100.61
        self.CI_CJ2_SK_CM_diangle = 70.1

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG2_length = 1.52
        self.SG_CB_CG2_angle = 113.68
        self.NB_SG_CB_CG2_diangle = -64.4

        self.CG2_SD_length = 1.81
        self.CB_CG2_SD_angle = 112.69
        self.SG_CB_CG2_SD_diangle = -179.6

        self.SD_CE_length = 1.79
        self.CG2_SD_CE_angle = 100.61
        self.CB_CG2_SD_CE_diangle = 70.1

        self.residue_name = "M"


class HisGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ2_length = 1.49
        self.CD1_CI_CJ2_angle = 113.74
        self.N_CD1_CI_CJ2_diangle = -63.2

        self.CJ2_NK1_length = 1.38
        self.CI_CJ2_NK1_angle = 122.85
        self.CD1_CI_CJ2_NK1_diangle = -75.7

        self.CJ2_CK2_length = 1.35
        self.CI_CJ2_CK2_angle = 130.61
        self.CD1_CI_CJ2_CK2_diangle = 180.0 + self.CD1_CI_CJ2_NK1_diangle

        self.NK1_CM2_length = 1.32
        self.CJ2_NK1_CM2_angle = 108.5
        self.CI_CJ2_NK1_CM2_diangle = 180.0

        self.CK2_NM2_length = 1.35
        self.CJ2_CK2_NM2_angle = 108.5
        self.CI_CJ2_CK2_NM2_diangle = 180.0

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG2_length = 1.49
        self.SG_CB_CG2_angle = 113.74
        self.NB_SG_CB_CG2_diangle = -63.2

        self.CG2_ND1_length = 1.38
        self.CB_CG2_ND1_angle = 122.85
        self.SG_CB_CG2_ND1_diangle = -75.7

        self.CG2_CD2_length = 1.35
        self.CB_CG2_CD2_angle = 130.61
        self.SG_CB_CG2_CD2_diangle = 180.0 + self.SG_CB_CG2_ND1_diangle

        self.ND1_CE2_length = 1.32
        self.CG2_ND1_CE2_angle = 108.5
        self.CB_CG2_ND1_CE2_diangle = 180.0

        self.CD2_NE2_length = 1.35
        self.CG2_CD2_NE2_angle = 108.5
        self.CB_CG2_CD2_NE2_diangle = 180.0

        self.residue_name = "H"


class ProGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ2_length = 1.49
        self.CD1_CI_CJ2_angle = 104.21
        self.N_CD1_CI_CJ2_diangle = 29.6

        self.CJ2_CK_length = 1.50
        self.CI_CJ2_CK_angle = 105.03
        self.CD1_CI_CJ2_CK_diangle = -34.8

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG2_length = 1.49
        self.SG_CB_CG2_angle = 104.21
        self.NB_SG_CB_CG2_diangle = 29.6

        self.CG2_CD_length = 1.50
        self.CB_CG2_CD_angle = 105.03
        self.SG_CB_CG2_CD_diangle = -34.8

        self.residue_name = "P"


class PheGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ2_length = 1.50
        self.CD1_CI_CJ2_angle = 113.85
        self.N_CD1_CI_CJ2_diangle = -64.7

        self.CJ2_CK2_length = 1.39
        self.CI_CJ2_CK2_angle = 120.0
        self.CD1_CI_CJ2_CK2_diangle = 93.3

        self.CJ2_CK3_length = 1.39
        self.CI_CJ2_CK3_angle = 120.0
        self.CD1_CI_CJ2_CK3_diangle = self.CD1_CI_CJ2_CK2_diangle - 180.0

        self.CK2_CM2_length = 1.39
        self.CJ2_CK2_CM2_angle = 120.0
        self.CI_CJ2_CK2_CM2_diangle = 180.0

        self.CK3_CM3_length = 1.39
        self.CJ2_CK3_CM3_angle = 120.0
        self.CI_CJ2_CK3_CM3_diangle = 180.0

        self.CM2_CN_length = 1.39
        self.CK2_CM2_CN_angle = 120.0
        self.CJ2_CK2_CM2_CN_diangle = 0.0

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG2_length = 1.50
        self.SG_CB_CG2_angle = 113.85
        self.NB_SG_CB_CG2_diangle = -64.7

        self.CG2_CD2_length = 1.39
        self.CB_CG2_CD2_angle = 120.0
        self.SG_CB_CG2_CD2_diangle = 93.3

        self.CG2_CD3_length = 1.39
        self.CB_CG2_CD3_angle = 120.0
        self.SG_CB_CG2_CD3_diangle = self.SG_CB_CG2_CD2_diangle - 180.0

        self.CD2_CE2_length = 1.39
        self.CG2_CD2_CE2_angle = 120.0
        self.CB_CG2_CD2_CE2_diangle = 180.0

        self.CD3_CE3_length = 1.39
        self.CG2_CD3_CE3_angle = 120.0
        self.CB_CG2_CD3_CE3_diangle = 180.0

        self.CE2_CZ_length = 1.39
        self.CD2_CE2_CZ_angle = 120.0
        self.CG2_CD2_CE2_CZ_diangle = 0.0

        self.residue_name = "F"



class TyrGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ2_length = 1.51
        self.CD1_CI_CJ2_angle = 113.8
        self.N_CD1_CI_CJ2_diangle = -64.3

        self.CJ2_CK2_length = 1.39
        self.CI_CJ2_CK2_angle = 120.98
        self.CD1_CI_CJ2_CK2_diangle = 93.1

        self.CJ2_CK3_length = 1.39
        self.CI_CJ2_CK3_angle = 120.82
        self.CD1_CI_CJ2_CK3_diangle = self.CD1_CI_CJ2_CK2_diangle + 180.0

        self.CK2_CM2_length = 1.39
        self.CJ2_CK2_CM2_angle = 120.0
        self.CI_CJ2_CK2_CM2_diangle = 180.0

        self.CK3_CM3_length = 1.39
        self.CJ2_CK3_CM3_angle = 120.0
        self.CI_CJ2_CK3_CM3_diangle = 180.0

        self.CM2_CN_length = 1.39
        self.CK2_CM2_CN_angle = 120.0
        self.CJ2_CK2_CM2_CN_diangle = 0.0

        self.CN_OQ_length = 1.39
        self.CM2_CN_OQ_angle = 119.78
        self.CK2_CM2_CN_OQ_diangle = 180.0

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG2_length = 1.51
        self.SG_CB_CG2_angle = 113.8
        self.NB_SG_CB_CG2_diangle = -64.3

        self.CG2_CD2_length = 1.39
        self.CB_CG2_CD2_angle = 120.98
        self.SG_CB_CG2_CD2_diangle = 93.1

        self.CG2_CD3_length = 1.39
        self.CB_CG2_CD3_angle = 120.82
        self.SG_CB_CG2_CD3_diangle = self.SG_CB_CG2_CD2_diangle + 180.0

        self.CD2_CE2_length = 1.39
        self.CG2_CD2_CE2_angle = 120.0
        self.CB_CG2_CD2_CE2_diangle = 180.0

        self.CD3_CE3_length = 1.39
        self.CG2_CD3_CE3_angle = 120.0
        self.CB_CG2_CD3_CE3_diangle = 180.0

        self.CE2_CZ_length = 1.39
        self.CD2_CE2_CZ_angle = 120.0
        self.CG2_CD2_CE2_CZ_diangle = 0.0

        self.CZ_OH_length = 1.39
        self.CE2_CZ_OH_angle = 119.78
        self.CD2_CE2_CZ_OH_diangle = 180.0

        self.residue_name = "Y"


class TrpGeo(Geo):

    def __init__(self):
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.51
        self.N_C_CA_NB_diangle = -146.25

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 120.45
        self.CD1_N_C_CA_diangle = -173.54

        self.CD1_CG_length = 1.54
        self.CG_CD1_N_angle = 112.261
        self.CG_CD1_N_C_diangle = -137.91

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 114.47
        self.NB_CG_CD1_N_diangle = 66.24

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 117.63
        self.CA_NB_CG_CD1_diangle = -121.13

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = 86.201

        self.C_O_length = 1.51
        self.CA_C_O_angle = 119.85
        self.NB_CA_C_O_diangle = 44.99

        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle = 115.052
        self.CD1_CG_NB_SG_diangle = 94.641

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        self.CG_NB_SG_OD1_diangle = -41.980

        self.OD2_SG_length = 1.42
        self.OD2_SG_NB_angle = 106.396
        self.CA_NB_SG_OD2_diangle = 38.915

        self.CI_CD1_length = 1.530
        self.CI_CD1_CG_angle = 108.23
        self.CI_CD1_CG_NB_diangle = -175.64

        self.CI_CJ2_length = 1.50
        self.CD1_CI_CJ2_angle = 114.10
        self.N_CD1_CI_CJ2_diangle = -66.4

        self.CJ2_CK2_length = 1.37
        self.CI_CJ2_CK2_angle = 127.07
        self.CD1_CI_CJ2_CK2_diangle = 96.3

        self.CJ2_CK3_length = 1.43
        self.CI_CJ2_CK3_angle = 126.66
        self.CD1_CI_CJ2_CK3_diangle = self.CD1_CI_CJ2_CK2_diangle - 180.0

        self.CK2_NM1_length = 1.38
        self.CJ2_CK2_NM1_angle = 108.5
        self.CI_CJ2_CK2_NM1_diangle = 180.0

        self.CK3_CM2_length = 1.40
        self.CJ2_CK3_CM2_angle = 108.5
        self.CI_CJ2_CK3_CM2_diangle = 180.0

        self.CK3_CM3_length = 1.40
        self.CJ2_CK3_CM3_angle = 133.83
        self.CI_CJ2_CK3_CM3_diangle = 0.0

        self.CM2_CN2_length = 1.40
        self.CK3_CM2_CN2_angle = 120.0
        self.CJ2_CK3_CM2_CN2_diangle = 180.0

        self.CM3_CN3_length = 1.40
        self.CK3_CM3_CN3_angle = 120.0
        self.CJ2_CK3_CM3_CN3_diangle = 180.0

        self.CN2_CQ2_length = 1.40
        self.CM2_CN2_CQ2_angle = 120.0
        self.CK3_CM2_CN2_CQ2_diangle = 0.0

        self.SG_CB_length = 1.52
        self.NB_SG_CB_angle = 109.5
        self.CG_NB_SG_CB_diangle = 72.68

        self.CB_CG2_length = 1.50
        self.SG_CB_CG2_angle = 114.10
        self.NB_SG_CB_CG2_diangle = -66.4

        self.CG2_CD2_length = 1.37
        self.CB_CG2_CD2_angle = 127.07
        self.SG_CB_CG2_CD2_diangle = 96.3

        self.CG2_CD3_length = 1.43
        self.CB_CG2_CD3_angle = 126.66
        self.SG_CB_CG2_CD3_diangle = self.SG_CB_CG2_CD2_diangle - 180.0

        self.CD2_NE1_length = 1.38
        self.CG2_CD2_NE1_angle = 108.5
        self.CB_CG2_CD2_NE1_diangle = 180.0

        self.CD3_CE2_length = 1.40
        self.CG2_CD3_CE2_angle = 108.5
        self.CB_CG2_CD3_CE2_diangle = 180.0

        self.CD3_CE3_length = 1.40
        self.CG2_CD3_CE3_angle = 133.83
        self.CB_CG2_CD3_CE3_diangle = 0.0

        self.CE2_CZ2_length = 1.40
        self.CD3_CE2_CZ2_angle = 120.0
        self.CG2_CD3_CE2_CZ2_diangle = 180.0

        self.CE3_CZ3_length = 1.40
        self.CD3_CE3_CZ3_angle = 120.0
        self.CG2_CD3_CE3_CZ3_diangle = 180.0

        self.CZ2_CH2_length = 1.40
        self.CE2_CZ2_CH2_angle = 120.0
        self.CD3_CE2_CZ2_CH2_diangle = 0.0

        self.residue_name = "W"




def geometry(AA: str) -> Geo:
    """Generates the geometry of the requested amino acid.
    The amino acid needs to be specified by its single-letter
    code. If an invalid code is specified, the function
    returns the geometry of Glycine."""
    if AA == "G":
        return GlyGeo()
    elif AA == "A":
        return AlaGeo()
    elif AA == "S":
        return SerGeo()
    elif AA == "C":
        return CysGeo()
    elif AA == "V":
        return ValGeo()
    elif AA == "I":
        return IleGeo()
    elif AA == "L":
        return LeuGeo()
    elif AA == "T":
        return ThrGeo()
    elif AA == "R":
        return ArgGeo()
    elif AA == "K":
        return LysGeo()
    elif AA == "D":
        return AspGeo()
    elif AA == "E":
        return GluGeo()
    elif AA == "N":
        return AsnGeo()
    elif AA == "Q":
        return GlnGeo()
    elif AA == "M":
        return MetGeo()
    elif AA == "H":
        return HisGeo()
    elif AA == "P":
        return ProGeo()
    elif AA == "F":
        return PheGeo()
    elif AA == "Y":
        return TyrGeo()
    elif AA == "W":
        return TrpGeo()
