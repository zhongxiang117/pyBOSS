from PyBOSS.init_solutes import InitSolutes, calc_dipoles
from PyBOSS.utils import read_ff_file, read_zmat, read_oplsaa_file

from PyBOSS.Tests.utils import allclose
from PyBOSS.Tests import *

import unittest


class Test_InitSolutes(unittest.TestCase):
    solutezmat = read_zmat(FILE_ZMAT)
    bond_pars, dihedral_pars = read_ff_file(FILE_PAR,msg=False)
    oplsaa_bond_pars, oplsaa_angle_pars = read_oplsaa_file(FILE_OPLS_PAR)
    kws = {
        'rc0':0.5,  'rc1':0.0,  'rc2':1.0,  'qmscale':1.2,
        'bond_pars':bond_pars,  'dihedral_pars':dihedral_pars,
        'oplsaa_bond_pars':oplsaa_bond_pars,  'oplsaa_angle_pars':oplsaa_angle_pars,
        'solutezmat':solutezmat,
        'ncent1':3,  'ncent2':3,
    }
    tol = 0.00001
    def test_get_xyzs(self):
        xyz_initial = [
            [0.000000 ,  0.000000 ,  0.000000    ],
            [1.000000 ,  0.000000 ,  0.000000    ],
            [0.000000 ,  1.000000 ,  0.000000    ],
            [-1.272211,   1.000000,   0.000000   ], 
            [-1.951191,   3.092603,   0.000000   ], 
            [-1.314016,   3.579578,   -1.111939  ], 
            [0.115010 ,  3.728199 ,  -1.017984   ], 
            [0.794482 ,  1.636670 ,  -1.079801   ], 
            [2.229705 ,  1.561497 ,  -0.960377   ], 
            [0.721120 ,  0.465778 ,  1.152697    ],
            [2.800125 ,  0.933158 ,  0.098094    ],
            [2.076044 ,  0.434700 ,  1.158720    ],
            [4.274558 ,  0.776325 ,  0.153827    ],
            [-1.541183,   3.224955,   1.010660   ], 
            [-3.022894,   2.874651,   -0.027937  ], 
            [-1.794869,   3.618911,   -2.098826  ], 
            [0.679429 ,  4.061883 ,  -1.891697   ], 
            [0.569278 ,  3.960164 ,  -0.048508   ], 
            [0.372282 ,  1.601078 ,  -2.104217   ], 
            [2.831750 ,  1.993660 ,  -1.775593   ], 
            [2.632522 ,  0.015774 ,  2.012175    ],
            [0.143385 ,  0.065823 ,  2.002397    ],
            [4.692515 ,  1.268777 ,  1.041707    ],
            [4.783623 ,  1.197427 ,  -0.723578   ], 
            [4.555057 ,  -0.283254,   0.216539   ], 
        ]
        xyz_final = [
            [0.000000 ,  0.000000 ,   0.000000      ],
            [1.000000 ,  0.000000 ,   0.000000      ],
            [0.000000 ,  1.000000 ,   0.000000      ],
            [-1.272211,   1.000000,    0.000000     ], 
            [-1.963537,   3.130650,    0.000000     ], 
            [-1.326361,   3.617625,    -1.111939    ], 
            [0.102665 ,  3.766247 ,   -1.017984     ], 
            [0.782137 ,  1.674717 ,   -1.079801     ], 
            [2.217360 ,  1.599545 ,   -0.960377     ], 
            [0.738509 ,  0.444898 ,   1.131629      ],
            [2.787737 ,  0.946889 ,   0.083300      ],
            [2.093559 ,  0.421894 ,   1.120248      ],
            [4.263760 ,  0.798441 ,   0.114560      ],
            [-1.553528,   3.263002,    1.010660     ], 
            [-3.035239,   2.912699,    -0.027937    ], 
            [-1.807214,   3.656958,    -2.098826    ], 
            [0.667084 ,  4.099930 ,   -1.891697     ], 
            [0.556933 ,  3.998211 ,   -0.048508     ], 
            [0.349160 ,  1.659925 ,   -2.100222     ], 
            [2.819540 ,  2.052534 ,   -1.764106     ], 
            [2.663045 ,  -0.012349,    1.957310     ], 
            [0.173732 ,  0.022840 ,   1.979347      ],
            [4.690656 ,  1.267461 ,   1.010835      ],
            [4.758124 ,  1.248510 ,   -0.756851     ], 
            [4.551810 ,  -0.260604,    0.142310     ], 
        ]
        cis = InitSolutes(**self.kws)
        cis.get_xyzs()
        self.assertTrue(allclose(cis.solutesdata['xyz_initial'], xyz_initial, self.tol))
        self.assertTrue(allclose(cis.solutesdata['xyz_final'], xyz_final, self.tol))

    def test_get_pars(self):
        cis = InitSolutes(**self.kws)
        cis.get_pars()
        atypes = [0,0,6,8,6,6,6,6,6,6,6,6,6,1,1,1,1,1,1,1,1,1,1,1,1]
        self.assertTrue(allclose(cis.solutesdata['reference']['atomtypes'], atypes, self.tol))
        self.assertTrue(allclose(cis.solutesdata['first']['atomtypes'], atypes, self.tol))
        self.assertTrue(allclose(cis.solutesdata['second']['atomtypes'], atypes, self.tol))

    def test_calc_neighbors(self):
        cis = InitSolutes(**self.kws)
        cis.get_xyzs()
        cis.get_pars()
        cis.calc_neighbors()
        nbors = {       # care! index starts from 1
           1:  [],
           2:  [],
           3:  [ 4,  8,  10],
           4:  [ 3         ],
           5:  [ 6, 14,  15],
           6:  [ 5,  7,  16],
           7:  [ 6, 17,  18],
           8:  [ 3,  9,  19],
           9:  [ 8, 11,  20],
          10:  [ 3, 12,  22],
          11:  [ 9, 12,  13],
          12:  [10, 11,  21],
          13:  [11, 23,  24,  25],
          14:  [ 5  ],
          15:  [ 5  ],
          16:  [ 6  ],
          17:  [ 7  ],
          18:  [ 7  ],
          19:  [ 8  ],
          20:  [ 9  ],
          21:  [12  ],
          22:  [10  ],
          23:  [13  ],
          24:  [13  ],
          25:  [13  ],
        }
        pnbors = [[i-1 for i in sorted(nbors[k])] for k in range(1,26)]
        calcnb = [sorted(list(cis.nbor[k])) for k in range(0,25)]
        self.assertTrue(allclose(pnbors, calcnb, self.tol))

    def test_symmetrize_charges(self):
        charges = [
            0.0,  0.0,      # two dummy atoms
            0.123782,  -0.324020,  -0.122514,  -0.131094,  -0.122398,  -0.090627,  -0.110795,
            -0.107795, -0.041283,  -0.098389,  -0.241108,  0.113690,  0.114137,  0.111586,
            0.105579,  0.109568,   0.110946,  0.103171,  0.110094,  0.112362,  0.090699,
            0.091562,  0.092955
        ]
        expects = [
            0.0,  0.0,
            0.148538,  -0.388824,  -0.147017,  -0.157313,  -0.146878,  -0.108752,  -0.132954,
            -0.129354,  -0.049540,  -0.118067,  -0.289330,  0.136696,  0.136696,   0.133903,
            0.129088,  0.129088,   0.133135,   0.123805,   0.132113,   0.134834,   0.110086,
            0.110086,  0.110086,
        ]
        cis = InitSolutes(**self.kws)
        cis.get_xyzs()
        cis.get_pars()
        cis.calc_neighbors()
        cis.calc_atoms_symmetry_list()
        cis.calc_coulombic_pars(charges)
        proexpects = [i*cis.dsqesq for i in expects]
        self.assertTrue(allclose(cis.solutesdata['reference']['Q'], proexpects, self.tol))

    def test_pert_bond_pars(self):
        bk = [0.0,  0.0, 0.0,       # for non-defined atoms
            570,  350,  350,  350,  350,  350,  350,  350 , 437.27,  437.27,  437.27,
            437.27,   437.27,  437.27,  437.27,  437.27,  437.27,  437.27,  437.27,  437.27
        ]
        br = [0.0, 0.0, 0.0,
            1.229,  1.51,  1.51,  1.51,  1.51,  1.51,  1.51,  1.51,  1.01,  1.01,  1.01,
            1.01,   1.01,  1.01,  1.01,  1.01,  1.01,  1.01,  1.01,  1.01
        ]
        bd = [0.05, 0.05, 0.05,
            0.0244,  0.0312,  0.0312,  0.0312,  0.0312,  0.0312,  0.0312,  0.0312,  0.0279,
            0.0279,  0.0279,  0.0279,  0.0279,  0.0279,  0.0279,  0.0279,  0.0279,  0.0279,
            0.0279,  0.0279
        ]
        cis = InitSolutes(**self.kws)
        cis.get_xyzs()
        cis.get_pars()
        cis.calc_neighbors()
        cis.calc_pert_bond_pars()
        self.assertTrue(allclose(cis.solutesdata['reference']['pert_bonds_r'], br, 0.01))
        self.assertTrue(allclose(cis.solutesdata['reference']['pert_bonds_k'], bk, 0.01))
        self.assertTrue(allclose(cis.solutesdata['pert_delta_bonds'], bd, 0.0001))

    def test_pert_angle_pars(self):
        aa = [
            112.4,  121.4,  118.03,  118.03,  118.03,  118.03,  118.03,  118.03,  118.03,
            118.03,  118.03,  118.03,  118.03,  118.03,  118.03,  118.03,  118.03,  118.03,
            118.03,  118.03,  118.03
        ]
        ak = [
            63,  80,  70.96,  70.96,  70.96,  70.96,  70.96,  70.96,  70.96,  70.96,  70.96,
            70.96,  70.96,  70.96,  70.96,  70.96,  70.96,  70.96,  70.96,  70.96,  70.96
        ]
        ad = [
            3.753,  3.331,  3.536,  3.536,  3.536,  3.536,  3.536,  3.536,  3.536,  3.536,
            3.536,  3.536,  3.536,  3.536,  3.536,  3.536,  3.536,  3.536,  3.536, 
            3.536,  3.536
        ]
        cis = InitSolutes(**self.kws)
        cis.get_xyzs()
        cis.get_pars()
        cis.calc_neighbors()
        cis.calc_pert_angle_pars()
        self.assertTrue(allclose(cis.solutesdata['reference']['pert_angles_a'], aa, 0.01))
        self.assertTrue(allclose(cis.solutesdata['reference']['pert_angles_k'], ak, 0.01))
        self.assertTrue(allclose(cis.solutesdata['pert_delta_angles'], ad, 0.01))

    def test_pert_dihedral_pars(self):
        dv = [
            [4.669,  5.124,  0,     0],
            [0,      0.5,    0,     0],
            [0.7,   -1.5,    0,     0],
            [0.7,   -1.5,    0,     0],
            [0.7,   -1.5,    0,     0],
            [0.7,   -1.5,    0,     0],
            [0,      0.5,    0,     0],
            [0.7,   -1.5,    0,     0],
            [0.8,   -0.76,   0,     0],
            [0.8,   -0.76,   0,     0],
            [0,      0,      0.3,   0],
            [0.8,   -0.76,   0,     0],
            [0.8,   -0.76,   0,     0],
            [0.8,   -0.76,   0,     0],
            [0.8,   -0.76,   0,     0],
            [0.8,   -0.76,   0,     0],
            [0,      0,      0,     0],
            [0.8,   -0.76,   0,     0],
            [0.8,   -0.76,   0,     0],
            [0.8,   -0.76,   0,     0],
        ]
        dp = [[0.0,0.0,0.0] for i in range(len(dv))]
        dd = [2, 15, 15, 5, 5, 5, 5, 5, 15, 15, 15, 10, 10, 5, 5, 5, 5, 10, 10, 10]
        cis = InitSolutes(**self.kws)
        cis.get_xyzs()
        cis.get_pars()
        cis.calc_neighbors()
        cis.calc_pert_dihedral_pars()
        self.assertTrue(allclose(cis.solutesdata['reference']['pert_dihedrals_v'], dv, self.tol))
        self.assertTrue(allclose(cis.solutesdata['reference']['pert_dihedrals_p'], dp, self.tol))
        self.assertTrue(allclose(cis.solutesdata['pert_delta_dihedrals'], dd, self.tol))

    def test_combo(self):
        before_centroid_xyz1 = [
            [0.000 ,    0.000 ,    0.000    ],
            [1.000 ,    0.000 ,    0.000    ],
            [0.000 ,    1.000 ,    0.000    ],
            [-1.272 ,    1.000 ,    0.000   ],
            [-1.957 ,    3.112 ,    0.000   ],
            [-1.320 ,    3.599 ,    -1.112  ],
            [0.109 ,    3.747 ,    -1.018   ],
            [0.788 ,    1.656 ,    -1.080   ],
            [2.224 ,    1.581 ,    -0.960   ],
            [0.730 ,    0.455 ,    1.142    ],
            [2.794 ,    0.940 ,    0.091    ],
            [2.085 ,    0.428 ,    1.139    ],
            [4.269 ,    0.787 ,    0.134    ],
            [-1.547 ,    3.244 ,    1.011   ],
            [-3.029 ,    2.894 ,    -0.028  ],
            [-1.801 ,    3.638 ,    -2.099  ],
            [0.673 ,    4.081 ,    -1.892   ],
            [0.563 ,    3.979 ,    -0.049   ],
            [0.361 ,    1.631 ,    -2.102   ],
            [2.826 ,    2.023 ,    -1.770   ],
            [2.648 ,    0.001 ,    1.985    ],
            [0.159 ,    0.044 ,    1.991    ],
            [4.692 ,    1.268 ,    1.026    ],
            [4.771 ,    1.223 ,    -0.740   ],
            [4.553 ,    -0.272 ,    0.179   ]
        ]
        before_centroid_xyz2 = [
            [0.000 ,    0.000 ,    0.000    ],
            [1.000 ,    0.000 ,    0.000    ],
            [0.000 ,    1.000 ,    0.000    ],
            [-1.272 ,    1.000 ,    0.000   ],
            [-1.951 ,    3.093 ,    0.000   ],
            [-1.314 ,    3.580 ,    -1.112  ],
            [0.115 ,    3.728 ,    -1.018   ],
            [0.794 ,    1.637 ,    -1.080   ],
            [2.230 ,    1.561 ,    -0.960   ],
            [0.721 ,    0.466 ,    1.153    ],
            [2.800 ,    0.933 ,    0.098    ],
            [2.076 ,    0.435 ,    1.159    ],
            [4.275 ,    0.776 ,    0.154    ],
            [-1.541 ,    3.225 ,    1.011   ],
            [-3.023 ,    2.875 ,    -0.028  ],
            [-1.795 ,    3.619 ,    -2.099  ],
            [0.679 ,    4.062 ,    -1.892   ],
            [0.569 ,    3.960 ,    -0.049   ],
            [0.372 ,    1.601 ,    -2.104   ],
            [2.832 ,    1.994 ,    -1.776   ],
            [2.633 ,    0.016 ,    2.012    ],
            [0.143 ,    0.066 ,    2.002    ],
            [4.693 ,    1.269 ,    1.042    ],
            [4.784 ,    1.197 ,    -0.724   ],
            [4.555 ,    -0.283 ,    0.217   ]
        ]
        before_centroid_xyz3 = [
            [0.000 ,    0.000 ,    0.000    ],
            [1.000 ,    0.000 ,    0.000    ],
            [0.000 ,    1.000 ,    0.000    ],
            [-1.272 ,    1.000 ,    0.000   ],
            [-1.964 ,    3.131 ,    0.000   ],
            [-1.326 ,    3.618 ,    -1.112  ],
            [0.103 ,    3.766 ,    -1.018   ],
            [0.782 ,    1.675 ,    -1.080   ],
            [2.217 ,    1.600 ,    -0.960   ],
            [0.739 ,    0.445 ,    1.132    ],
            [2.788 ,    0.947 ,    0.083    ],
            [2.094 ,    0.422 ,    1.120    ],
            [4.264 ,    0.798 ,    0.115    ],
            [-1.554 ,    3.263 ,    1.011   ],
            [-3.035 ,    2.913 ,    -0.028  ],
            [-1.807 ,    3.657 ,    -2.099  ],
            [0.667 ,    4.100 ,    -1.892   ],
            [0.557 ,    3.998 ,    -0.049   ],
            [0.349 ,    1.660 ,    -2.100   ],
            [2.820 ,    2.053 ,    -1.764   ],
            [2.663 ,    -0.012 ,    1.957   ],
            [0.174 ,    0.023 ,    1.979    ],
            [4.691 ,    1.267 ,    1.011    ],
            [4.758 ,    1.249 ,    -0.757   ],
            [4.552 ,    -0.261 ,    0.142   ]
        ]
        beg_xyz1 = [
            [0.000 ,    -1.000 ,    0.000  ],
            [1.000 ,    -1.000 ,    0.000  ],
            [0.000 ,    0.000 ,    0.000   ],
            [-1.272 ,    0.000 ,    0.000  ],
            [-1.957 ,    2.112 ,    0.000  ],
            [-1.320 ,    2.599 ,    -1.112 ],
            [0.109 ,    2.747 ,    -1.018  ],
            [0.788 ,    0.656 ,    -1.080  ],
            [2.224 ,    0.581 ,    -0.960  ],
            [0.730 ,    -0.545 ,    1.142  ],
            [2.794 ,    -0.060 ,    0.091  ],
            [2.085 ,    -0.572 ,    1.139  ],
            [4.269 ,    -0.213 ,    0.134  ],
            [-1.547 ,    2.244 ,    1.011  ],
            [-3.029 ,    1.894 ,    -0.028 ],
            [-1.801 ,    2.638 ,    -2.099 ],
            [0.673 ,    3.081 ,    -1.892  ],
            [0.563 ,    2.979 ,    -0.049  ],
            [0.361 ,    0.631 ,    -2.102  ],
            [2.826 ,    1.023 ,    -1.770  ],
            [2.648 ,    -0.999 ,    1.985  ],
            [0.159 ,    -0.956 ,    1.991  ],
            [4.692 ,    0.268 ,    1.026   ],
            [4.771 ,    0.223 ,    -0.740  ],
            [4.553 ,    -1.272 ,    0.179  ]
        ]
        beg_xyz2 = [
            [-0.001 ,    -0.991 ,    -0.014], 
            [0.999 ,    -0.991 ,    -0.016 ],
            [-0.002 ,    0.009 ,    -0.008 ],
            [-1.274 ,    0.008 ,    -0.005 ],
            [-1.954 ,    2.100 ,    0.008  ],
            [-1.319 ,    2.594 ,    -1.103 ],
            [0.110 ,    2.743 ,    -1.011  ],
            [0.790 ,    0.652 ,    -1.086  ],
            [2.226 ,    0.576 ,    -0.970  ],
            [0.722 ,    -0.532 ,    1.140  ],
            [2.799 ,    -0.058 ,    0.084  ],
            [2.077 ,    -0.562 ,    1.143  ],
            [4.273 ,    -0.214 ,    0.136  ],
            [-1.542 ,    2.227 ,    1.018  ],
            [-3.025 ,    1.882 ,    -0.019 ],
            [-1.802 ,    2.638 ,    -2.089 ],
            [0.673 ,    3.081 ,    -1.884  ],
            [0.566 ,    2.969 ,    -0.041  ],
            [0.366 ,    0.622 ,    -2.110  ],
            [2.826 ,    1.013 ,    -1.784  ],
            [2.635 ,    -0.986 ,    1.993  ],
            [0.146 ,    -0.937 ,    1.989  ],
            [4.693 ,    0.274 ,    1.026   ],
            [4.780 ,    0.212 ,    -0.740  ],
            [4.554 ,    -1.274 ,    0.192  ]
        ]
        beg_xyz3 = [
            [0.001 ,    -1.009 ,    0.014  ],
            [1.001 ,    -1.009 ,    0.016  ],
            [0.002 ,    -0.009 ,    0.008  ],
            [-1.271 ,    -0.008 ,    0.005 ],
            [-1.961 ,    2.123 ,    -0.008 ],
            [-1.321 ,    2.603 ,    -1.121 ],
            [0.108 ,    2.752 ,    -1.025  ],
            [0.786 ,    0.660 ,    -1.074  ],
            [2.221 ,    0.585 ,    -0.951  ],
            [0.738 ,    -0.558 ,    1.144  ],
            [2.789 ,    -0.063 ,    0.097  ],
            [2.093 ,    -0.582 ,    1.136  ],
            [4.265 ,    -0.212 ,    0.133  ],
            [-1.553 ,    2.261 ,    1.003  ],
            [-3.033 ,    1.905 ,    -0.037 ],
            [-1.800 ,    2.638 ,    -2.109 ],
            [0.674 ,    3.081 ,    -1.899  ],
            [0.560 ,    2.989 ,    -0.056  ],
            [0.356 ,    0.640 ,    -2.095  ],
            [2.825 ,    1.033 ,    -1.756  ],
            [2.660 ,    -1.012 ,    1.976  ],
            [0.171 ,    -0.975 ,    1.993  ],
            [4.690 ,    0.262 ,    1.027   ],
            [4.761 ,    0.233 ,    -0.740  ],
            [4.553 ,    -1.271 ,    0.167  ]
        ]
        end_xyz1 = [
            [0.000 ,    -0.999 ,    -0.046  ],
            [-0.153 ,    -1.044 ,    0.941  ],
            [0.000 ,    0.000 ,    0.000    ],
            [0.195 ,    0.058 ,    -1.256   ],
            [0.300 ,    2.199 ,    -1.835   ],
            [1.301 ,    2.648 ,    -1.013   ],
            [0.989 ,    2.732 ,    0.390    ],
            [0.946 ,    0.611 ,    0.974    ],
            [0.608 ,    0.472 ,    2.369    ],
            [-1.241 ,    -0.569 ,    0.520  ],
            [-0.518 ,    -0.187 ,    2.741  ],
            [-1.446 ,    -0.658 ,    1.857  ],
            [-0.787 ,    -0.406 ,    4.184  ],
            [-0.761 ,    2.319 ,    -1.579  ],
            [0.492 ,    2.029 ,    -2.899   ],
            [2.350 ,    2.702 ,    -1.335   ],
            [1.766 ,    3.034 ,    1.096    ],
            [-0.038 ,    2.950 ,    0.701   ],
            [2.022 ,    0.599 ,    0.707    ],
            [1.316 ,    0.881 ,    3.107    ],
            [-2.367 ,    -1.104 ,    2.264  ],
            [-1.992 ,    -0.948 ,    -0.192 ],
            [-1.734 ,    0.061 ,    4.486   ],
            [0.000 ,    0.000 ,    4.833    ],
            [-0.875 ,    -1.477 ,    4.409  ]
        ]
        end_xyz2 = [
            [0.014 ,    -0.990 ,    -0.045  ],
            [-0.138 ,    -1.035 ,    0.943  ],
            [0.008 ,    0.009 ,    0.000    ],
            [0.201 ,    0.066 ,    -1.256   ],
            [0.292 ,    2.187 ,    -1.833   ],
            [1.292 ,    2.643 ,    -1.014   ],
            [0.982 ,    2.727 ,    0.390    ],
            [0.952 ,    0.607 ,    0.976    ],
            [0.617 ,    0.468 ,    2.372    ],
            [-1.237 ,    -0.556 ,    0.514  ],
            [-0.512 ,    -0.184 ,    2.747  ],
            [-1.448 ,    -0.648 ,    1.849  ],
            [-0.789 ,    -0.407 ,    4.188  ],
            [-0.770 ,    2.302 ,    -1.575  ],
            [0.483 ,    2.018 ,    -2.897   ],
            [2.340 ,    2.703 ,    -1.337   ],
            [1.759 ,    3.034 ,    1.094    ],
            [-0.046 ,    2.940 ,    0.702   ],
            [2.029 ,    0.589 ,    0.713    ],
            [1.329 ,    0.871 ,    3.109    ],
            [-2.374 ,    -1.091 ,    2.251  ],
            [-1.988 ,    -0.928 ,    -0.203 ],
            [-1.733 ,    0.067 ,    4.488   ],
            [-0.001 ,    -0.011 ,    4.842  ],
            [-0.888 ,    -1.479 ,    4.408  ]
        ]
        end_xyz3 = [
            [-0.014 ,    -1.008 ,    -0.047 ],
            [-0.169 ,    -1.054 ,    0.939  ],
            [-0.008 ,    -0.009 ,    0.000  ],
            [0.189 ,    0.050 ,    -1.255   ],
            [0.308 ,    2.210 ,    -1.837   ],
            [1.310 ,    2.653 ,    -1.012   ],
            [0.996 ,    2.737 ,    0.390    ],
            [0.941 ,    0.616 ,    0.971    ],
            [0.599 ,    0.476 ,    2.365    ],
            [-1.244 ,    -0.583 ,    0.527  ],
            [-0.524 ,    -0.189 ,    2.735  ],
            [-1.443 ,    -0.668 ,    1.865  ],
            [-0.785 ,    -0.405 ,    4.180  ],
            [-0.753 ,    2.336 ,    -1.582  ],
            [0.501 ,    2.041 ,    -2.900   ],
            [2.360 ,    2.702 ,    -1.332   ],
            [1.773 ,    3.033 ,    1.098    ],
            [-0.031 ,    2.960 ,    0.699   ],
            [2.016 ,    0.608 ,    0.701    ],
            [1.302 ,    0.891 ,    3.106    ],
            [-2.361 ,    -1.118 ,    2.277  ],
            [-1.996 ,    -0.968 ,    -0.182 ],
            [-1.734 ,    0.055 ,    4.485   ],
            [0.001 ,    0.011 ,    4.824    ],
            [-0.863 ,    -1.475 ,    4.410  ]
        ]
        cis = InitSolutes(**self.kws)
        cis.get_xyzs()
        cis.get_pars()
        cis.calc_neighbors()
        cis.calc_pert_bond_pars()
        cis.calc_pert_angle_pars()
        cis.calc_pert_xyzs()

        self.assertTrue(allclose(cis.solutesdata['reference']['xyz'],before_centroid_xyz1,0.008))
        self.assertTrue(allclose(cis.solutesdata['first']['xyz'],before_centroid_xyz2,0.008))
        self.assertTrue(allclose(cis.solutesdata['second']['xyz'],before_centroid_xyz3,0.008))

        cis.maxovl = True
        cis.calc_maximum_overlap()
        cis.recentroid()

        self.assertTrue(allclose(cis.solutesdata['reference']['xyz'],beg_xyz1,0.008))
        self.assertTrue(allclose(cis.solutesdata['first']['xyz'],beg_xyz2,0.008))
        self.assertTrue(allclose(cis.solutesdata['second']['xyz'],beg_xyz3,0.008))

        cis.izlong = True
        cis.calc_onto_z_axis()

        self.assertTrue(allclose(cis.solutesdata['reference']['xyz'],end_xyz1,0.08))
        self.assertTrue(allclose(cis.solutesdata['first']['xyz'],end_xyz2,0.08))
        self.assertTrue(allclose(cis.solutesdata['second']['xyz'],end_xyz3,0.08))

        # expdipole = [3.34488, 3.31363, 3.37612]
        # expquadrupole = [
        #     #   xx      yy       zz       xy     xz      yz
        #     [11.9151, 14.7701, 13.6605, 4.8538, 0.0660, 1.0411],
        #     [11.9084, 14.5355, 13.7303, 4.7106, 0.0572, 1.0109],
        #     [11.9224, 15.0028, 13.5953, 4.9960, 0.0753, 1.0712],
        # ]
        # for i,k in enumerate(['reference','first','second']):
        #     xyz = cis.solutesdata[k]['xyz']
        #     ql = cis.solutesdata[k]['Q']
        #     dip, qud = calc_dipoles(xyz,ql)
        #     self.assertTrue(dip, expdipole[i], self.tol)
        #     self.assertTrue(qud, expdipole[i], 0.0002)

    def test_calc_nonbond_pars(self):
        sa = [
            0.0,  0.0,  1.992370135,  1.683693072,  1.992370135,  1.992370135,  1.992370135,
            1.992370135,  1.992370135,  1.992370135,  1.992370135,  1.992370135,  1.992370135,
            1.380628319,  1.380628319,  1.380628319,  1.380628319,  1.380628319,  1.380628319,
            1.380628319,  1.380628319,  1.380628319,  1.380628319,  1.380628319,  1.380628319
        ]
        a = [
            0.0,  0.0,  1059.129669,  601.1488002,  1059.129669,  1059.129669,  1059.129669,
            1059.129669,  1059.129669,  1059.129669,  1059.129669,  1059.129669,  1059.129669,
            76.77171911,  76.77171911,  76.77171911,  76.77171911,  76.77171911,  76.77171911,
            76.77171911,  76.77171911,  76.77171911,  76.77171911,  76.77171911,  76.77171911
        ]
        b = [
            0.0,  0.0,  23.67358744,  22.26477038,  23.67358744,  23.67358744,  23.67358744,
            23.67358744,  23.67358744,  23.67358744,  23.67358744,  23.67358744,  23.67358744,
            5.156985904,  5.156985904,  5.156985904,  5.156985904,  5.156985904,  5.156985904,
            5.156985904,  5.156985904,  5.156985904,  5.156985904,  5.156985904,  5.156985904
        ]
        cis = InitSolutes(**self.kws)
        cis.get_xyzs()
        cis.get_pars()
        cis.calc_neighbors()
        cis.calc_nonbond_pars()
        self.assertTrue(allclose(cis.solutesdata['reference']['A'],a,self.tol))
        self.assertTrue(allclose(cis.solutesdata['reference']['B'],b,self.tol))

        cis.calc_sasa()
        self.assertTrue(allclose(cis.solutesdata['radius_sa'],sa,self.tol))

    def test_calc_total_molweights(self):
        cis = InitSolutes(**self.kws)
        cis.get_xyzs()
        cis.get_pars()
        cis.calc_total_molweights()
        self.assertTrue(allclose(cis.total_molecular_weights,148.2042,0.0001))


if __name__ == '__main__':
    unittest.main()


