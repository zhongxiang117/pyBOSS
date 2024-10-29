from PyBOSS.init_pars import get_solvents_data
from PyBOSS.utils import read_ff_file

from PyBOSS.Tests.utils import allclose
from PyBOSS.Tests import *

import unittest


class Test_get_solvents_data(unittest.TestCase):
    def test_data(self):
        A1 = [774.598765  ,    0          ,    0          ,    0            ]
        B1 = [24.69655193 ,    0          ,    0          ,    0            ]
        Q1 = [0           ,    9.475707045,    9.475707045,    -18.95141409 ]
        A2 = [690.3744484 ,    0          ,    2633.400442,    0            ]
        B2 = [23.85995294 ,    0          ,    48.95150361,    0            ]
        Q2 = [-12.75575948,    7.926793393,    4.82896609 ,    0            ]
        QAB = [
            [0           ,   600003.2467,     609.9196774],
            [0           ,   0           ,    0          ],
            [0           ,   0           ,    0          ],
            [0           ,   0           ,    0          ],
            [0           ,   0           ,    0          ],
            [89.789024   ,   0           ,    0          ],
            [89.789024   ,   0           ,    0          ],
            [-179.578048 ,   0           ,    0          ],
            [0           ,   0           ,    0          ],
            [89.789024   ,   0           ,    0          ],
            [89.789024   ,   0           ,    0          ],
            [-179.578048 ,   0           ,    0          ],
            [0           ,   0           ,    0          ],
            [-179.578048 ,   0           ,    0          ],
            [-179.578048 ,   0           ,    0          ],
            [359.156096  ,   0           ,    0          ],
        ]
        bond_pars, _ = read_ff_file(FILE_PAR,msg=False)
        pars = get_solvents_data('tip4p', 'ch3oh', bond_pars)
        tol = 0.00001
        self.assertTrue(allclose(pars['1']['AW'], A1, tol))
        self.assertTrue(allclose(pars['1']['BW'], B1, tol))
        self.assertTrue(allclose(pars['1']['QW'], Q1, tol))
        self.assertTrue(allclose(pars['2']['AW'], A2, tol))
        self.assertTrue(allclose(pars['2']['BW'], B2, tol))
        self.assertTrue(allclose(pars['2']['QW'], Q2, tol))
        self.assertTrue(
            allclose([pars['1-1'][f'{i}-{j}']
                for i in range(4) for j in range(4)], QAB, 0.0001)
        )


if __name__ == '__main__':
    unittest.main()


