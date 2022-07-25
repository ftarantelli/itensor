#! /usr/bin/env wolframscript

mat = { {0, 0, 0, 2*d*Sin[k]}, {0, -m  - 2*Cos[k], 0, 0}, {0, 0, - m  - 2*Cos[k], 0}, {2*d*Sin[k], 0, 0, -2*m - 4*Cos[k]} }


Export["out.dat", Eigenvectors[mat] ]
Export["out.dat", Eigenvalues[mat] ]

F[k_, m_, d_] :=  Eigenvalues[{ {0, 0, 0, 2*d*Sin[k]}, {0, -m  - 2*Cos[k], 0, 0}, {0, 0, - m  - 2*Cos[k], 0}, {2*d*Sin[k], 0, 0, -2*m - 4*Cos[k]} }]
G[k_, m_, d_] :=  Eigenvectors[{ {0, 0, 0, 2*d*Sin[k]}, {0, -m  - 2*Cos[k], 0, 0}, {0, 0, - m  - 2*Cos[k], 0}, {2*d*Sin[k], 0, 0, -2*m - 4*Cos[k]} }]

Export["out.dat", G[0.06283185307179587,-2., 1.] ]

Export["out.dat", F[0.06283185307179587,-2., 1.] ]
