Answer to question2:

2.1
/*
double A_ip = a[i+p*m];
double B_pj = b[p+j*k];
double C_ij = c[i+j*m];
C_ij = C_ij + A_ip * B_pj;
c[i+j*m] = C_ij;
*/

The whole block contains 10 float calculation and 4 memory access.


flops = 10 * m * n * k * NREPEATS / (time * 1e9)   (G/s)

bandwidth = 4 * 8 * m * n * k * NREPEATS  / (time * 1e9) (GB / s)

2.2
The optimal value for BLOCK_SIZE is the BLOCK_SIZE that allows the L1 cache stores three blocks
in the it.

We need to have
  BLOCK_SIZE * BLOCK_SIZE * 3 * 8 <= L1 cache capacity

In my laptop, L1 cache capacity = 4K

so BLOCK_SIZE <= 13

so we can choose BLOCK_SIZE = 12 in my laptop to get the best result.


