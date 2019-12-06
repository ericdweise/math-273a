import numpy as np
from scipy import sparse


def conjugate_gradient(a, b, x):
       r = (b - np.matmul(a, x));
       p = r;
       rsold = np.matmul(r, r);

       for i in range(len(b)):
           a_p = np.matmul(a, p);
           alpha = rsold / np.matmul(p, a_p);
           x = x + (alpha * p);
           r = r - (alpha * a_p);
           rsnew = np.matmul(r, r);
           if (np.sqrt(rsnew) < (10 ** -7)):
               break;
           p = r + ((rsnew / rsold) * p);
           rsold = rsnew;
       return x
