import numpy as np

# input subroutine: read from data file: MM, tend, factor, dtout, D, a, b
def INPUT(filename):
    data = np.loadtxt(filename)
    return (data)

# output subroutine:: print data to file
def OUTPUT(filename, x, U):
    # print x(i), U(i) into output file
    np.savetxt(filename, np.transpose([x, U]))
