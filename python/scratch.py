import math
def z2p(z):
    return 0.5 * (1 + math.erf(z / math.sqrt(2)))

print(int((z2p(4)*100)+0.5))
