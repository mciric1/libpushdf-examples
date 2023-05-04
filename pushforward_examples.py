#!/usr/bin/env python3
#
#   Further examples for libpushdf
#

from libpushdf import PushforwardCalculator, PushforwardCalculatorMP
from gen_se import scattering_equations, worldsheet_associahedron

class PushforwardExample(PushforwardCalculatorMP):
    """
    Pushforward of
        \frac{1}{z1 z2} dz1 ^ dz2
    Through
        I = < a1 z1 + a2 z2 , a3 z2^2 - 1 >
    """

    def __init__(self):
        super().__init__(['z1', 'z2'], ['a1', 'a2', 'a3'], (0,1))

    def ideal(self, Z, A):
        z1, z2 = Z
        a1, a2, a3 = A

        return [
            a1 * z1 + a2 * z2,
            a3 * z2**2 - 1,
        ]

    def form(self, Z):
        z1, z2 = Z
        return 1 / (z1 * z2)

class PushforwardNewtonSquare(PushforwardCalculatorMP):
    """
    Pushforward of Newton Square form
    """

    def __init__(self, a, b):
        super().__init__(['X1', 'X2'], ['x', 'y'], (0,1))
        self.a, self.b = a, b

    def ideal(self, Z, A):
        X1, X2 = Z
        x, y = A
        return [
            (X1 * X2 - 1) * y - X1 - self.a,
            (X1 * X2 - 1) * x - X2 - self.b,
        ]

    def form(self, Z):
        X1, X2 = Z
        return 1 / (X1 * X2)

class PushforwardPizzaSlice(PushforwardCalculatorMP):
    """
    Pushforward of the "pizza slice" form from Example 6.1 in
        "Positive Geometries and Canonical Forms" (https://arxiv.org/abs/1703.04541)
    """

    def __init__(self, z1, z2):
        super().__init__(['z', 't'], ['x', 'y'], (0,1))
        self.z1, self.z2 = z1, z2

    def ideal(self, Z, A):
        z, t = Z
        x, y = A
        return [
            (1 + z ** 2) * (1 + t) * x - 1 + z ** 2,
            (1 + z ** 2) * (1 + t) * y - 2 * z,
        ]

    def form(self, Z):
        z, t = Z
        return (self.z2 - self.z1) / ((self.z2 - z) * (z - self.z1) * t)

class PushforwardABHY(PushforwardCalculatorMP):
    """
    Pushforward of Worldsheet Associahedron forms through the Scattering Equations
    """

    def __init__(self, n : int):
        self.N = n
        nz = n - 3
        na = n * (n - 3) // 2
        super().__init__([f"z{i+2}" for i in range(nz)], [f"a{i+1}" for i in range(na)], tuple(range(nz)))

    def ideal(self, Z, A):
        return scattering_equations(self.N, Z[0].ring, A[0].ring)

    def form(self, Z):
        return worldsheet_associahedron(self.N, Z)

class PushforwardMomentumAmp4(PushforwardCalculatorMP):
    """
    Pushforward of n=4 k=2 Reduced Momentum Amplituhedron form
        \frac{1}{z1} dz1
    """

    def __init__(self):
        super().__init__(['z1', 'z2'], ['a1', 'a2'], (0,))

    def ideal(self, Z, A):
        z1, z2 = Z
        a1, a2 = A

        return [
            a1 - z1 * z2,
            a2 - z2
        ]

    def form(self, Z):
        z1, z2 = Z
        return 1 / z1

class PushforwardMomentumAmp_5_2(PushforwardCalculatorMP):
    """
    Pushforward of n=5 k=2 Reduced Momentum Amplituhedron form
        \frac{1}{z1 z2} dz1 ^ dz2
    """

    def __init__(self):
        super().__init__(['z1', 'z2', 'z3', 'z4', 'z5'], ['a1', 'a2', 'a3', 'a4', 'a5'], (0,1))

    def ideal(self, Z, A):
        z1, z2, z3, z4, z5 = Z
        a1, a2, a3, a4, a5 = A

        return [
            a1 - z1 * z2 * z3 - z1 * z3 - z1 * z2 * z4 + z1 * z4 * z5,
            a2 - z1 * z2 * z4,
            a3 - z3 + z4 * z5,
            a4 - z3 - z1 * z4 - z4,
            a5 + z1 * z4 * z5,
        ]

    def form(self, Z):
        z1, z2, z3, z4, z5 = Z
        return 1 / (z1 * z2)



if __name__ == '__main__':
    import time

    TIME = time.time()
    pf = PushforwardNewtonSquare(2, 4)
    results = pf.run()
    print(time.time() - TIME)
    pf.print_tex(results)

    TIME = time.time()
    pf = PushforwardPizzaSlice(1, -1)
    results = pf.run()
    print(time.time() - TIME)
    pf.print_tex(results)

    TIME = time.time()
    pf = PushforwardExample()
    results = pf.run()
    print(time.time() - TIME)
    pf.print_tex(results)

    TIME = time.time()
    pf = PushforwardMomentumAmp4()
    results = pf.run()
    print(time.time() - TIME)
    pf.print_tex(results)

    TIME = time.time()
    pf = PushforwardMomentumAmp_5_2()
    results = pf.run()
    print(time.time() - TIME)
    pf.print_tex(results)

    TIME = time.time()
    pf = PushforwardABHY(4)
    results = pf.run()
    print(time.time() - TIME)
    pf.print_tex(results)

    TIME = time.time()
    pf = PushforwardABHY(5)
    results = pf.run()
    print(time.time() - TIME)
    pf.print_tex(results)

    TIME = time.time()
    pf = PushforwardABHY(6)
    results = pf.run()
    print(time.time() - TIME)
    pf.print_tex(results)

    """
    TIME = time.time()
    pf = PushforwardABHY(7)
    results = pf.run()
    print(time.time() - TIME)
    pf.print_tex(results)
    """
