import numpy as np

from math import sqrt, fabs, acos, acosh, isnan, atan, atan2, sin, cos, pi
from cmath import sqrt as sqrt_complex

import pymeeus

from pymeeus.JupiterMoons import JupiterMoons
from pymeeus.Jupiter import Jupiter
from pymeeus.Earth import Earth

import plotly.graph_objects as go
import plotly.express as px

import plotly.io as pio
pio.renderers.default = "browser"

import datetime

class solverClass :

    def solveP3(self, a, b, c):
        
        """
        
        Solves a cubic equation algebraically, can also tolerate complex results

        :param a: coefficient
        :type a: float, complex float
        :param b: coefficient
        :type b: float, complex float
        :param c: coefficient
        :type c: float, complex float

        :returns: array containing result, as well as number of nonzero solutions
        :rtype: float/complex-float numpy array, int
        
        """


        a2 = a ** 2
        q = (a2 - 3 * b) / 9
        r = (a * (2 * a2 - 9 * b) + 27 * c) / 54
        r2 = r ** 2
        q3 = q ** 3
        M_2PI = pi * 2
        eps = 1.0e-12

        x = np.zeros(3, dtype=np.double)

        if r2 < q3:

            t = r / sqrt(q3)

            if t < -1:
                t = -1

            if t > 1:
                t = 1

            t = acos(t)
            a = a / 3
            q = -2 * sqrt(q)

            x[0] = q * cos(t / 3) - a
            x[1] = q * cos((t + M_2PI) / 3) - a
            x[2] = q * cos((t - M_2PI) / 3) - a

            return x, 3

        else:

            A = -pow(fabs(r) + sqrt(r2 - q3), 1 / 3)
            if r < 0:
                A = -A

            if A == 0:

                B = 0
            else:
                B = q / A

            a = a / 3
            x[0] = (A + B) - a
            x[1] = -0.5 * (A + B) - a
            x[2] = 0.5 * sqrt(3.) * (A - B)
            if fabs(x[2]) < eps:
                x[2] = x[1]
                return x, 2

            return x, 1

    def solve_quartic(self, a, b, c, d):

        """
        Solves a quartic equation algebraically, can also tolerate complex results

        :param a: coefficient
        :type a: float, complex float
        :param b: coefficient
        :type b: float, complex float
        :param c: coefficient
        :type c: float, complex float
        :param d: coefficient
        :type d: float, complex float

        :returns: array containing result, as well as number of nonzero solutions
        :rtype: foat/complex-float numpy array, int
        """


        eps = 1 * 10 ** -5

        num_real = 0
        retval = np.zeros(4, dtype=np.complex)
        x = np.zeros(4, dtype=np.double)

        A = b - 3.0 / 8.0 * a * a
        B = c - a * b / 2.0 + a * a * a / 8.0
        C = d - a * c / 4.0 + a * a * b / 16.0 - 3 * a * a * a * a / 256

        coef3 = np.zeros(3)
        coef3[0] = -A / 2.0
        coef3[1] = -C
        coef3[2] = (4 * A * C - B * B) / 8

        x3, iZeroes = self.solveP3(coef3[0], coef3[1], coef3[2])

        s = x3[0]

        # The essence - choosing Y with maximal absolute value.

        if iZeroes != 1:

            if fabs(x3[1]) > fabs(s):
                s = x3[1]
            if fabs(x3[2]) > fabs(s):
                s = x3[2]

        # y^4 + Ay^2 + C = 0
        if fabs(B) < eps:
            D = A * A - 4 * C
            x1 = (-A + sqrt_complex(D)) * 0.5
            x2 = (-A - sqrt_complex(D)) * 0.5

            # no abs originally
            retval[0] = sqrt_complex(x1) - a / 4
            retval[1] = -sqrt_complex(x1) - a / 4
            retval[2] = sqrt_complex(x2) - a / 4
            retval[3] = -sqrt_complex(x2) - a / 4



        else:

            sqD1 = np.lib.scimath.sqrt(2 * s - A)

            sqD2 = np.lib.scimath.sqrt(-2 * s - A + 2 * B / sqD1)

            sqD3 = np.lib.scimath.sqrt(-2 * s - A - 2 * B / sqD1)

            retval[0] = -0.5 * sqD1 + 0.5 * sqD2 - a / 4.0
            retval[1] = -0.5 * sqD1 - 0.5 * sqD2 - a / 4.0
            retval[2] = 0.5 * sqD1 + 0.5 * sqD3 - a / 4.0
            retval[3] = 0.5 * sqD1 - 0.5 * sqD3 - a / 4.0

        for i in range(0, 4):

            if fabs(np.imag(retval[i])) < eps:
                x[num_real] = np.real(retval[i]) #np.real
                num_real = num_real + 1

        return x, num_real

    def area_hyperbola(self, Q1, Q2, Os, Oe, Rs, a, b, in_out, x_axis):

        """
        the projection of the shadow cone might be an hyperbola, in that case compute the area of the hyperbola

        :param Q1: coefficient
        :type Q1: float, complex float
        :param Q2: coefficient
        :type Q2: float, complex float
        :param Os: position of the sun, real vector [3,1]
        :type Os: float, real, vector
        :param Oe: position of the shader, real vector [3,1]
        :type Oe:float, real, vector
        :param Rs: Radius of the sun (spherical model)
        :type Rs:float, real

        :returns: area of the section of the hyperbola
        :rtype: float
        """

        area = 0.0
        a_s = np.array([Q1[0] - Os[0], Q1[1] - Os[1]])
        b_s = np.array([Q2[0] - Os[0], Q2[1] - Os[1]])

        a_e = np.array([Q1[0] - Oe[0], Q1[1] - Oe[1]])
        b_e = np.array([Q2[0] - Oe[0], Q2[1] - Oe[1]])

        TQ1Q2Os = 0.5 * fabs(a_s[0] * b_s[1] - b_s[0] * a_s[1])
        TQ1Q2Oe = 0.5 * fabs(a_e[0] * b_e[1] - b_e[0] * a_e[1])

        cts = (a_s[0] * b_s[0] + a_s[1] * b_s[1]) / sqrt(
            (a_s[0] * a_s[0] + a_s[1] * a_s[1]) * (b_s[0] * b_s[0] + b_s[1] * b_s[1]))
        SQ1Q2Os = 0.5 * acos(cts) * Rs * Rs
        # print("cts : " + str(cts))

        if np.abs(cts) > 1:
            print("Math domain error")
            print("Cannot calculate the area of the hyperbolic section, now dropping out of the computation")
            exit(0)

        # calculate the area of hyperbolic section

        if not x_axis:
            tmp = 0.0
            tmp = Q1[0]

            Q1[0] = Q1[1]
            Q1[1] = -tmp

            tmp = Q2[0]
            Q2[0] = Q2[1]
            Q2[1] = -tmp

            tmp = a
            a = b
            b = tmp

        # ref: http:#www.rob_ertobigoni.eu/Matematica/Conics/segmentHyp/segmentHyp.html
        xQ1 = fabs(Q1[0])
        xQ2 = fabs(Q2[0])

        s1 = b * (xQ1 * sqrt_complex((xQ1 * xQ1 / a / a) - 1.0) - a * acosh(xQ1 / a))

        s2 = b * (xQ2 * sqrt_complex((xQ2 * xQ2 / a / a) - 1.0) - a * acosh(xQ2 / a))

        if isnan(s1) or isnan(s2):
            testc = 0

        if s1 <= s2:
            ss = s1
        else:
            ss = s2

        if s1 >= s2:
            sl = s1
        else:
            sl = s2

        S2 = 0.0
        if Q1[1] * Q2[1] < 0.0:

            x_m = 0.0

            if fabs(Q2[0] - Q1[0]) > 1.0E-10:

                k = (Q2[1] - Q1[1]) / (Q2[0] - Q1[0])
                m = Q1[1] - k * Q1[0]
                x_m = - m / k

            else:

                x_m = Q1[0]

            if Q1[0] > Q2[0]:

                S2 = ss + (sl - ss) / 2.0 - 0.5 * fabs((x_m - Q1[0]) * Q1[1]) + 0.5 * fabs((Q2[0] - x_m) * Q2[1])

            else:

                S2 = ss + (sl - ss) / 2.0 - 0.5 * fabs((Q2[0] - x_m) * Q2[1]) + 0.5 * fabs((x_m - Q1[0]) * Q1[1])




        else:

            s_trapezium = 0.5 * fabs(Q1[0] - Q2[0]) * (fabs(Q1[1]) + fabs(Q2[1]))
            S2 = (sl - ss - 2 * s_trapezium) / 2.0

        if in_out:

            area = SQ1Q2Os - TQ1Q2Os - S2

        else:  # the centre of the sun is outside the ellipse

            area_shadow = SQ1Q2Os - TQ1Q2Os + S2
            area = pi * Rs * Rs - area_shadow

        return area

    def area_ellipse(self, Q1, Q2, Os, Oe, Rs, a, b, in_out):

        """
        If the projection of the shadow cone is an ellipse, in that case compute the area of the ellipse

        :param Q1: coefficient
        :type Q1: float, complex float
        :param Q2: coefficient
        :type Q2: float, complex float
        :param Os: position of the sun, real vector [3,1]
        :type Os: float, real, vector
        :param Oe: position of the shader, real vector [3,1]
        :type Oe:float, real, vector
        :param Rs: Radius of the sun (spherical model)
        :type Rs:float, real

        :returns: area of the section of the hyperbola
        :rtype: float
        """


        area = 0.0

        Q1 = (Q1)
        Q2 = (Q2)
        Os = (Os)

        a_s = np.array([Q1[0] - Os[0], Q1[1] - Os[1]])
        b_s = np.array([Q2[0] - Os[0], Q2[1] - Os[1]])

        a_e = np.array([Q1[0], Q1[1]])
        b_e = np.array([Q2[0], Q2[1]])

        TQ1Q2Os = 0.5 * np.abs(a_s[0] * b_s[1] - b_s[0] * a_s[1])
        TQ1Q2Oe = 0.5 * np.abs(a_e[0] * b_e[1] - b_e[0] * a_e[1])

        cts = (a_s[0] * b_s[0] + a_s[1] * b_s[1]) / np.sqrt(
            (a_s[0] * a_s[0] + a_s[1] * a_s[1]) * (b_s[0] * b_s[0] + b_s[1] * b_s[1]))

        cts = (cts)

        # print("cts : " + str(cts))

        if np.abs(cts) > 1:
            print("Math domain error \n aborting execution with exit code -1")
            return (-1)

        SQ1Q2Os = 0.5 * acos(cts) * Rs ** 2

        SQ1Q2Oe = 0.0
        if a > b:
            aa = a
        else:
            aa = b

        q1 = 0.0
        q2 = 0.0
        if a_e[0] > 0 and a_e[1] > 0:  # the first quadrant

            q1 = atan(a_e[1] / a_e[0])

        elif a_e[0] < 0 and a_e[1] > 0:

            q1 = np.pi - atan(-a_e[1] / a_e[0])

        elif a_e[0] < 0 and a_e[1] < 0:

            q1 = np.pi + atan(a_e[1] / a_e[0])

        elif a_e[0] > 0 and a_e[1] < 0:

            q1 = np.pi * 2 - atan(-a_e[1] / a_e[0])

        if b_e[0] > 0 and b_e[1] > 0:  # the first quadrant

            q2 = atan(b_e[1] / b_e[0])

        elif b_e[0] < 0 and b_e[1] > 0:

            q2 = np.pi - atan(-b_e[1] / b_e[0])

        elif b_e[0] < 0 and b_e[1] < 0:

            q2 = np.pi + atan(b_e[1] / b_e[0])

        elif b_e[0] > 0 and b_e[1] < 0:

            q2 = np.pi * 2 - atan(-b_e[1] / b_e[0])

        s_a = a * b / 2.0 * (q1 - atan2((b - a) * sin(2.0 * q1), (a + b) + (b - a) * cos(2 * q1)))
        s_b = a * b / 2.0 * (q2 - atan2((b - a) * sin(2.0 * q2), (a + b) + (b - a) * cos(2 * q2)))

        SQ1Q2Oe = fabs(s_b - s_a)

        S1 = SQ1Q2Oe - TQ1Q2Oe

        if in_out:  # the centre of the sun is inside the ellipse

            area = SQ1Q2Os - TQ1Q2Os - S1

        else:  # the centre of the sun is outside the ellipse

            area_shadow = SQ1Q2Os - TQ1Q2Os + S1
            area = np.pi * Rs * Rs - area_shadow

        return area

class ShaderClass :

    """
    Class for objects casting a shadow

    methods :
        __init__
        perspectiveProjection()
        shadow_factor_parametric()
    """

    def __init__(self, polarRadius, equatorialRadius, heightAtmosphere, Solver):
        self.polarRadius = polarRadius
        self.equatorialRadius = equatorialRadius
        self.heightAtmosphere = heightAtmosphere
        self.Solver = Solver

    def perspectiveProjection(self, a, b, sunpos_ecef, satpos_ecef, r_solar, area_bright, dis_boundary, dis_circle, Solver):
        state = -1
        Rs = 695700.0  # Radius of the Sun in km , we assume the sun is a sphere in this model
        a2 = a ** 2
        b2 = b * b

        r = satpos_ecef
        rs = sunpos_ecef

        f = 1000.0
        n = r - rs

        dis_sat_earth = np.linalg.norm(r)
        dis_sat_sun = np.linalg.norm(n)
        n = n / np.linalg.norm(n)

        A = np.diag(np.diag(np.array([1 / a2, 1 / a2, 1 / b2])))

        nAin = n[0] ** 2 * a2 + n[1] ** 2 * a2 + + n[2] ** 2 * b2
        rtn = np.dot(r, n)
        ntn = np.dot(n, n)

        s1 = sqrt_complex(1.0 / nAin)
        s2 = - s1

        t1 = (rtn - 1.0 / s1) / ntn
        t2 = (rtn - 1.0 / s2) / ntn

        xs1 = r - t1 * n
        xs2 = r - t2 * n

        ds1 = np.linalg.norm((xs1 - rs))
        ds2 = np.linalg.norm((xs2 - rs))

        if ds1 <= ds2:
            t = ds1
        else:
            t = ds2

        if ds1 >= ds2:
            ds2 = ds1
        else:
            ds2 = ds2

        nAn = n[0] * n[0] / a2 + n[1] * n[1] / a2 + n[2] * n[2] / b2
        rsAn = n[0] * rs[0] / a2 + n[1] * rs[1] / a2 + n[2] * rs[2] / b2
        rsArs = rs[0] * rs[0] / a2 + rs[1] * rs[1] / a2 + rs[2] * rs[2] / b2
        Delta = rsAn * rsAn - nAn * (rsArs - 1.0)

        if Delta > 0:  # sun-sat line intersects the Earth ellipsoid

            t1 = (-2.0 * rsAn + sqrt_complex(Delta)) / 2.0 / nAn
            t2 = (-2.0 * rsAn - sqrt_complex(Delta)) / 2.0 / nAn

            if t1 <= t2:
                ds1 = t1
            else:
                ds1 = t2


        else:

            # normal vector to the plane sat,sun and Earth
            p = np.cross(r, n)
            # normal at the ellipsoid that is perpendicular to n
            ns = np.cross(n, p)
            ns = ns / np.linalg.norm(ns)

            t = ns[0] * a2 * ns[0] + ns[1] * a2 * ns[1] + ns[2] * b2 * ns[2]
            lam1 = sqrt(1.0 / t)
            lam2 = -lam1

            s1 = np.dot(ns, rs) - lam1 * t
            s2 = np.dot(ns, rs) - lam2 * t

            if fabs(s1) < fabs(s2):
                lam = lam1
            else:
                lam = lam2

            dis = lam * (n[0] * ns[0] * a2 + n[1] * ns[1] * a2 + n[2] * ns[2] * b2) - np.dot(n, rs)

            ds1 = dis

        if dis_sat_sun < ds1:  # full phase

            state = 0

            return state, r_solar, area_bright, dis_boundary, dis_circle

        elif dis_sat_sun >= ds2:  # ellipse

            state = 2

        elif not (not (dis_sat_sun >= ds1) or not (dis_sat_sun <= ds2)):  # hyperbola

            state = 1  # hyperbola

        # B: build the ISF coordinate system
        u = n - 1.0 / rtn * r

        u = u / np.linalg.norm(u)
        v = np.cross(n, u)

        # the projection of earth's centre PEC in ECEF
        t_dis = sqrt((f * dis_sat_earth / rtn) * (f * dis_sat_earth / rtn) - f * f)
        PEC = t_dis * u  # vector from PSC to PEC
        PEC_isf = np.array([0, 0, 0], dtype=np.double)
        # convert PEC into ISF
        PEC_isf[0] = np.dot(u, PEC)
        PEC_isf[1] = np.dot(v, PEC)
        PEC_isf[2] = np.dot(n, PEC)

        t = (r[0] * r[0] / a2 + r[1] * r[1] / a2 + r[2] * r[2] / b2 - 1.0)

        M = np.zeros([3, 3])
        M[0, 0] = r[0] * r[0] / a2 / a2 - t / a2
        M[0, 1] = r[0] * r[1] / a2 / a2
        M[0, 2] = r[0] * r[2] / a2 / b2
        M[1, 0] = M[0, 1]
        M[1, 1] = r[1] * r[1] / a2 / a2 - t / a2
        M[1, 2] = r[1] * r[2] / a2 / b2
        M[2, 0] = M[0, 2]
        M[2, 1] = M[1, 2]
        M[2, 2] = r[2] * r[2] / b2 / b2 - t / b2

        # D: calculate K
        K = np.zeros([6, 1], dtype=np.double)

        K[0] = (u[0] * M[0, 0] + u[1] * M[1, 0] + u[2] * M[2, 0]) * u[0] + (
                u[0] * M[0, 1] + u[1] * M[1, 1] + u[2] * M[2, 1]) * u[1] + (
                       u[0] * M[0, 2] + u[1] * M[1, 2] + u[2] * M[2, 2]) * u[2]

        K[1] = (u[0] * M[0, 0] + u[1] * M[1, 0] + u[2] * M[2, 0]) * v[0] + (
                u[0] * M[0, 1] + u[1] * M[1, 1] + u[2] * M[2, 1]) * v[1] + (
                       u[0] * M[0, 2] + u[1] * M[1, 2] + u[2] * M[2, 2]) * v[2]

        K[2] = (v[0] * M[0, 0] + v[1] * M[1, 0] + v[2] * M[2, 0]) * v[0] + (
                v[0] * M[0, 1] + v[1] * M[1, 1] + v[2] * M[2, 1]) * v[1] + (
                       v[0] * M[0, 2] + v[1] * M[1, 2] + v[2] * M[2, 2]) * v[2]

        K[3] = -((u[0] * M[0, 0] + u[1] * M[1, 0] + u[2] * M[2, 0]) * n[0] + (
                u[0] * M[0, 1] + u[1] * M[1, 1] + u[2] * M[2, 1]) * n[1] + (
                         u[0] * M[0, 2] + u[1] * M[1, 2] + u[2] * M[2, 2]) * n[2]) * f * 2.0

        K[4] = -((v[0] * M[0, 0] + v[1] * M[1, 0] + v[2] * M[2, 0]) * n[0] + (
                v[0] * M[0, 1] + v[1] * M[1, 1] + v[2] * M[2, 1]) * n[1] + (
                         v[0] * M[0, 2] + v[1] * M[1, 2] + v[2] * M[2, 2]) * n[2]) * f * 2.0

        K[5] = ((n[0] * M[0, 0] + n[1] * M[1, 0] + n[2] * M[2, 0]) * n[0] + (
                n[0] * M[0, 1] + n[1] * M[1, 1] + n[2] * M[2, 1]) * n[1] + (
                        n[0] * M[0, 2] + n[1] * M[1, 2] + n[2] * M[2, 2]) * n[2]) * f * f

        # just make the numbers larger for visualisation
        for i in range(6):
            K[i] = K[i] * 1.0E6

        # E: calculate the eigen value and eigen vector of matrix B
        t = sqrt((K[0] + K[2]) * (K[0] + K[2]) - 4.0 * (K[0] * K[2] - K[1] * K[1])) #sqrt
        lambda1 = (K[0] + K[2] + t) / 2.0
        lambda2 = (K[0] + K[2] - t) / 2.0

        r1 = np.array([0.0, 1.0])
        r2 = np.array([1.0, 0.0])

        if fabs(K[1]) < 1.0E-12:

            r1[0] = 1
            r1[1] = 0

            r2[0] = 0
            r2[1] = 1

        else:

            r1[0] = lambda1 - K[2]
            r1[1] = K[1]
            r2[0] = lambda2 - K[2]
            r2[1] = K[1]

        # get the unit vector
        t = sqrt(r1[0] * r1[0] + r1[1] * r1[1])
        r1[0] = r1[0] / t
        r1[1] = r1[1] / t

        t = sqrt(r2[0] * r2[0] + r2[1] * r2[1])
        r2[0] = r2[0] / t
        r2[1] = r2[1] / t

        OM = (K[2] * K[3] * K[3] - 2.0 * K[1] * K[3] * K[4] + K[0] * K[4] * K[4]) / (4.0 * (K[0] * K[2] - K[1] * K[1]))
        # the translation parameters
        tx = (r1[0] * K[3] + r1[1] * K[4]) / lambda1

        ty = (r2[0] * K[3] + r2[1] * K[4]) / lambda2

        # semi-major axis and the semi-minor axis
        AA = (OM - K[5]) / lambda1
        BB = (OM - K[5]) / lambda2
        R0 = f * Rs / dis_sat_sun
        r_solar = R0

        # F: Solve the quartic in the translated and rotated frame

        A = lambda1 * (R0 - 0.5 * tx) ** 2 + 0.25 * lambda2 * ty ** 2 - OM + K[5]
        B = 2 * lambda2 * R0 * ty
        C = lambda1 * (0.5 * tx ** 2 - 2 * R0 ** 2) + lambda2 * (0.5 * ty ** 2 + 4 * R0 ** 2) - 2 * OM + 2 * K[5]
        D = 2 * lambda2 * R0 * ty
        E = lambda1 * (R0 + 0.5 * tx) * (R0 + 0.5 * tx) + 0.25 * lambda2 * ty * ty - OM + K[5]

        ce = np.array([B / A, C / A, D / A, E / A], dtype=np.double)  # object

        XX, num_of_solution = Solver.solve_quartic(ce[0], ce[1], ce[2], ce[3])

        if num_of_solution == 3 or num_of_solution == 4:
            print("WARNING: shadowfactor, intersections -->" + str(num_of_solution))

        # calculate the sun's projection centre Os (PSC) in the transformed frame
        Os = np.array([0.5 * tx, 0.5 * ty])
        PEC_new = np.array([0, 0], dtype=np.double)
        # the projection of the earth's centre Pe in the transformed frame is obtained by converting PEC_isf into the rotated and translated frame
        PEC_new[0] = PEC_isf[0] * r1[0] + r1[1] * PEC_isf[1] + 0.5 * tx
        PEC_new[1] = PEC_isf[0] * r2[0] + r2[1] * PEC_isf[1] + 0.5 * ty

        # calculate the coordinates of the intersections between the circle and the conical curve in transformed frame
        # in the new frame, the origin is at the center of the projection of the Earth
        Q1 = np.array([(1.0 - XX[0] * XX[0]) / (1.0 + XX[0] * XX[0]) * R0 + 0.5 * tx,
                       2 * XX[0] / (1.0 + XX[0] * XX[0]) * R0 + 0.5 * ty])
        Q2 = np.array([(1.0 - XX[1] * XX[1]) / (1.0 + XX[1] * XX[1]) * R0 + 0.5 * tx,
                       2 * XX[1] / (1.0 + XX[1] * XX[1]) * R0 + 0.5 * ty])

        # figure out the intersection between the line from PEC_new to Os and the conical curve.
        # parameter equation of line from Oe to Os in transformed frame, x = td
        PEC_Os_len = sqrt((Os[0] - PEC_new[0]) * (Os[0] - PEC_new[0]) + (Os[1] - PEC_new[1]) * (Os[1] - PEC_new[1]))
        d = np.array([(Os[0] - PEC_new[0]) / PEC_Os_len, (Os[1] - PEC_new[1]) / PEC_Os_len])
        aa = lambda1 * d[0] * d[0] + lambda2 * d[1] * d[1]
        bb = 2 * (lambda1 * PEC_new[0] * d[0] + lambda2 * PEC_new[1] * d[1])
        cc = lambda1 * PEC_new[0] * PEC_new[0] + lambda2 * PEC_new[1] * PEC_new[1] + K[5] - OM
        t1 = (-bb + sqrt_complex(bb * bb - 4.0 * aa * cc)) / aa / 2.0
        t2 = (-bb - sqrt_complex(bb * bb - 4.0 * aa * cc)) / aa / 2.0

        if t1 * t2 < 0.0:

            if t1 > 0:
                t = t1
            else:
                t = t2

        else:

            if t1 <= t2:
                t = t1
            else:
                t = t2

        dis_boundary = t

        # figure out the intersection between the line from PEC_new to Os and the solar circle
        aa = d[0] * d[0] + d[1] * d[1]
        bb = 2 * (PEC_new[0] * d[0] + PEC_new[1] * d[1]) - (d[0] * tx + d[1] * ty)
        cc = PEC_new[0] * PEC_new[0] + PEC_new[1] * PEC_new[1] - (PEC_new[0] * tx + PEC_new[1] * ty) + 0.25 * (
                tx * tx + ty * ty) - R0 * R0

        t1 = (-bb + sqrt_complex(bb * bb - 4.0 * aa * cc)) / aa / 2.0
        t2 = (-bb - sqrt_complex(bb * bb - 4.0 * aa * cc)) / aa / 2.0

        # t should be the smaller one, which means closer to the solid Earth

        if t1 * t2 < 0.0:

            if t1 > 0:
                t = t1
            else:
                t = t2

        else:

            if t1 <= t2:
                t = t1
            else:
                t = t2

        dis_circle = t

        t = K[0] * K[2] - K[1] * K[1]

        if t > 0.0:  # ellipse

            if OM / (OM - K[5]) <= 1.0:

                in_out = True

            else:

                in_out = False

            if in_out == True and num_of_solution < 2:  # umbra

                state = -1
                return state, r_solar, area_bright, dis_boundary, dis_circle

            if in_out == False and num_of_solution < 2:  # full phase

                state = 0
                return state, r_solar, area_bright, dis_boundary, dis_circle

            _ = 0
            # penumbra
            area_bright = Solver.area_ellipse(Q1, Q2, Os, _, R0, sqrt(AA), sqrt(BB), in_out)


        elif t < 0.0:  # hyperbola

            if OM / (OM - K[5]) <= 1.0:

                in_out = False

            else:

                in_out = True

            if in_out == True and num_of_solution < 2:
                state = -1  # umbra
                return state, r_solar, area_bright, dis_boundary, dis_circle

            if in_out == False and num_of_solution < 2:
                state = 0  # full phase
                return state, r_solar, area_bright, dis_boundary, dis_circle

            # penumbra
            a = 0
            b = 0
            x_axis = True

            if AA > 0:
                a = sqrt(AA)

            else:
                a = sqrt(-AA)
                x_axis = False

            if BB > 0:

                b = sqrt(BB)
                x_axis = False

            else:

                b = sqrt(-BB)

            area_bright = Solver.area_hyperbola(Q1, Q2, Os, PEC_new, R0, a, b, in_out, x_axis)

        return state, r_solar, area_bright, dis_boundary, dis_circle

    def shadow_factor_parametric(self, sunpos_eci, satpos_eci,a, b, hgt_atm, Solver):

        # considering atmosphere
        a_atm = a + hgt_atm  # km
        b_atm = b * (a_atm / a)  # km  6356.7523142

        dis0 = 0.0
        dis1 = 0.0
        dis2 = 0.0

        r_sun = 0.0

        if hgt_atm == 0.0:

            Area_earth = 0.0

            state2, r_sun, Area_earth, dis1, dis0 = self.perspectiveProjection(a, b, sunpos_eci, satpos_eci, r_sun,
                                                                            Area_earth,
                                                                            dis1, dis0, Solver)
            Area_solar = np.pi * r_sun * r_sun

            if state2 == 0:

                factor = 1.0

            elif state2 == -1:

                factor = 0.0

            else:

                factor = Area_earth / Area_solar

        else:  # considering the atmosphere

            Area_earth = 0
            Area_atm = 0

            state1, r_sun, Area_atm, dis1, dis0 = self.perspectiveProjection(a_atm, b_atm, sunpos_eci, satpos_eci, r_sun,
                                                                          Area_atm, dis1, dis0, Solver)
            state2, r_sun, Area_earth, dis2, dis0 = self.perspectiveProjection(a, b, sunpos_eci, satpos_eci, r_sun,
                                                                            Area_earth,
                                                                            dis2, dis0, Solver)
            mu1 = 0.0  # distance = 0
            mu2 = 1.0  # distance = 1

            Area_solar = np.pi * r_sun ** 2

            # the thickness of atmosphere
            thickness = dis1 - dis2

            if state1 == 0:  # full phase

                factor = 1.0

            # partly in the atmosphere
            if (state1 == 1 or state1 == 2) and state2 == 0:
                x1 = 1.0
                x2 = (thickness - (dis1 - dis0)) / thickness

                u1 = (mu2 - mu1) * x1 + mu1
                u2 = (mu2 - mu1) * x2 + mu1

                factor = (Area_atm + (Area_solar - Area_atm) * ((u1 + u2) / 2.0)) / Area_solar

            # totally in the atmosphere
            if state1 == -1 and state2 == 0:

                u1 = 0
                u2 = 0
                if dis2 != 0.0:

                    x1 = (dis0 - dis2) / thickness
                    x2 = (dis0 - dis2 + 2.0 * r_sun) / thickness

                    u1 = (mu2 - mu1) * x1 + mu1
                    u2 = (mu2 - mu1) * x2 + mu1

                elif (
                        fabs(dis2) < 1.0E-9):  # no projection for the solid Earth (outside of the solid Earth), but the satellite is in the atmosphere's projection

                    print(
                        "no projection for the solid Earth (outside of the solid Earth), but the satellite is in the atmosphere's projection")

                    u1 = 0.5
                    u2 = 0.5

                factor = (u1 + u2) / 2.0

            # partially in the earth and atmosphere
            if state1 == -1 and (state2 == 1 or state2 == 2):
                x1 = (2.0 * r_sun - (dis2 - dis0)) / thickness
                x2 = 0.0

                u1 = (mu2 - mu1) * x1 + mu1
                u2 = (mu2 - mu1) * x2 + mu1

                factor = (Area_earth * ((u1 + u2) / 2.0)) / Area_solar

            # totally in the shadow, umbra
            if state2 == -1:
                factor = 0.0

            # the sun is bigger than the atmospheric thickness
            if (state1 == 1 or state1 == 2) and (state2 == 1 or state2 == 2):
                area1 = Area_atm
                area2 = Area_earth - Area_atm

                u1 = (mu2 - mu1) + mu1
                u2 = mu1

                factor = (area1 * 1.0 + area2 * (u1 + u2) / 2.0) / Area_solar

        return factor

class ShadedClass :
    """
    Class of objects getting a shadow cast on
    """
    def __init__(self, radius, n_points, height_atmosphere, delta, Shader, SelfShader):
        self.SelfShader = SelfShader
        self.Shader = Shader
        self.height_atmosphere = height_atmosphere
        self.radius = radius
        self.n_points = n_points
        self.phi_arr = np.linspace(0, 2 * np.pi, n_points)
        self.theta_arr = np.linspace(0, +np.pi, n_points)
        self.delta = delta

        self.phi_grid, self.theta_grid = np.meshgrid(self.phi_arr, self.theta_arr)

    def ComputeMutualShadow(self, sunpos, centerpos, Shader):

        shadowMesh = np.zeros([self.n_points, self.n_points])

        for i in range(self.n_points):

            for e in range(self.n_points):

                #position on the surface of the sphere

                x = self.radius * np.cos(self.phi_arr[i]) * np.sin(self.theta_arr[e])
                y = self.radius * np.sin(self.phi_arr[i]) * np.sin(self.theta_arr[e])
                z = self.radius * np.cos(self.theta_arr[e])

                surfacePos = np.array([x, y, z]) + centerpos

                shadowMesh[e, i] = np.real(Shader.shadow_factor_parametric(sunpos, surfacePos, Shader.polarRadius, Shader.equatorialRadius, Shader.heightAtmosphere,Shader.Solver))

        return shadowMesh

    def ComputeSelfShadow(self, sunpos, centerpos, SelfShader):

        selfShadowMesh = np.zeros([self.n_points, self.n_points])

        for i in range(self.n_points):

            for e in range(self.n_points):
                # position on the surface of the sphere

                x = (self.radius+self.delta) * np.cos(self.phi_arr[i]) * np.sin(self.theta_arr[e])
                y = (self.radius+self.delta) * np.sin(self.phi_arr[i]) * np.sin(self.theta_arr[e])
                z = (self.radius+self.delta) * np.cos(self.theta_arr[e])

                surfacePos = np.array([x, y, z]) + centerpos

                #print("surface pos : " + str(surfacePos - centerpos))
                #print(np.linalg.norm(surfacePos - centerpos))
                #print("Polar radius : " + str(SelfShader.polarRadius))
                #print("Polar radius : " + str(SelfShader.equatorialRadius))

                selfShadowMesh[e, i] = self.ComputeReflectedLight(sunpos, surfacePos, 1)#SelfShader.shadow_factor_parametric(sunpos, surfacePos, SelfShader.polarRadius, SelfShader.equatorialRadius, SelfShader.heightAtmosphere, SelfShader.Solver)

                print(selfShadowMesh[e, i])

                #print (str(e)+ " , " +str(i) + " : " + str(selfShadowMesh[e, i]))

        return selfShadowMesh

    def ComputeReflectedLight(self, centerpos, ObserverPos, albedo):

        ReflectivityMask = np.ones([self.n_points, self.n_points]) * albedo

        for i in range(self.n_points):

            for e in range(self.n_points):
                # position on the surface of the sphere

                x = (self.radius + self.delta) * np.cos(self.phi_arr[i]) * np.sin(self.theta_arr[e])
                y = (self.radius + self.delta) * np.sin(self.phi_arr[i]) * np.sin(self.theta_arr[e])
                z = (self.radius + self.delta) * np.cos(self.theta_arr[e])

                surfacePos = np.array([x, y, z]) + centerpos

                a= self.radius + self.delta
                b = np.linalg.norm(centerpos - ObserverPos)
                c = np.linalg.norm(surfacePos - ObserverPos)

                gamma = np.arccos((a**2+b**2-c**2)/(2*a*b))
                delta = np.pi - gamma

                ReflectivityMask[e, i] = ReflectivityMask[e, i] * np.cos(delta)

        return np.clip(ReflectivityMask, 0, 1)

    def FilterInstabilityOut(self, shadow_array):

        max_pos = np.shape(shadow_array)[2]-1
        filter_width = 6

        derivative = shadow_array[:, :, 0:max_pos-5] - shadow_array[:, :, 1:max_pos-4] - shadow_array[:, :, 2:max_pos-3] \
                     - shadow_array[:,:,3:max_pos-2] - shadow_array[:, :,4:max_pos-1] - shadow_array[:,:,5:max_pos]

        x_dim = np.shape(shadow_array)[0]
        y_dim = np.shape(shadow_array)[1]

        print(np.shape(derivative))

        filtered = np.zeros([x_dim, y_dim, max_pos-filter_width])

        for t in range(0, max_pos-filter_width-1):

            for i in range(x_dim):
                for a in range(y_dim):

                    if (derivative[i, a, t]) != 0:
                        filtered[i, a, t] = 1

        max_pos = np.shape(filtered)[2]-1

        derivative_2nd_deg = filtered[:, :, 0:max_pos-6] - filtered[:, :, 1:max_pos-5] - filtered[:, :, 2:max_pos-4] \
                           - filtered[:,:,3:max_pos-3] - filtered[:, :,4:max_pos-2] - filtered[:,:,5:max_pos-1]

        filtered_2nd_deg = np.zeros([x_dim, y_dim, max_pos-filter_width])

        for t in range(max_pos-filter_width):

            for i in range(x_dim):
                for a in range(y_dim):

                    if (derivative_2nd_deg[i, a, t]) != 0:
                        filtered_2nd_deg[i, a, t] = 1

        return filtered_2nd_deg #, 2*filter_width

class CoordinatesTransformator :
    """
    Convert all tp jovian equatorial
    """
    def __init__(self, Omega, i, I):
        self.Omega = Omega
        self.i = i
        self.I = I
    def HelioEcJovEq(self, JupiterPosHeliocentricEcliptic, TargetPosHeliocentricEcliptic):

        #translate system to jovian coordinatas
        TargetPosHeliocentricEcliptic  = TargetPosHeliocentricEcliptic - JupiterPosHeliocentricEcliptic

        #rotate away from the vernal equinox
        R_1 = np.array([[np.cos(-self.Omega), -np.sin(-self.Omega), 0],
                        [np.sin(-self.Omega),  np.cos(-self.Omega), 0],
                        [                  0,                    0, 1]])

        TargetPosHeliocentricEcliptic = TargetPosHeliocentricEcliptic @ R_1

        R_2 = np.array([[1,               0,                0],
                        [0, np.cos(-self.i), -np.sin(-self.i)],
                        [0, np.sin(-self.i), -np.cos(-self.i)]])

        TargetPosHeliocentricEcliptic = TargetPosHeliocentricEcliptic @ R_2

        R_3 = np.array([[1,               0,                0],
                        [0, np.cos(-self.I), -np.sin(-self.I)],
                        [0, np.sin(-self.I), -np.cos(-self.I)]])

        TargetPosHeliocentricEcliptic = TargetPosHeliocentricEcliptic @ R_3

        return TargetPosHeliocentricEcliptic

    def BesselianEq(self, BesselianPos, Psi, Omega):

        Phi = Psi - Omega

        R = np.array([[np.cos(Phi), -np.sin(Phi), 0],
                      [np.sin(Phi),  np.cos(Phi), 0],
                      [0,                      0, 1]])

        return BesselianPos@R

def test():

    # general simulation parameters, for this model to work all coordinates need to be in km
    s_jd = 1.157401129603386e-05
    au_to_km = 149597870700.0
    jupiter_eq_radius_to_km = 7.1492 * 10 ** 4
    n_points = 5

    dt = datetime.datetime(2021, 6, 7, 0, 0, 0)

    n_timesteps = 5*60*24 # 3 h worth of data
    delta = 10 ** -3  # km --> 1 m
    io_radius = 1821.6 # km
    europa_radius = 3121.6 # km
    ganymede_radius = 5268.2 # km
    callisto_radius = 4820.6 # km
    

    # set the epoch to the first of january 2020, this should be ok for this explicative test
    epoch = pymeeus.Epoch.Epoch()
    epoch.set(dt)
    #epoch = epoch + 5200 *s_jd

    # there are two other reference epochs we need to compute, to be able to calculate the angles I and i for the
    # coordinate transform we will be doing later on, in order to have all coordinates in teh jovian equatorial system
    ref_epoch_1900 = pymeeus.Epoch.Epoch()
    ref_epoch_1900.set(1900,1,1,0)

    ref_epoch_2000 = pymeeus.Epoch.Epoch()
    ref_epoch_2000.set(2000, 1, 1, 0)

    T_1900 = (epoch - ref_epoch_1900)/36525
    T_2000 = (epoch - ref_epoch_2000)/36525

    # declare solver, needed for solving the cubic and quartic equations, for the computation of the are
    defaultSolver = solverClass()

    # characterizing jupiter as an ellipsoid
    a_jupiter = 71492 # km , equatorial radius
    b_jupiter = 66854 # km , polar radius
    h_atm_jupiter = 100 # km
    jupiter = ShaderClass(a_jupiter, b_jupiter, h_atm_jupiter, defaultSolver)

    # get the position of jupiter and of the earth for the test epoch, those position are heliocentric ecliptic spherical
    satellite_positions = np.array(JupiterMoons.rectangular_positions_jovian_equatorial(epoch))
    jupiter_position = Jupiter.geometric_heliocentric_position(epoch)
    earth_position = Earth.geometric_heliocentric_position(epoch)

    # convert everything to cartesian and scale to km
    x_jupiter = jupiter_position[2] * np.cos(jupiter_position[1].rad()) * np.cos(jupiter_position[0].rad())
    y_jupiter = jupiter_position[2] * np.cos(jupiter_position[1].rad()) * np.sin(jupiter_position[0].rad())
    z_jupiter = jupiter_position[2] * np.sin(jupiter_position[1].rad())

    jupiter_pos_rectangular_heliocentric = np.array([x_jupiter, y_jupiter, z_jupiter])
    jupiter_pos_rectangular_heliocentric = jupiter_pos_rectangular_heliocentric * au_to_km

    x_earth = earth_position[2] * np.cos(earth_position[1].rad()) * np.cos(earth_position[0].rad())
    y_earth = earth_position[2] * np.cos(earth_position[1].rad()) * np.sin(earth_position[0].rad())
    z_earth = earth_position[2] * np.sin(earth_position[1].rad())

    earth_pos_rectangular_heliocentric = np.array([x_earth, y_earth, z_earth])
    earth_pos_rectangular_heliocentric = earth_pos_rectangular_heliocentric * au_to_km


    # the two angles needed for the coordinate conversion
    I = np.deg2rad(3.120262  + 0.0006*T_1900)
    i = np.deg2rad(1.303267 - 0.0054965*T_2000 + 0.00000466*T_2000**2 + 0.000000002*T_2000**3)
    Omega = np.deg2rad(100.464407 + 1.0209774*T_2000 + 0.00040315*T_2000**2 + 0.000000404*T_2000**3)

    ref_epoch_1976 = pymeeus.Epoch.Epoch()
    ref_epoch_1976.set(243000.5)
    t_1976 = epoch - ref_epoch_1976

    Psi = np.deg2rad(316.5182 - 0.0000208*t_1976)

    # now create an instance of the coordinate transformer object
    coordTrans = CoordinatesTransformator(Omega, i, I)

    # convert ad scale the coordinates to the target frame in km
    sun_pos_eq_jvc = coordTrans.HelioEcJovEq(jupiter_pos_rectangular_heliocentric, -earth_pos_rectangular_heliocentric)
    earth_pos_eq_jvc = coordTrans.HelioEcJovEq(jupiter_pos_rectangular_heliocentric, earth_pos_rectangular_heliocentric)
    satellite_positions = satellite_positions * jupiter_eq_radius_to_km

    # define io as shader and shaded
    io_shader = ShaderClass(io_radius, io_radius, delta, defaultSolver)
    io_shaded = ShadedClass(io_radius, n_points, 0, delta, jupiter, io_shader)

    # define europa as shader and shaded
    europa_shader = ShaderClass(europa_radius, europa_radius, delta, defaultSolver)
    europa_shaded = ShadedClass(europa_radius, n_points, 0, delta, jupiter, europa_shader)
    europa_shaded_from_io = ShadedClass(europa_radius, n_points, 0, delta, io_shader, europa_shader)

    # define ganymede as shader and shaded
    ganymede_shader = ShaderClass(ganymede_radius, ganymede_radius, delta, defaultSolver)
    ganymede_shaded = ShadedClass(ganymede_radius, n_points, 0, delta, jupiter, ganymede_shader)

    # define callisto as shader and shaded
    callisto_shader = ShaderClass(callisto_radius, callisto_radius, delta, defaultSolver)
    callisto_shaded = ShadedClass(callisto_radius, n_points, 0, delta, jupiter, callisto_shader)

    # initialize arrays
    io_shadow_array = np.zeros([n_points, n_points, n_timesteps])
    reflected_light_array = np.zeros([n_points, n_points, n_timesteps])

    # get the evolution of the positions for debugging purposes
    io_position_array = np.zeros([n_timesteps, 3])
    europa_position_array = np.zeros([n_timesteps, 3])
    diffpos_iocentric_jovian_eq = np.zeros([n_timesteps, 3])

    # compute the jovian shadow array and the quantity of reflected light towards earth
    # the following returns an array containing the shadow cast by io on the selected europa
    for t in range(0, n_timesteps):

        # set epochs and get coordinates
        epoch = epoch + s_jd*60  # advance in 1min steps
        satellite_positions = np.array(JupiterMoons.rectangular_positions_jovian_equatorial(epoch, solar=True))
        io_pos = satellite_positions[0, [0,2,1]] * jupiter_eq_radius_to_km
        #io_pos = coordTrans.BesselianEq(io_pos, Psi, Omega)
        europa_pos = satellite_positions[1, [0,2,1]] * jupiter_eq_radius_to_km
        #europa_pos = coordTrans.BesselianEq(europa_pos, Psi, Omega)

        # the model works only if the shader is at the origin of the coordinates system
        diff_pos = europa_pos - io_pos

        #log the positions
        diffpos_iocentric_jovian_eq[t, :] = diff_pos
        io_position_array[t, :] = io_pos
        europa_position_array[t, :] = europa_pos

        jupiter_position = Jupiter.geometric_heliocentric_position(epoch)
        earth_position = Earth.geometric_heliocentric_position(epoch)

        # convert everything to cartesian and scale to km
        x_jupiter = jupiter_position[2] * np.cos(jupiter_position[1].rad()) * np.cos(jupiter_position[0].rad())
        y_jupiter = jupiter_position[2] * np.cos(jupiter_position[1].rad()) * np.sin(jupiter_position[0].rad())
        z_jupiter = jupiter_position[2] * np.sin(jupiter_position[1].rad())

        jupiter_pos_rectangular_heliocentric = np.array([x_jupiter, y_jupiter, z_jupiter])
        jupiter_pos_rectangular_heliocentric = jupiter_pos_rectangular_heliocentric * au_to_km

        x_earth = earth_position[2] * np.cos(earth_position[1].rad()) * np.cos(earth_position[0].rad())
        y_earth = earth_position[2] * np.cos(earth_position[1].rad()) * np.sin(earth_position[0].rad())
        z_earth = earth_position[2] * np.sin(earth_position[1].rad())

        earth_pos_rectangular_heliocentric = np.array([x_earth, y_earth, z_earth])
        earth_pos_rectangular_heliocentric = earth_pos_rectangular_heliocentric * au_to_km

        # transform all coordinates to the jovian equatorial frame
        sun_pos_eq_jvc = - jupiter_pos_rectangular_heliocentric#coordTrans.HelioEcJovEq(jupiter_pos_rectangular_heliocentric, -earth_pos_rectangular_heliocentric)
        earth_pos_eq_jvc = coordTrans.HelioEcJovEq(jupiter_pos_rectangular_heliocentric,
                                                   earth_pos_rectangular_heliocentric)

        # translate to io as new coordinates origin
        sun_pos_eq_ioc = sun_pos_eq_jvc - io_pos
        earth_pos_eq_ioc = earth_pos_eq_jvc - io_pos

        # compute the Io shadow on Europa and the quantify of light reflected toward Earth
        io_shadow_array[:, :, t] = europa_shaded.ComputeMutualShadow(sun_pos_eq_ioc, diff_pos, io_shader)
        reflected_light_array[:, :, t] = europa_shaded.ComputeReflectedLight(diff_pos, earth_pos_eq_ioc, 1.0)

    sun_illumination_array_unfiltered = np.zeros([n_points, n_points, n_timesteps])

    '''
    epoch.set(dt)
    #epoch = epoch + 5200 * s_jd

    # the following returns the shadow caused by the moon itself
    for t in range(0, n_timesteps):

        # set epochs and get coordinates
        epoch = epoch + s_jd*60 # advance in 1s steps

        satellite_positions = np.array(JupiterMoons.rectangular_positions_jovian_equatorial(epoch))
        io_pos = satellite_positions[0, :] * jupiter_eq_radius_to_km
        io_pos = coordTrans.BesselianEq(io_pos, Psi, Omega)
        europa_pos = satellite_positions[1, :] * jupiter_eq_radius_to_km
        europa_pos = coordTrans.BesselianEq(europa_pos, Psi, Omega)

        jupiter_position = Jupiter.geometric_heliocentric_position(epoch)
        earth_position = Earth.geometric_heliocentric_position(epoch)

        # convert everything to cartesian
        x_jupiter = jupiter_position[2] * np.cos(jupiter_position[1].rad()) * np.cos(jupiter_position[0].rad())
        y_jupiter = jupiter_position[2] * np.cos(jupiter_position[1].rad()) * np.sin(jupiter_position[0].rad())
        z_jupiter = jupiter_position[2] * np.sin(jupiter_position[1].rad())

        jupiter_pos_rectangular_heliocentric = np.array([x_jupiter, y_jupiter, z_jupiter])
        jupiter_pos_rectangular_heliocentric = jupiter_pos_rectangular_heliocentric * au_to_km

        x_earth = earth_position[2] * np.cos(earth_position[1].rad()) * np.cos(earth_position[0].rad())
        y_earth = earth_position[2] * np.cos(earth_position[1].rad()) * np.sin(earth_position[0].rad())
        z_earth = earth_position[2] * np.sin(earth_position[1].rad())

        earth_pos_rectangular_heliocentric = np.array([x_earth, y_earth, z_earth])
        earth_pos_rectangular_heliocentric = earth_pos_rectangular_heliocentric * au_to_km

        #convert to jovian equatorial
        sun_pos_eq_jvc = coordTrans.HelioEcJovEq(jupiter_pos_rectangular_heliocentric,
                                                 -earth_pos_rectangular_heliocentric)
        earth_pos_eq_jvc = coordTrans.HelioEcJovEq(jupiter_pos_rectangular_heliocentric,
                                                   earth_pos_rectangular_heliocentric)

        # actually I should use the self shadowing method, nevertheless it is better to use a simpler approach to have
        # a lower computational toll and is a good approximation if we consider io being pretty small when compared to
        # the sun
        sun_illumination_array_unfiltered[:, :, t] = europa_shaded.ComputeReflectedLight(europa_pos, sun_pos_eq_jvc, 1)

    '''

    '''
    sun_illumination_array_filtered = io_shaded.FilterInstabilityOut(sun_illumination_array_unfiltered)
    sun_illumination_array_filtered = sun_illumination_array_filtered[:, :, filter_depth:(n_timesteps+filter_depth)]

    fig_2d_1 = px.imshow(sun_illumination_array_unfiltered[:, :, :], animation_frame=2)
    fig_2d_1.update_layout(title="Unfiltered sun Illumination array")
    fig_2d_1.show()

    fig_2d_1 = px.imshow(sun_illumination_array_filtered[:, :, :], animation_frame=2)
    fig_2d_1.update_layout(title="Filtered sun Illumination array")
    fig_2d_1.show()
    
    illumination_array = np.clip((io_shadow_array[:, :, 0:(np.shape(io_shadow_array)[2]-14)] + (1-sun_illumination_array_filtered)), 0, 1)
    '''

    # add the array of the shadows caused by the sun and by jupiter itself

    sun_illumination_array_filtered = sun_illumination_array_unfiltered
    sun_illumination_array_filtered[sun_illumination_array_filtered > 0] = 1

    fig_sun_filtered = px.imshow(sun_illumination_array_filtered[:, :, :], animation_frame=2)
    fig_sun_filtered.update_layout(title="Combination of sun and jovian shadows on Io")
    fig_sun_filtered.show()


    illumination_array = np.clip((io_shadow_array[:, :, 0:(np.shape(io_shadow_array)[2])] + (1-sun_illumination_array_filtered)), 0, 1)
    illumination_array = illumination_array * reflected_light_array[:, :, 0:(np.shape(io_shadow_array)[2])]


    fig_2d = px.imshow(illumination_array[:, :, :], animation_frame=2)
    fig_2d.update_layout(title="Combination of sun and jovian shadows on Io")
    fig_2d.show()

    #fig_2d_2 = px.imshow(io_shadow_array[:, :, :], animation_frame=2)
    #fig_2d_2.update_layout(title="Jovian shadow on Io")
    #fig_2d_2.show()

    sum_array = np.apply_over_axes(np.sum, illumination_array, [0, 1])
    sum_array = np.reshape(sum_array, [n_timesteps, ])

    fig_2_1 = px.line(y=sum_array, x=np.linspace(1, n_timesteps, n_timesteps))
    fig_2_1.update_layout(title= "Shadow cast by jupiter on Io and shadow caused by the sun")
    fig_2_1.show()

    '''
    sum_array = np.apply_over_axes(np.sum, io_shadow_array, [0,1]) #jovian shadow array
    sum_array = np.reshape(sum_array, [1500,])

    print(np.shape(sum_array))

    fig_2 = px.line(y=sum_array, x=np.linspace(1, n_timesteps, n_timesteps))
    fig_2.update_layout(title= "Shadow cast by jupiter on Io")
    fig_2.show()

    sum_array_sun = np.apply_over_axes(np.sum, sun_illumination_array_filtered* reflected_light_array[:, :, 0:(np.shape(io_shadow_array)[2]-14)], [0, 1])  # jovian shadow array
    sum_array_sun = np.reshape(sum_array_sun, [n_timesteps-14, ])

    print(np.shape(sum_array))

    final_light_array = np.apply_over_axes(np.sum, illumination_array, [0, 1])  # jovian shadow array
    final_light_array = np.reshape(final_light_array, [n_timesteps - 14, ])

    print(np.shape(sum_array))

    fig_2 = px.line(y=final_light_array, x=np.linspace(1, n_timesteps - 14, n_timesteps - 14))
    fig_2.update_layout(title="Light reflected thowars earth")
    fig_2.show()
    '''

    # scatter plot the coordinates of io and europa

    fig = go.Figure(data=[go.Scatter3d(x=io_position_array[:, 0], y=io_position_array[:, 1], z=io_position_array[:, 2])])
    fig.add_scatter3d(x = europa_position_array[:, 0], y = europa_position_array[:, 1], z = europa_position_array[:, 2])
    fig.show()

    min_distance = np.min(np.sqrt((io_position_array[:, 0]-europa_position_array[:, 0])**2 +
                                    (io_position_array[:, 2]-europa_position_array[:, 1])**2 +
                                    (io_position_array[:, 2]-europa_position_array[:, 2])**2))

    print("Mean Distance : " + str(min_distance))

    # create fictive sphere

    dt = datetime.datetime(2021, 6, 7, 1, 37, 0)

    epoch.set(dt)

    n_fictive_points = 30
    r = 10**4 # min_distance - 10**4
    phi_arr = np.linspace(0, 2 * np.pi, n_fictive_points)
    theta_arr = np.linspace(0, +np.pi, n_fictive_points)

    phi_grid, theta_grid = np.meshgrid(phi_arr, theta_arr)

    x_arr = (min_distance - 10**4) * np.cos(phi_grid) * np.sin(theta_grid)
    y_arr = (min_distance - 10**4) * np.sin(phi_grid) * np.sin(theta_grid)
    z_arr = (min_distance - 10**4) * np.cos(theta_grid)

    shadow_sphere_grid = np.zeros([n_points, n_points])

    jupiter_position = Jupiter.geometric_heliocentric_position(epoch)
    earth_position = Earth.geometric_heliocentric_position(epoch)

    # convert everything to cartesian
    x_jupiter = jupiter_position[2] * np.cos(jupiter_position[1].rad()) * np.cos(jupiter_position[0].rad())
    y_jupiter = jupiter_position[2] * np.cos(jupiter_position[1].rad()) * np.sin(jupiter_position[0].rad())
    z_jupiter = jupiter_position[2] * np.sin(jupiter_position[1].rad())

    jupiter_pos_rectangular_heliocentric = np.array([x_jupiter, y_jupiter, z_jupiter])
    jupiter_pos_rectangular_heliocentric = jupiter_pos_rectangular_heliocentric * au_to_km

    x_earth = earth_position[2] * np.cos(earth_position[1].rad()) * np.cos(earth_position[0].rad())
    y_earth = earth_position[2] * np.cos(earth_position[1].rad()) * np.sin(earth_position[0].rad())
    z_earth = earth_position[2] * np.sin(earth_position[1].rad())

    earth_pos_rectangular_heliocentric = np.array([x_earth, y_earth, z_earth])
    earth_pos_rectangular_heliocentric = earth_pos_rectangular_heliocentric * au_to_km

    sun_pos_eq_jvc = coordTrans.HelioEcJovEq(jupiter_pos_rectangular_heliocentric,
                                             -earth_pos_rectangular_heliocentric)

    diff_pos = np.array([0,0,0])

    fictive_sphere_shader = ShaderClass(r, r, delta, defaultSolver)
    fictive_sphere_shaded = ShadedClass(r, n_fictive_points, 0, delta, io_shader, fictive_sphere_shader)
    shadow_sphere_grid = fictive_sphere_shaded.ComputeMutualShadow(sun_pos_eq_jvc, diff_pos, io_shader)

    print("Io position : " + str(europa_position_array[int(len(europa_position_array) -1), :]))

    fig_shadow_sphere = go.Figure(go.Surface(x=x_arr, y=y_arr, z=z_arr, surfacecolor=shadow_sphere_grid.T))
    fig_shadow_sphere.add_scatter3d(x=diffpos_iocentric_jovian_eq[:, 0],
                                    y=diffpos_iocentric_jovian_eq[:, 1],
                                    z=diffpos_iocentric_jovian_eq[:, 2])
    fig_shadow_sphere.show()



test()
