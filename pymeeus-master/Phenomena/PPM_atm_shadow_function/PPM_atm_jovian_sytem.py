import numpy as np
import math

import matplotlib.pyplot as plt

from math import sqrt, fabs, acos, acosh, isnan, atan, atan2, sin, cos, pi

import cmath
from cmath import sqrt as sqrt_complex

# solve cubic equation x^3 + a*x^2 + b*x + c
# x - array of size 3
# In case 3 real roots: => x[0], x[1], x[2], return 3
#         2 real roots: x[0], x[1],          return 2
#         1 real root : x[0], x[1] ± i*x[2], return 1

def solveP3(a, b, c):
    a2 = a**2
    q = (a2 - 3 * b) / 9
    r = (a * (2 * a2 - 9 * b) + 27 * c) / 54
    r2 = r**2
    q3 = q**3
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
        a = a/3
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

        a = a/3
        x[0] = (A + B) - a
        x[1] = -0.5 * (A + B) - a
        x[2] = 0.5 * sqrt(3.) * (A - B)
        if fabs(x[2]) < eps:
            x[2] = x[1]
            return x, 2

        return x, 1


# solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d
def solve_quartic(a, b, c, d):
    eps = 1*10**-5

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

    x3, iZeroes = solveP3(coef3[0], coef3[1], coef3[2])

    s = x3[0]

    # The essence - choosing Y with maximal absolute value.

    if iZeroes != 1:

        if fabs(x3[1]) > fabs(s):
            s = x3[1]
        if fabs(x3[2]) > fabs(s):
            s = x3[2]

    # what is the solution ? A= 2s, i.e. B==0
    # y^4 + Ay^2 + C = 0
    if fabs(B) < eps:


        D = A * A - 4 * C
        x1 = (-A + sqrt_complex(D)) * 0.5
        x2 = (-A - sqrt_complex(D)) * 0.5

        #no abs originally
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
            x[num_real] = np.real(retval[i])
            num_real = num_real + 1

    return x, num_real


def area_hyperbola(Q1, Q2, Os, Oe, Rs, a, b, in_out, x_axis):

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
    #print("cts : " + str(cts))

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

        s_trapizium = 0.5 * fabs(Q1[0] - Q2[0]) * (fabs(Q1[1]) + fabs(Q2[1]))
        S2 = (sl - ss - 2 * s_trapizium) / 2.0

    if in_out:

        area = SQ1Q2Os - TQ1Q2Os - S2

    else:  # the centre of the sun is outside the ellipse

        area_shadow = SQ1Q2Os - TQ1Q2Os + S2
        area = pi * Rs * Rs - area_shadow

    return area


def area_ellispe(Q1, Q2, Os, Oe, Rs, a, b, in_out):
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

    #print("cts : " + str(cts))

    if np.abs(cts) > 1:
        print("Math domain error \n aborting execution with exit code -1")
        return(-1)

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

        q1 = 3.14159265357 - atan(-a_e[1] / a_e[0])

    elif a_e[0] < 0 and a_e[1] < 0:

        q1 = 3.14159265357 + atan(a_e[1] / a_e[0])

    elif a_e[0] > 0 and a_e[1] < 0:

        q1 = 3.14159265357 * 2 - atan(-a_e[1] / a_e[0])

    if b_e[0] > 0 and b_e[1] > 0:  # the first quadrant

        q2 = atan(b_e[1] / b_e[0])

    elif b_e[0] < 0 and b_e[1] > 0:

        q2 = 3.14159265357 - atan(-b_e[1] / b_e[0])

    elif b_e[0] < 0 and b_e[1] < 0:

        q2 = 3.14159265357 + atan(b_e[1] / b_e[0])

    elif b_e[0] > 0 and b_e[1] < 0:

        q2 = 3.14159265357 * 2 - atan(-b_e[1] / b_e[0])

    s_a = a * b / 2.0 * (q1 - atan2((b - a) * sin(2.0 * q1), (a + b) + (b - a) * cos(2 * q1)))
    s_b = a * b / 2.0 * (q2 - atan2((b - a) * sin(2.0 * q2), (a + b) + (b - a) * cos(2 * q2)))

    SQ1Q2Oe = fabs(s_b - s_a)

    S1 = SQ1Q2Oe - TQ1Q2Oe

    if in_out:  # the centre of the sun is inside the ellipse

        area = SQ1Q2Os - TQ1Q2Os - S1

    else:  # the centre of the sun is outside the ellipse

        area_shadow = SQ1Q2Os - TQ1Q2Os + S1
        area = 3.14159265357 * Rs * Rs - area_shadow

    return area


def myperspectiveProjection(a, b, sunpos_ecef, satpos_ecef, r_solar, area_bright, dis_boundary, dis_circle):
    state = -1
    Rs = 695700.0  # Radius of the Sun in km , we assume the sun is a sphere in this model
    ab = a * b
    a2 = a * a
    b2 = b * b
    ab2 = ab * ab

    r = satpos_ecef
    rs = sunpos_ecef

    t = 0.0
    t1 = 0.0
    t2 = 0.0
    dis = 0
    s1 = 0
    s2 = 0
    ds1 = 0.0
    ds2 = 0.0
    # A: first test if the satellite is in the front of earth,
    # if it is in the front of earth, it is always full phase
    # otherwise,using the photogrammetry method
    # A= diag{1/a2,1/a2, 1/b2 }, A^{-1} = diag{a2, a2, b2}

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

    ds1 = t

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
    o = r - f * n  # the origin of the photo coordinate system (PSC)

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
    t = sqrt((K[0] + K[2]) * (K[0] + K[2]) - 4.0 * (K[0] * K[2] - K[1] * K[1]))
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

    #print("A, B, C, D, E : ")

    #print(A)
    #print(B)
    #print(C)
    #print(D)
    #print(E)

    ce = np.array([B / A, C / A, D / A, E / A], dtype=np.double)  # object

    XX, num_of_solution = solve_quartic(ce[0], ce[1], ce[2], ce[3])

    # XX = np.roots(ce)
    #print(XX)
    # num_of_solution = XX.size

    #print(num_of_solution)

    if num_of_solution == 3 or num_of_solution == 4:
        print("WARNING: shadowfactor, intersections -->" + str(num_of_solution))
        # exit(0) #-->commented because there could be more than one equal soliution

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

    #print("T : " + str(t) + "out of : " + str(t1) + " , " + str(t2))

    dis_circle = t

    t = K[0] * K[2] - K[1] * K[1]

    #print("t : " + str(t))

    if t > 0.0:  # ellipse

        #print("We should compute an elliptic area")
        in_out = False

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
        area_bright = area_ellispe(Q1, Q2, Os, _, R0, sqrt(AA), sqrt(BB), in_out)

        #print("Area Bright")
        #print(area_bright)

    elif t < 0.0:  # hyperbola

        #print("We should compute an hyperbolic area")

        in_out = False

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

        area_bright = area_hyperbola(Q1, Q2, Os, PEC_new, R0, a, b, in_out, x_axis)

    return state, r_solar, area_bright, dis_boundary, dis_circle


# the main entry of the shadow function model
def myshadowFactor(sunpos_eci, satpos_eci):
    factor = 1.0

    a = 7.1492*10**4  # km
    b = 6.6854*10**4  # km

    # for sphere test
    # a = 6371.0 #km
    # b = 6371.0 #km  6356.7523142

    hgt_atm = 10**3  # the hight of atmosphere, this parameter is very important

    # considering atmosphere
    a_atm = a + hgt_atm  # km
    b_atm = b * (a_atm / a)  # km  6356.7523142

    Area_earth = 0.0
    Area_atm = 0.0
    Area_solar = 0.0
    state1 = -1
    state2 = -1

    x1 = 0.0
    x2 = 0.0
    dis0 = 0.0
    dis1 = 0.0
    dis2 = 0.0

    r_sun = 0.0

    if hgt_atm == 0.0:

        state2, r_sun, Area_earth, dis1, dis0 = myperspectiveProjection(a, b, sunpos_eci, satpos_eci, r_sun, Area_earth,
                                                                        dis1, dis0)
        Area_solar = 3.14159265357 * r_sun * r_sun

        if state2 == 0:

            factor = 1.0

        elif state2 == -1:

            factor = 0.0

        else:

            factor = Area_earth / Area_solar



    else:  # considering the atmosphere

        state1, r_sun, Area_atm, dis1, dis0 = myperspectiveProjection(a_atm, b_atm, sunpos_eci, satpos_eci, r_sun,
                                                                      Area_atm, dis1, dis0)
        state2, r_sun, Area_earth, dis2, dis0 = myperspectiveProjection(a, b, sunpos_eci, satpos_eci, r_sun, Area_earth,
                                                                        dis2, dis0)

        mu1 = 0.0  # distance = 0
        mu2 = 1.0  # distance = 1

        a_linear = mu1
        b_linear = mu2

        # so the linear model would be y = 0.6x + 0.3
        # a_exp = mu1
        # b_exp = np.log(mu2 / mu1)

        # a_log = mu1
        # b_log = (mu2 - mu1) / np.log(2.0)

        Area_solar = 3.14159265357 * r_sun * r_sun

        #print("Area Solar")
        #print(Area_solar)

        # the thickness of atmosphere
        thickness = dis1 - dis2

        #print("Thickness" + str(thickness))

        if state1 == 0:  # full phase

            factor = 1.0

        # partly in the atmosphere
        if (state1 == 1 or state1 == 2) and state2 == 0:
            x1 = 1.0
            x2 = (thickness - (dis1 - dis0)) / thickness

            # log
            # u1 = a_log + b_log*log(x1+1.0)
            # u2 = a_log + b_log*log(x2+1.0)

            # linear
            u1 = (mu2 - mu1) * x1 + mu1
            u2 = (mu2 - mu1) * x2 + mu1

            coeff = (u1 + u2) / 2.0
            factor = (Area_atm + (Area_solar - Area_atm) * coeff) / Area_solar

            #print("factor" + str(factor))

        # totally in the atmosphere
        if state1 == -1 and state2 == 0:

            u1 = 0
            u2 = 0
            if dis2 != 0.0:

                x1 = (dis0 - dis2) / thickness
                x2 = (dis0 - dis2 + 2.0 * r_sun) / thickness


                # linear
                u1 = (mu2 - mu1) * x1 + mu1
                u2 = (mu2 - mu1) * x2 + mu1

            elif (
                    fabs(dis2) < 1.0E-9):  # no projection for the solid Earth (outside of the solid Earth), but the satellite is in the atmosphere's projection

                print("no projection for the solid Earth (outside of the solid Earth), but the satellite is in the atmosphere's projection")

                u1 = 0.5
                u2 = 0.5

            factor = (u1 + u2) / 2.0

        # partly in the earth and atmosphere
        if state1 == -1 and (state2 == 1 or state2 == 2):
            x1 = (2.0 * r_sun - (dis2 - dis0)) / thickness
            x2 = 0.0

            # log
            # u1 = a_log + b_log*log(x1+1.0)
            # u2 = a_log + b_log*log(x2+1.0)

            # linear
            u1 = (mu2 - mu1) * x1 + mu1
            u2 = (mu2 - mu1) * x2 + mu1

            coeff = (u1 + u2) / 2.0

            factor = (Area_earth * coeff) / Area_solar

        # totally in the earth's shadow, umbra
        # in umbra, the CERES data will capture the atmosphere effect
        if state2 == -1:
            factor = 0.0

        # the sun is very big, bigger than the atmosphere thickness
        if (state1 == 1 or state1 == 2) and (state2 == 1 or state2 == 2):
            #print("hello there")

            area1 = Area_atm
            area2 = Area_earth - Area_atm

            #print("areas : " + str(area1) + str(area2))

            x1 = 1.0
            x2 = 0.0

            # linear
            u1 = (mu2 - mu1) * x1 + mu1
            u2 = (mu2 - mu1) * x2 + mu1

            factor = (area1 * 1.0 + area2 * (u1 + u2) / 2.0) / Area_solar

    return factor


def test():
    # test time: UTC: 2015 1 11 18 33 34
    # shadow function should give: 0.474838

    print("Beginning test")

    satpos_x =-13205.655784525363
    satpos_y = 21522.519302073124
    satpos_z = 15446.72240793841

    satpos = np.array([satpos_x, satpos_y, satpos_z], dtype=np.double)
    sunpos = np.array([52727703.80386541, -126017147.89721917, -54630443.258015752], dtype=np.double)

    shadow_factor = myshadowFactor(sunpos, satpos)
    print("light factor : " + str(shadow_factor))

    #now test it with several points

    #for i in range(1000) :

    #    satpos_y = satpos_y - i * 0.1

    #    satpos = np.array([satpos_x, satpos_y, satpos_z], dtype=np.double)

    #    shadow_factor = myshadowFactor(sunpos, satpos)

    #    print("pos : " + str(i) + " || light factor : " + str(shadow_factor))

    return 0





test()
