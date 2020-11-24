from math import sin, cos, sqrt, tan, atan, atan2, radians, degrees, pi

from pymeeus.Epoch import Epoch
from pymeeus.JupiterMoons import JupiterMoons
from pymeeus.Jupiter import Jupiter
from pymeeus.Earth import Earth

import numpy as np

class Phenomena(object):

    def round_base_cone_param(epoch, ellipsoid=False):
        """This method constructs a cone modelling jupiter's shadow

        :param epoch: Epoch that should be checked
        :type epoch: :py:class:`Epoch`
        :param ellipsoid: weather or not to distort jupiter and the cone concurrently
        :type bool

        :returns: alpha_cone_rad : aperture of the cone in radians measured from the base
        :rtype: float
        :returns cone_vertex_jupiter_radii : distance of the umbral cone's sharpest point in jupiter-radii, always behind
        jupiter if viewed from the sun
        :rtype: float
        :returns: beta_cone_rad : aperture of the penumbral cone in radians measured from jupiter
        :rtype: float
        :returns: cone_beta_vertex_jupiter_radii : distance of the penumbral cone's sharpest point in jupiter-radii, always
        before jupiter if viewed from the sun
        :rtype: float
        """

        # define radius of the sun and of jupiter
        sun_radius_au = 0.00465047
        jupiter_radius_au = 0.10045 * sun_radius_au

        if (ellipsoid == True):
            jupiter_radius_multiplicator = 1.071374
        else:
            jupiter_radius_multiplicator = 1.0

        jupiter_radius_au = jupiter_radius_au * 1.071374

        # Check if Epoch is given
        if epoch is not None:
            # Check types
            if isinstance(epoch, Epoch):
                # calcuate the position of jupiter in solar-spherical coordinates
                l, b, r = Jupiter.geometric_heliocentric_position(epoch)

                # alpha is the umbra defining angle
                alpha_cone_rad = atan(r / (sun_radius_au - jupiter_radius_au))

                # beta is the penumbra defing angle
                beta_cone_rad = atan(r / (sun_radius_au + jupiter_radius_au))

                # compute distance of the sharpest pint behind jupiter in jupiter radii
                cone_vertex_jupiter_radii = tan(alpha_cone_rad)
                cone_penumbra_basis_jupiter_radii = cone_vertex_jupiter_radii * tan(beta_cone_rad) + 1  # the +1 compensates for the missing jupiter radius
                cone_beta_vertex_jupiter_radii = (-1 * jupiter_radius_multiplicator * tan(beta_cone_rad))

                return (alpha_cone_rad, cone_vertex_jupiter_radii, beta_cone_rad, cone_beta_vertex_jupiter_radii)
            else:
                raise TypeError("Invalid input type")
                return -1
    def time_correction(epoch):
        '''
        corrects a bug in pymeeus for the differential Jupiter-Earth light travel time when computing eclipses

        :param epoch: epoch to be corrcted
        :return: t_correctd : corrected time
        '''

        time_step = 1.157401129603386e-05

        l_earth, b_earth, r_earth = Earth.geometric_heliocentric_position(epoch, tofk5=True)
        l_jupiter, b_jupiter, r_jupiter = Jupiter.geometric_heliocentric_position(epoch, tofk5=True)

        l_earth = float(l_earth) * np.pi / 180
        b_earth = float(b_earth) * np.pi / 180

        l_jupiter = float(l_jupiter) * np.pi / 180
        b_jupiter = float(b_jupiter) * np.pi / 180

        x_e = r_earth * np.cos(b_earth) * np.cos(l_earth)
        y_e = r_earth * np.cos(b_earth) * np.sin(l_earth)
        z_e = r_earth * np.sin(b_earth)

        x_j = r_jupiter * np.cos(b_jupiter) * np.cos(l_jupiter)
        y_j = r_jupiter * np.cos(b_jupiter) * np.sin(l_jupiter)
        z_j = r_jupiter * np.sin(b_jupiter)

        Delta = np.sqrt((x_j - x_e) ** 2 + (y_j - y_e) ** 2 + (z_j - z_e) ** 2)

        correction_term_au = Delta - r_jupiter

        time_correction = correction_term_au * (149597870700) / (299792458)  # (150*(10**9))/(3*(10**8))

        time_correction_jd = time_step * time_correction

        t_corrected = epoch + time_correction_jd

        print("t_corretto: \t" + str(t_corrected.get_full_date()))

        return t_corrected
    def crossing_point(epoch_start, time_span, time_step=1.157401129603386e-05):

        '''
        Checks weather any kind of event happens in the time span beginning at epoch_start and containing time_span steps of
        with time_step

        :param epoch_start: begining date for scan
        :param time_span: numeber of samples
        :param time_step: time step duration in JD
        :return: Returns a matrix of tristate values indicating weather an event took place in the given time span. The moon
        number is indicated by the column of the matrix, the row indicates a Penumbral, Umbral, Shadow-Crossing, Occultation
        and Overlay events. If the tristate matrix at pos i,j has 0 then no event was detected in the time-span, if it has a
        1 then an event begun in the time span, if it has a -1 then an event ended in the time span. A matrix containing the
        timing offsets is also returned.
        '''

        epoch = epoch_start

        penumbra_matrix = [0, 0, 0, 0]
        umbra_matrix = [0, 0, 0, 0]
        occultation_matrix = [0, 0, 0, 0]
        overlap_matrix = [0, 0, 0, 0]
        shadow_matrix = [0, 0, 0, 0]

        penumbra_matrix_new = [0, 0, 0, 0]
        umbra_matrix_new = [0, 0, 0, 0]
        occultation_matrix_new = [0, 0, 0, 0]
        overlap_matrix_new = [0, 0, 0, 0]
        shadow_matrix_new = [0, 0, 0, 0]

        time_offset_matrix = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

        for i in range(time_span):

            epoch += time_step

            position_matrix = (np.array(JupiterMoons.rectangular_positions(epoch, solar=True))).reshape(4, 3)

            position_matrix[:, 1] = position_matrix[:, 1] * 1.071347

            # calculate the form of the shadow
            alpha_rad, vertex_jr, beta_rad, beta_vertex = JupiterMoons.round_base_cone_param(epoch, True)  # True

            alpha_comp = np.pi / 2 - alpha_rad

            r = np.sqrt(position_matrix[:, 0] ** 2 + position_matrix[:, 1] ** 2 + position_matrix[:, 2] ** 2)

            intersection_point = vertex_jr - r

            r_width = intersection_point * np.tan(alpha_comp)

            r_penumbra = 1 + (1 / np.tan(beta_rad)) * r

            for i_sat in range(0, 4):

                satellite_position = position_matrix[i_sat, :]

                if (satellite_position[2] > 0):

                    # check for umbra and penumbra

                    r_std = np.sqrt(position_matrix[:, 0] ** 2 + (position_matrix[:, 1]) ** 2)

                    if (r_std[i_sat] <= r_width[i_sat]):

                        if (umbra_matrix[i_sat] == umbra_matrix_new[i_sat]):

                            if (umbra_matrix[i_sat] == 0):
                                umbra_matrix_new[i_sat] = -1

                            else:
                                umbra_matrix_new[i_sat] = 1

                        if (umbra_matrix[i_sat] != umbra_matrix_new[i_sat]):
                            print("Umbral event detected on satellite : " + str(i_sat) + " @ " + str(epoch.get_full_date()))
                            time_offset_matrix[1][i_sat] = i * time_step
                            t_centro = epoch

                    umbra_matrix[i_sat] = umbra_matrix_new[i_sat]

                    if (r_std[i_sat] <= r_penumbra[i_sat]):

                        if (penumbra_matrix[i_sat] == penumbra_matrix_new[i_sat]):

                            if (penumbra_matrix[i_sat] == 0):
                                penumbra_matrix_new[i_sat] = -1
                            else:
                                penumbra_matrix_new[i_sat] = 1

                        if (penumbra_matrix[i_sat] != penumbra_matrix_new[i_sat]):
                            print("Penumbral event detected on satellite : " + str(i_sat) + " @ " + str(
                                epoch.get_full_date()))
                            time_offset_matrix[0][i_sat] = i * time_step
                            t_centro = epoch

                    penumbra_matrix[i_sat] = penumbra_matrix_new[i_sat]

                else:

                    # check for shadow
                    r_std = np.sqrt(position_matrix[:, 0] ** 2 + (position_matrix[:, 1]) ** 2)

                    if (r_std[i_sat] <= 1):

                        if (shadow_matrix[i_sat] == shadow_matrix_new[i_sat]):

                            if (shadow_matrix[i_sat] == 0):
                                shadow_matrix_new[i_sat] = -1
                            else:
                                shadow_matrix_new[i_sat] = 1

                        if (shadow_matrix[i_sat] != shadow_matrix_new[i_sat]):
                            print("Shadow event detected on satellite : " + str(i_sat) + " @ " + str(epoch.get_full_date()))
                            time_offset_matrix[4][i_sat] = i * time_step
                            t_centro = epoch

                    shadow_matrix[i_sat] = shadow_matrix_new[i_sat]

                # now consider earth based phenomena

                position_matrix_eb = (np.array(JupiterMoons.rectangular_positions(epoch, solar=False))).reshape(4, 3)

                position_matrix_eb[:, 1] = position_matrix_eb[:, 1] * 1.071347

                satellite_position_eb = position_matrix_eb[i_sat, :]

                if (satellite_position_eb[2] > 0):
                    # check for occultation

                    r_std = np.sqrt(position_matrix_eb[:, 0] ** 2 + (position_matrix_eb[:, 1]) ** 2)

                    if (r_std[i_sat] <= 1):

                        if (occultation_matrix[i_sat] == occultation_matrix_new[i_sat]):

                            if (occultation_matrix[i_sat] == 0):
                                occultation_matrix_new[i_sat] = -1
                            else:
                                occultation_matrix_new[i_sat] = 1

                        if (occultation_matrix[i_sat] != occultation_matrix_new[i_sat]):
                            print("Occultation event detected on satellite : " + str(i_sat) + " @ " + str(
                                epoch.get_full_date()))
                            time_offset_matrix[2][i_sat] = i * time_step
                            t_centro = epoch

                    occultation_matrix[i_sat] = occultation_matrix_new[i_sat]

                else:
                    # check for overlay
                    r_std = np.sqrt(position_matrix_eb[:, 0] ** 2 + (position_matrix_eb[:, 1]) ** 2)

                    if (r_std[i_sat] <= 1):

                        if (overlap_matrix[i_sat] == overlap_matrix_new[i_sat]):

                            if (overlap_matrix[i_sat] == 0):
                                overlap_matrix_new[i_sat] = -1
                            else:
                                overlap_matrix_new[i_sat] = 1

                        if (overlap_matrix[i_sat] != overlap_matrix_new[i_sat]):
                            print(
                                "Overlap event detected on satellite : " + str(i_sat) + " @ " + str(epoch.get_full_date()))
                            time_offset_matrix[3][i_sat] = i * time_step
                            t_centro = epoch

                    overlap_matrix[i_sat] = overlap_matrix_new[i_sat]

        return [[penumbra_matrix], [umbra_matrix], [occultation_matrix], [overlap_matrix], [shadow_matrix]], time_offset_matrix
    def crossing_point_avanzo(epoch_start, time_span, time_step=1.157401129603386e-05):
        '''
        Checks weather any kind of event happens in the time span beginning at epoch_start and containing time_span steps of
        with time_step, as in crossing point, this time considering the time variation determied by the radius of the moons
        of jupiter (outermost point)
        '''

        epoch = epoch_start

        # phenomena matrix :

        penumbra_matrix = [0, 0, 0, 0]
        umbra_matrix = [0, 0, 0, 0]
        occultation_matrix = [0, 0, 0, 0]
        overlap_matrix = [0, 0, 0, 0]
        shadow_matrix = [0, 0, 0, 0]

        penumbra_matrix_new = [0, 0, 0, 0]
        umbra_matrix_new = [0, 0, 0, 0]
        occultation_matrix_new = [0, 0, 0, 0]
        overlap_matrix_new = [0, 0, 0, 0]
        shadow_matrix_new = [0, 0, 0, 0]

        time_offset_matrix = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

        moons_radii = np.array(
            [[3643.2 / (2 * 69911), 3121.6 / (2 * 69911), 5262.4 / (2 * 69911), 4820.6 / (2 * 69911)], [0, 0, 0, 0],
             [0, 0, 0, 0]])

        for i in range(time_span):

            epoch += time_step

            position_matrix = (np.array(JupiterMoons.rectangular_positions(epoch, solar=True))).reshape(4, 3)

            position_matrix += moons_radii.T

            position_matrix[:, 1] = position_matrix[:, 1] * 1.071347

            # calculate the form of the shadow
            alpha_rad, vertex_jr, beta_rad, beta_vertex = JupiterMoons.round_base_cone_param(epoch, True)  # True

            alpha_comp = np.pi / 2 - alpha_rad

            r = np.sqrt(position_matrix[:, 0] ** 2 + position_matrix[:, 1] ** 2 + position_matrix[:, 2] ** 2)

            intersection_point = vertex_jr - r

            r_width = intersection_point * np.tan(alpha_comp)

            r_penumbra = 1 + (1 / np.tan(beta_rad)) * r

            for i_sat in range(0, 4):

                satellite_position = position_matrix[i_sat, :]

                if (satellite_position[2] > 0):

                    # check for umbra and penumbra

                    r_std = np.sqrt(position_matrix[:, 0] ** 2 + (position_matrix[:, 1]) ** 2)

                    if (r_std[i_sat] <= r_width[i_sat]):

                        if (umbra_matrix[i_sat] == umbra_matrix_new[i_sat]):

                            if (umbra_matrix[i_sat] == 0):
                                umbra_matrix_new[i_sat] = -1

                            else:
                                umbra_matrix_new[i_sat] = 1

                        if (umbra_matrix[i_sat] != umbra_matrix_new[i_sat]):
                            print("Umbral event detected on satellite : " + str(i_sat) + " @ " + str(epoch.get_full_date()))
                            time_offset_matrix[1][i_sat] = i * time_step
                            t_centro = epoch

                    umbra_matrix[i_sat] = umbra_matrix_new[i_sat]

                    if (r_std[i_sat] <= r_penumbra[i_sat]):

                        if (penumbra_matrix[i_sat] == penumbra_matrix_new[i_sat]):

                            if (penumbra_matrix[i_sat] == 0):
                                penumbra_matrix_new[i_sat] = -1
                            else:
                                penumbra_matrix_new[i_sat] = 1

                        if (penumbra_matrix[i_sat] != penumbra_matrix_new[i_sat]):
                            print("Penumbral event detected on satellite : " + str(i_sat) + " @ " + str(
                                epoch.get_full_date()))
                            time_offset_matrix[0][i_sat] = i * time_step
                            t_centro = epoch

                    penumbra_matrix[i_sat] = penumbra_matrix_new[i_sat]

                else:

                    # check for shadow
                    r_std = np.sqrt(position_matrix[:, 0] ** 2 + (position_matrix[:, 1]) ** 2)

                    if (r_std[i_sat] <= 1):

                        if (shadow_matrix[i_sat] == shadow_matrix_new[i_sat]):

                            if (shadow_matrix[i_sat] == 0):
                                shadow_matrix_new[i_sat] = -1
                            else:
                                shadow_matrix_new[i_sat] = 1

                        if (shadow_matrix[i_sat] != shadow_matrix_new[i_sat]):
                            print("Shadow event detected on satellite : " + str(i_sat) + " @ " + str(epoch.get_full_date()))
                            time_offset_matrix[4][i_sat] = i * time_step
                            t_centro = epoch

                    shadow_matrix[i_sat] = shadow_matrix_new[i_sat]

                # now consider earth based phenomena

                position_matrix_eb = (np.array(JupiterMoons.rectangular_positions(epoch, solar=False))).reshape(4, 3)

                position_matrix_eb += moons_radii.T

                position_matrix_eb[:, 1] = position_matrix_eb[:, 1] * 1.071347

                satellite_position_eb = position_matrix_eb[i_sat, :]

                if (satellite_position_eb[2] > 0):
                    # check for occultation

                    r_std = np.sqrt(position_matrix_eb[:, 0] ** 2 + (position_matrix_eb[:, 1]) ** 2)

                    if (r_std[i_sat] <= 1):

                        if (occultation_matrix[i_sat] == occultation_matrix_new[i_sat]):

                            if (occultation_matrix[i_sat] == 0):
                                occultation_matrix_new[i_sat] = -1
                            else:
                                occultation_matrix_new[i_sat] = 1

                        if (occultation_matrix[i_sat] != occultation_matrix_new[i_sat]):
                            print("Occultation event detected on satellite : " + str(i_sat) + " @ " + str(
                                epoch.get_full_date()))
                            time_offset_matrix[2][i_sat] = i * time_step
                            t_centro = epoch

                    occultation_matrix[i_sat] = occultation_matrix_new[i_sat]

                else:
                    # check for overlay
                    r_std = np.sqrt(position_matrix_eb[:, 0] ** 2 + (position_matrix_eb[:, 1]) ** 2)

                    if (r_std[i_sat] <= 1):

                        if (overlap_matrix[i_sat] == overlap_matrix_new[i_sat]):

                            if (overlap_matrix[i_sat] == 0):
                                overlap_matrix_new[i_sat] = -1
                            else:
                                overlap_matrix_new[i_sat] = 1

                        if (overlap_matrix[i_sat] != overlap_matrix_new[i_sat]):
                            print(
                                "Overlap event detected on satellite : " + str(i_sat) + " @ " + str(epoch.get_full_date()))
                            time_offset_matrix[3][i_sat] = i * time_step
                            t_centro = epoch

                    overlap_matrix[i_sat] = overlap_matrix_new[i_sat]

        return [[penumbra_matrix], [umbra_matrix], [occultation_matrix], [overlap_matrix],[shadow_matrix]], time_offset_matrix
    def crossing_point_regresso(epoch_start, time_span, time_step=1.157401129603386e-05):

        '''
        Checks weather any kind of event happens in the time span beginning at epoch_start and containing time_span steps of
        with time_step, as in crossing point, this time considering the time variation determied by the radius of the moons
        of jupiter (innermost point)
        '''

        epoch = epoch_start

        # phenomena matrix :

        penumbra_matrix = [0, 0, 0, 0]
        umbra_matrix = [0, 0, 0, 0]
        occultation_matrix = [0, 0, 0, 0]
        overlap_matrix = [0, 0, 0, 0]
        shadow_matrix = [0, 0, 0, 0]

        penumbra_matrix_new = [0, 0, 0, 0]
        umbra_matrix_new = [0, 0, 0, 0]
        occultation_matrix_new = [0, 0, 0, 0]
        overlap_matrix_new = [0, 0, 0, 0]
        shadow_matrix_new = [0, 0, 0, 0]

        time_offset_matrix = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

        moons_radii = np.array(
            [[3643.2 / (2 * 69911), 3121.6 / (2 * 69911), 5262.4 / (2 * 69911), 4820.6 / (2 * 69911)], [0, 0, 0, 0],
             [0, 0, 0, 0]])

        for i in range(time_span):

            epoch += time_step

            position_matrix = (np.array(JupiterMoons.rectangular_positions(epoch, solar=True))).reshape(4, 3)

            position_matrix += moons_radii.T

            position_matrix[:, 1] = position_matrix[:, 1] * 1.071347

            # calculate the form of the shadow
            alpha_rad, vertex_jr, beta_rad, beta_vertex = JupiterMoons.round_base_cone_param(epoch, True)  # True

            alpha_comp = np.pi / 2 - alpha_rad

            r = np.sqrt(position_matrix[:, 0] ** 2 + position_matrix[:, 1] ** 2 + position_matrix[:, 2] ** 2)

            intersection_point = vertex_jr - r

            r_width = intersection_point * np.tan(alpha_comp)

            r_penumbra = 1 + (1 / np.tan(beta_rad)) * r

            for i_sat in range(0, 4):

                satellite_position = position_matrix[i_sat, :]

                if (satellite_position[2] > 0):

                    # check for umbra and penumbra

                    r_std = np.sqrt(position_matrix[:, 0] ** 2 + (position_matrix[:, 1]) ** 2)

                    if (r_std[i_sat] <= r_width[i_sat]):

                        if (umbra_matrix[i_sat] == umbra_matrix_new[i_sat]):

                            if (umbra_matrix[i_sat] == 0):
                                umbra_matrix_new[i_sat] = -1

                            else:
                                umbra_matrix_new[i_sat] = 1

                        if (umbra_matrix[i_sat] != umbra_matrix_new[i_sat]):
                            print("Umbral event detected on satellite : " + str(i_sat) + " @ " + str(epoch.get_full_date()))
                            time_offset_matrix[1][i_sat] = i * time_step
                            t_centro = epoch

                    umbra_matrix[i_sat] = umbra_matrix_new[i_sat]

                    if (r_std[i_sat] <= r_penumbra[i_sat]):

                        if (penumbra_matrix[i_sat] == penumbra_matrix_new[i_sat]):

                            if (penumbra_matrix[i_sat] == 0):
                                penumbra_matrix_new[i_sat] = -1
                            else:
                                penumbra_matrix_new[i_sat] = 1

                        if (penumbra_matrix[i_sat] != penumbra_matrix_new[i_sat]):
                            print("Penumbral event detected on satellite : " + str(i_sat) + " @ " + str(
                                epoch.get_full_date()))
                            time_offset_matrix[0][i_sat] = i * time_step
                            t_centro = epoch

                    penumbra_matrix[i_sat] = penumbra_matrix_new[i_sat]

                else:

                    # check for shadow
                    r_std = np.sqrt(position_matrix[:, 0] ** 2 + (position_matrix[:, 1]) ** 2)

                    if (r_std[i_sat] <= 1):

                        if (shadow_matrix[i_sat] == shadow_matrix_new[i_sat]):

                            if (shadow_matrix[i_sat] == 0):
                                shadow_matrix_new[i_sat] = -1
                            else:
                                shadow_matrix_new[i_sat] = 1

                        if (shadow_matrix[i_sat] != shadow_matrix_new[i_sat]):
                            print("Shadow event detected on satellite : " + str(i_sat) + " @ " + str(epoch.get_full_date()))
                            time_offset_matrix[4][i_sat] = i * time_step
                            t_centro = epoch

                    shadow_matrix[i_sat] = shadow_matrix_new[i_sat]

                # now consider earth based phenomena

                position_matrix_eb = (np.array(JupiterMoons.rectangular_positions(epoch, solar=False))).reshape(4, 3)

                position_matrix_eb += moons_radii.T

                position_matrix_eb[:, 1] = position_matrix_eb[:, 1] * 1.071347

                satellite_position_eb = position_matrix_eb[i_sat, :]

                if (satellite_position_eb[2] > 0):
                    # check for occultation

                    r_std = np.sqrt(position_matrix_eb[:, 0] ** 2 + (position_matrix_eb[:, 1]) ** 2)

                    if (r_std[i_sat] <= 1):

                        if (occultation_matrix[i_sat] == occultation_matrix_new[i_sat]):

                            if (occultation_matrix[i_sat] == 0):
                                occultation_matrix_new[i_sat] = -1
                            else:
                                occultation_matrix_new[i_sat] = 1

                        if (occultation_matrix[i_sat] != occultation_matrix_new[i_sat]):
                            print("Occultation event detected on satellite : " + str(i_sat) + " @ " + str(
                                epoch.get_full_date()))
                            time_offset_matrix[2][i_sat] = i * time_step
                            t_centro = epoch

                    occultation_matrix[i_sat] = occultation_matrix_new[i_sat]

                else:
                    # check for overlay
                    r_std = np.sqrt(position_matrix_eb[:, 0] ** 2 + (position_matrix_eb[:, 1]) ** 2)

                    if (r_std[i_sat] <= 1):

                        if (overlap_matrix[i_sat] == overlap_matrix_new[i_sat]):

                            if (overlap_matrix[i_sat] == 0):
                                overlap_matrix_new[i_sat] = -1
                            else:
                                overlap_matrix_new[i_sat] = 1

                        if (overlap_matrix[i_sat] != overlap_matrix_new[i_sat]):
                            print(
                                "Overlap event detected on satellite : " + str(i_sat) + " @ " + str(epoch.get_full_date()))
                            time_offset_matrix[3][i_sat] = i * time_step
                            t_centro = epoch

                    overlap_matrix[i_sat] = overlap_matrix_new[i_sat]

        return [[penumbra_matrix], [umbra_matrix], [occultation_matrix], [overlap_matrix],
                [shadow_matrix]], time_offset_matrix

#usage example
epoch = Epoch()

epoch.set(2020, 1, 2, 12, 30, 0)

base_epoch = epoch

event_matrix, time_matrix = Phenomena.crossing_point(epoch, 600)

eclipse_vector = time_matrix[1][:]
penumbra_vector = time_matrix[0][:]
shadow_vector = time_matrix[4][:]

for i in eclipse_vector :
    if(i != 0):
        i = Phenomena.time_correction(base_epoch + i)

for i in penumbra_vector :
    if (i != 0):
        i = Phenomena.time_correction(base_epoch + i)

for i in shadow_vector:
    if (i != 0):
        i = Phenomena.time_correction(base_epoch + i)