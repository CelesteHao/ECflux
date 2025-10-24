"""
contact: Hao Zhang (hao.zhang1@wsu.edu)
date: 19/10/2025
    To calculate fluxes over Hn-1 site.

update: 20/10/2025
    test update.

update: 24/10/2025
         to modify a problem in /2019 data, just select data from begin in 2019 to 03/25/2019
"""
# calculate 30min flux from 10hz nc file for INL site, change input and output directory in the main function, output calculated flux data saved in excel
import netCDF4 as nc
import pandas as pd
import numpy as np
import glob
import os
import traceback
import warnings

warnings.filterwarnings('ignore')

# ------------------------------------
## Defined functions
# ------------------------------------
def despike_std(data, threshold):
    time_block_start = np.arange(0, 24 * 60 * 60 * 10, 30 * 60 * 10)
    time_block_end = time_block_start + 30 * 60 * 10
    for row_num in range(0, len(time_block_start)):
        data_block = data[time_block_start[row_num]:time_block_end[row_num]]
        std_num = np.nanstd(data_block)
        limit_up = np.nanmean(data_block) + threshold * std_num
        limit_dn = np.nanmean(data_block) - threshold * std_num
        data_block[data_block > limit_up] = np.nan
        data_block[data_block < limit_dn] = np.nan
    return data


def despike_MADstd(data, threshold):
    time_block_start = np.arange(0, 24 * 60 * 60 * 10, 30 * 60 * 10)
    time_block_end = time_block_start + 30 * 60 * 10
    for row_num in range(0, len(time_block_start)):
        data_block = data[time_block_start[row_num]:time_block_end[row_num]]
        MAD = np.nanmedian(abs(data_block - np.nanmedian(data_block)))  # use block median
        MAD_derived_std = 1.4826 * MAD  # assume normal distribution
        limit_up = np.nanmean(data_block) + threshold * MAD_derived_std
        limit_dn = np.nanmean(data_block) - threshold * MAD_derived_std
        data_block[data_block > limit_up] = np.nan
        data_block[data_block < limit_dn] = np.nan
    return data


def double_rotation(ux, uy, uz):
    # double rotation
    # ref: Wilczak (2001): sonic anemometer tilt correction algorithms
    # Input: ux,uy,uz: 30min 10hz despiked wind speed measuremnets (m/s)
    # Output: ux_DR/uy_DR/uz_DR: rotated wind speed (m/s)
    theta = np.arctan2(np.nanmean(uy), np.nanmean(ux))
    ux1 = ux * np.cos(theta) + uy * np.sin(theta)
    uy1 = -ux * np.sin(theta) + uy * np.cos(theta)
    uz1 = uz
    phi = np.arctan2(np.nanmean(uz1), np.nanmean(ux1))
    ux_DR = ux1 * np.cos(phi) + uz1 * np.sin(phi)
    uy_DR = uy1
    uz_DR = -ux1 * np.sin(phi) + uz1 * np.cos(phi)
    return ux_DR, uy_DR, uz_DR


def rho_a_correct(T, rho_v, P):  # T: temp K, rho_v water vapor density kg/m^3, P air pressure kPa
    e = rho_v / 18.016 * 8.3145 * T  # unit: kPa
    Tv = T / (1 - e / P * (1 - 0.622))  # virtural temp, ref: atmospheric science an introductory survey, P66
    rho_a = P / (287 * Tv) * 1000  # unit: kg/m^3
    rho_d = rho_a - rho_v
    return rho_a, rho_d


def corerct_T_SND(Ts_30min, qc, rho_d):
    # sonic temperature correction,
    # ref:  Schotanus (1983): Temperature measurement with a sonic anemometer
    # and its application to heat and moisture fluxes
    # Input: Ts:10hz sonic temperature measurements (K)
    # qc: corrected or uncorrected water vapor density (kg/m^3)
    # rho_d: dry air density (kg/m^3)
    # Outputs: Tc: corrected air temperature (K)
    Ts_mean = np.nanmean(Ts_30min)
    Ts_prime = Ts_30min - Ts_mean
    rho_a = rho_d + qc  # unit kg/m^3
    Tc_mean = Ts_mean / (1 + 0.51 * np.nanmean(qc / rho_a))
    Tc_prime = Ts_prime - 0.51 * (qc / rho_a - np.nanmean(qc / rho_a)) * Tc_mean
    Tc = Tc_mean + Tc_prime  # unit: K
    return Tc


def correct_q_WPL(Tc, q_30min, rho_d):
    # WPL correction for water vapor density,
    # ref: Webb et al., (1980): Correction of flux measurements for density
    # effects due to heat and water vapour transfer
    # Input: Tc: corrected air temperature (K)
    # q_30min: 30min 10hz measured water vapor density (kg/m^3)
    # rho_d: dry air density (kg/m^3)
    # Outputs: qc: corrected water vapor density (kg/m^3)
    Tc_mean = np.nanmean(Tc)
    Tc_prime = Tc - Tc_mean
    q_mean = np.nanmean(q_30min)
    q_prime = q_30min - q_mean
    qc_prime = q_prime + 1.61 * (q_mean / np.nanmean(rho_d)) * q_prime + (1 + 1.61 * (q_mean / np.nanmean(rho_d))) * (
                q_mean / Tc_mean) * Tc_prime
    qc = q_mean + qc_prime  # unit kg/m^3
    return qc


def correct_co2_WPL(Tc, qc, co2_30min, rho_d):
    # WPL correction for CO2
    # ref: Webb et al., (1980): Correction of flux measurements for density
    # effects due to heat and water vapour transfer
    # Input: Tc: corrected air temperature (K)
    # qc: corrected water vapor density (kg/m^3)
    # co2_30min: 30min 10hz measured CO2 density (kg/m^3)
    # rho_d: dry air density (kg/m^3)
    # Output: co2c: corrected CO2 density (kg/m^3)
    Tc_mean = np.nanmean(Tc)
    Tc_prime = Tc - Tc_mean
    co2_mean = np.nanmean(co2_30min)
    co2_prime = co2_30min - co2_mean
    q_mean = np.nanmean(qc)
    q_prime = qc - q_mean
    co2c_prime = co2_prime + 1.61 * (co2_mean / np.nanmean(rho_d)) * q_prime + (
                1 + 1.61 * (q_mean / np.nanmean(rho_d))) * (co2_mean / Tc_mean) * Tc_prime
    co2c = co2_mean + co2c_prime  # unit kg/m^3
    return co2c


# low frequency gas correction, using low frequency measurmenet as true value
def correct_T_SND_lowfreq(Ts_30min, qc, rho_d, q_mean_true):
    # sonic temperature correction,
    # ref:  Schotanus (1983): Temperature measurement with a sonic anemometer
    # and its application to heat and moisture fluxes
    # Input: Ts:10hz sonic temperature measurements (K)
    # qc: corrected or uncorrected water vapor density (kg/m^3)
    # rho_d: dry air density (kg/m^3)
    # q_mean_true: mean of low frequency water vapor density (kg/m^3)
    # Outputs: Tc: corrected air temperature (K)
    Ts_mean = np.nanmean(Ts_30min)
    Ts_prime = Ts_30min - Ts_mean
    rho_a = rho_d + q_mean_true  # unit kg/m^3
    Tc_mean = Ts_mean / (1 + 0.51 * q_mean_true / np.nanmean(rho_a))
    Tc_prime = Ts_prime - 0.51 * (qc / rho_a - np.nanmean(qc / rho_a)) * Tc_mean
    Tc = Tc_mean + Tc_prime
    return Tc


def correct_q_WPL_lowfreq(Tc, q_30min, rho_d_true, q_mean_true, T_mean_true):
    # WPL correction for water vapor density,
    # ref: Webb et al., (1980): Correction of flux measurements for density
    # effects due to heat and water vapour transfer
    # Input: Tc: corrected air temperature (K)
    # q_30min: 30min 10hz measured water vapor density (kg/m^3)
    # rho_d: dry air density (kg/m^3)
    # Outputs: qc: corrected water vapor density (kg/m^3)
    Tc_mean = np.nanmean(Tc)
    Tc_prime = Tc - Tc_mean
    q_mean = np.nanmean(q_30min)
    q_prime = q_30min - q_mean
    qc_prime = q_prime + 1.61 * (q_mean_true / np.nanmean(rho_d_true)) * q_prime + (
                1 + 1.61 * (q_mean_true / np.nanmean(rho_d_true))) * (q_mean_true / T_mean_true) * Tc_prime
    qc = q_mean_true + qc_prime  # unit kg/m^3
    return qc


def correct_co2_WPL_lowfreq(Tc, qc, co2_30min, rho_d_true, q_mean_true, T_mean_true):
    # WPL correction for CO2
    # ref: Webb et al., (1980): Correction of flux measurements for density
    # effects due to heat and water vapour transfer
    # Input: Tc: corrected air temperature (K)
    # qc: corrected water vapor density (kg/m^3)
    # co2_30min: 30min 10hz measured CO2 density (kg/m^3)
    # rho_d_true: mean dry air density (kg/m^3)
    # q_mean_true: mean low frequency water vapor measurement (kg/m^3)
    # T_mean_true: mean low frequency temperature measurement (K)
    # Output: co2c: corrected CO2 density (kg/m^3)
    Tc_mean = np.nanmean(Tc)
    Tc_prime = Tc - Tc_mean
    co2_mean = np.nanmean(co2_30min)
    co2_prime = co2_30min - co2_mean
    q_mean = np.nanmean(qc)
    q_prime = qc - q_mean
    # co2c_prime = co2_prime + 1.61*(co2_mean/np.nanmean(rho_d))*q_prime+(1+1.61*(q_mean/np.nanmean(rho_d)))*(co2_mean/Tc_mean)*Tc_prime;
    co2c_prime = co2_prime + 1.61 * (co2_mean / np.nanmean(rho_d_true)) * q_prime + (
                1 + 1.61 * (q_mean_true / np.nanmean(rho_d_true))) * (co2_mean / T_mean_true) * Tc_prime
    co2c = co2_mean + co2c_prime  # unit kg/m^3
    return co2c


def cal_flux(uz_30min_doublerotate, Ts_correct, q_correct, co2_correct, rho_air):
    # calculate corrected fluxes
    # Input: uz_30min_doublerotate: rotated w (m/s)
    # Ts_correct: corrected air temperature (K)
    # q_correct: corrected water vapor density (kg/m^3)
    # co2_correct: corrected CO2 density (kg/m^3)
    # rho_air: moist air density (kg/m^3)
    # Output: SHc: corrected sensible heat flux (W/m^2)
    # Ec: corrected latent heat flux (W/m^2)
    # Fcc: corrected co2 flux (umol/m^2/s)
    uz_mean = np.nanmean(uz_30min_doublerotate)
    uz_prime = uz_30min_doublerotate - uz_mean
    qc_mean = np.nanmean(q_correct)  # unit: kg/m^3
    qc_prime = q_correct - qc_mean
    rho_air_mean = np.nanmean(rho_air)
    Tc_mean = np.nanmean(Ts_correct)
    Tc_prime = Ts_correct - Tc_mean
    CPd = 1004.67  # J/(kg K)
    CP = CPd * (1 + 0.84 * qc_mean / rho_air_mean)  # q_mean (kg/m^3); rho_a (kg/m^3).

    wTc = np.nanmean(uz_prime * Tc_prime)
    SHc = wTc * CP * rho_air_mean  # corrected sensible heat flux

    LV = (2.501 - 0.00237 * (Tc_mean - 273.15)) * 1.0E+06  # unit J/kg, Tc_mean (K).
    wqc = np.nanmean(uz_prime * qc_prime)
    Ec = wqc * LV  # corrected latent heat flux, unit: W/m^2
    # eff_co2 = 22.72;% Convert the unit of CO2 flux from mg/(m^2 s) into vmol/(m^2 s)
    co2c_mean = np.nanmean(co2_correct)
    co2c_prime = co2_correct - co2c_mean
    wcc = np.nanmean(uz_prime * co2c_prime)
    Fcc = wcc * 1000 / 44 * 10 ** 6  # corrected co2 flux,convert unit to umol/m^2 s
    return SHc, Ec, Fcc


def cal_winddir(ux_mean, uy_mean, angle):
    # calculate 30min wind direction
    # Input: ux_mean/uy_mean: mean wind speed without rotation
    # angle: anemometer direction (true north, degree)
    # Output: direction: wind direction (degree, 0 is north)
    direction = 360 + angle - np.arctan2(uy_mean, ux_mean) * 180 / np.pi  # +360 to prevent negative angle input
    while direction > 360:
        direction = direction - 360
    # speed = (ux_mean.^2+uy_mean.^2).^(1/2);
    return direction


def cal_meanprime(data):
    data_mean = np.nanmean(data)
    data_prime = data - data_mean
    return data_mean, data_prime


def cal_flux_from1daync(ncfile, dz, date_str):

    # create dataframe for saving output flux result
    # date_1day_start = ncfile.variables['time'][0]
    # date_1day_start = pd.Timestamp(date_1day_start.item(), unit='s')
    # flux_output_time = pd.date_range(start=pd.Timestamp(date_1day_start) - pd.Timedelta(minutes=0),
    #                                  end=pd.Timestamp(date_1day_start) + pd.Timedelta(days=1), freq='30min',
    #                                  inclusive='right')

    flux_output_time = pd.date_range(start=date_str, periods=48, freq='30min')

    column_name = ['SH', 'LE', 'Fc', 'SHc', 'LEc', 'Fcc', 'SHc_low', 'LEc_low', 'Fcc_low', 'SHc_FW', 'LEc_FW', 'Fcc_FW',
                   'Bo', 'TKE', 'tau', 'u_star', 'T_star', 'q_star', 'zL', 'L', 'WS_ave', 'Winddirc', 'rho_a_mean',
                   'ux_DR_mean', 'uy_DR_mean', 'uz_DR_mean', 'Tc_mean', 'Ts_mean', 'qc_mean', 'q_mean', 'co2c_mean',
                   'co2_mean', 'Press_mean', 'FW_mean', 'Ta_mean', 'q_low_mean',
                   'ux_DR_var', 'uy_DR_var', 'uz_DR_var', 'Tc_var', 'Ts_var', 'qc_var', 'q_var', 'co2c_var', 'co2_var',
                   'Press_var', 'FW_var',
                   'cor_wu', 'cor_wv', 'cor_wTc', 'cor_wqc', 'cor_wco2c', 'cor_wFW',
                   'wwu', 'wwv', 'wwTc', 'wwqc', 'wwco2c', 'wwFW',
                   'SHc_p', 'LEc_p', 'SHc_v', 'LEc_v', 'SHc_vp', 'LEc_vp',
                   'Tc_p_mean', 'Tc_v_mean', 'Tc_vp_mean', 'Tc_p_var', 'Tc_v_var', 'Tc_vp_var'
                   ]
    units = ['W/m^2', 'W/m^2', 'umol/(m^2 s)', 'W/m^2', 'W/m^2', 'umol/(m^2 s)', 'W/m^2', 'W/m^2', 'umol/(m^2 s)',
             'W/m^2', 'W/m^2', 'umol/(m^2 s)',
             '#', 'm^2/s^2', 'kg/(m s^2)', 'm/s', 'K', 'g/m^3', '#', 'm', 'm/s', 'degree', 'kg/m^3',
             'm/s', 'm/s', 'm/s', 'K', 'K', 'g/m^3', 'g/m^3', 'mg/m^3', 'mg/m^3', 'kPa', 'K', 'K', 'g/m^3',
             '(m/s)^2', '(m/s)^2', '(m/s)^2', 'K^2', 'k^2', '(g/m^3)^2', '(g/m^3)^2', '(mg/m^3)^2', '(mg/m^3)^2',
             'kPa^2', 'K^2',
             '#', '#', '#', '#', '#', '#',
             '(m/s)^3', '(m/s)^3', '(m/s)^2*K', '(m/s)^2*(g/m^3)', '(m/s)^2*(mg/m^3)', '(m/s)^2*K',
             'W/m^2', 'W/m^2', 'W/m^2', 'W/m^2', 'W/m^2', 'W/m^2',
             'K', 'K', 'K', 'K^2', 'K^2', 'K^2'
             ]
    # df_flux = pd.DataFrame(np.nan, dtype='object', index=flux_output_time,
    #                        columns=[col + layer_name for col in column_name])
    # header = pd.MultiIndex.from_tuples(zip(df_flux.keys(), units))
    # df_flux.columns = header

    df_flux = pd.DataFrame(np.nan, index=flux_output_time, columns=column_name, dtype='float64')
    header = pd.MultiIndex.from_tuples(list(zip(column_name, units)), names=['var', 'unit'])
    df_flux.columns = header
    # read data
    data_length = len(ncfile.variables['time'][:])
    z1 = 3 #set for Hn1
    angle = 297 #set for Hn1

    ux_all = ncfile.variables['Ux'][:]  # m/s
    ux_all[ux_all > 50] = np.nan
    ux_all[ux_all < -50] = np.nan

    uy_all = ncfile.variables['Uy'][:]  # m/s
    uy_all[uy_all > 50] = np.nan
    uy_all[uy_all < -50] = np.nan

    uz_all = ncfile.variables['Uz'][:]  # m/s
    uz_all[uz_all > 10] = np.nan
    uz_all[uz_all < -10] = np.nan

    Ts_all = ncfile.variables['Ts'][:]  # C
    Ts_all[Ts_all > 50] = np.nan
    Ts_all[Ts_all < -20] = np.nan

    q_all = ncfile.variables['h2o'][:]  # g/m^3
    q_all[q_all > 30] = np.nan
    q_all[q_all < 0] = np.nan

    co2_all = ncfile.variables['co2'][:]  # mg/m^3
    co2_all[co2_all > 900] = np.nan
    co2_all[co2_all < 500] = np.nan

    Press_all = ncfile.variables['Press'][:]  # kPa
    Press_all[Press_all > 110] = np.nan
    Press_all[Press_all < 80] = np.nan

    # no Fw data, when it has ,change the line
    FW_all = ncfile.variables['Tair'][:]  # C
    FW_all[FW_all > 50] = np.nan
    FW_all[FW_all < -20] = np.nan

    Ta_all = ncfile.variables['Tair'][:]  # C
    Ta_all[Ta_all > 50] = np.nan
    Ta_all[Ta_all < -20] = np.nan

    RH_all = ncfile.variables['RH'][:]*100  # %
    RH_all[RH_all > 100] = np.nan
    RH_all[RH_all < 0] = np.nan

    diag_sonic = ncfile.variables['diag_sonic'][:]

    sig_irga = ncfile.variables['agc'][:]# this is the Hn1 special name
    # substitute masked value as np.nan
    ux_all.fill_value = np.nan
    uy_all.fill_value = np.nan
    uz_all.fill_value = np.nan
    q_all.fill_value = np.nan
    co2_all.fill_value = np.nan
    Press_all.fill_value = np.nan
    FW_all.fill_value = np.nan
    Ta_all.fill_value = np.nan
    RH_all.fill_value = np.nan

    ux_all = ux_all.filled()
    uy_all = uy_all.filled()
    uz_all = uz_all.filled()
    Ts_all = Ts_all.filled()
    q_all = q_all.filled()
    co2_all = co2_all.filled()
    Press_all = Press_all.filled()
    FW_all = FW_all.filled()
    Ta_all = Ta_all.filled()
    RH_all = RH_all.filled()

    # despike
    # first use one day std to filter bad signal strength for sonic and IRGA
    diag_sonic_std = np.nanstd(diag_sonic)
    diag_sonic[diag_sonic > np.nanmean(diag_sonic) + 3.5 * diag_sonic_std] = np.nan
    diag_sonic[diag_sonic < np.nanmean(diag_sonic) - 3.5 * diag_sonic_std] = np.nan

    sig_irga_std = np.nanstd(sig_irga)
    sig_irga[sig_irga > np.nanmean(sig_irga) + 3.5 * sig_irga_std] = np.nan
    sig_irga[sig_irga < np.nanmean(sig_irga) - 3.5 * sig_irga_std] = np.nan

    ux_all[np.isnan(diag_sonic)] = np.nan
    uy_all[np.isnan(diag_sonic)] = np.nan
    uz_all[np.isnan(diag_sonic)] = np.nan
    Ts_all[np.isnan(diag_sonic)] = np.nan

    q_all[np.isnan(sig_irga)] = np.nan
    co2_all[np.isnan(sig_irga)] = np.nan
    Press_all[np.isnan(sig_irga)] = np.nan

    threshold = 3.5
    # despike using the MAD method
    ux_all = despike_MADstd(ux_all, threshold)
    uy_all = despike_MADstd(uy_all, threshold)
    uz_all = despike_MADstd(uz_all, threshold)
    Ts_all = despike_MADstd(Ts_all, threshold)
    q_all = despike_MADstd(q_all, threshold)
    co2_all = despike_MADstd(co2_all, threshold)
    Press_all = despike_MADstd(Press_all, threshold)
    FW_all = despike_MADstd(FW_all, threshold)
    Ta_all = despike_MADstd(Ta_all, threshold)
    RH_all = despike_MADstd(RH_all, threshold)

    # convert unit
    Ts_all = Ts_all + 273.15  # convert from C to K
    Ta_all = Ta_all + 273.15  # convert from C to K
    FW_all = FW_all + 273.15  # convert from C to K
    q_all = q_all / 1000  # convert from g/m^3 to kg/m^3
    co2_all = co2_all / 10 ** 6  # convert from mg/m^3 to kg/m^3

    # calculate flux for each 30 min
    time_block_start_array = np.arange(0, 24 * 60 * 60 * 10, 30 * 60 * 10)
    time_block_end_array = time_block_start_array + 30 * 60 * 10

    for row_num in range(0, len(time_block_start_array)):
        ux = ux_all[time_block_start_array[row_num]:time_block_end_array[row_num]]
        uy = uy_all[time_block_start_array[row_num]:time_block_end_array[row_num]]
        uz = uz_all[time_block_start_array[row_num]:time_block_end_array[row_num]]
        Ts = Ts_all[time_block_start_array[row_num]:time_block_end_array[row_num]]
        q = q_all[time_block_start_array[row_num]:time_block_end_array[row_num]]
        co2 = co2_all[time_block_start_array[row_num]:time_block_end_array[row_num]]
        Press = Press_all[time_block_start_array[row_num]:time_block_end_array[row_num]]
        FW = FW_all[time_block_start_array[row_num]:time_block_end_array[row_num]]
        Ta = Ta_all[time_block_start_array[row_num]:time_block_end_array[row_num]]
        RH = RH_all[time_block_start_array[row_num]:time_block_end_array[row_num]]
        # double rotation
        ux_DR, uy_DR, uz_DR = double_rotation(ux, uy, uz)
        # calculate flux without correction
        rho_a, rho_d = rho_a_correct(Ts, q, Press)
        SH, LE, Fc = cal_flux(uz_DR, Ts, q, co2, rho_a)

        # use only high frequency loop correction
        Tc = corerct_T_SND(Ts, q, rho_d)
        qc = correct_q_WPL(Tc, q, rho_d)
        for loop in np.arange(0, 5):
            rho_a, rho_d = rho_a_correct(Tc, qc, Press)
            Tc = corerct_T_SND(Ts, qc, rho_d)
            qc = correct_q_WPL(Tc, q, rho_d)
        rho_a, rho_d = rho_a_correct(Tc, qc, Press)
        co2c = correct_co2_WPL(Tc, qc, co2, rho_d)
        SHc, LEc, Fcc = cal_flux(uz_DR, Tc, qc, co2c, rho_a)

        # calculate flux using potential temperature
        Tc_p = Tc * (Press / 100) ** 0.286
        SHc_p, LEc_p, Fcc_p = cal_flux(uz_DR, Tc_p, qc, co2c, rho_a)
        Tc_p_mean, Tc_p_prime = cal_meanprime(Tc_p)
        Tc_p_var = np.nanmean(Tc_p_prime * Tc_p_prime)
        # calculate flux using virtual temperature
        e = qc / 18.016 * 8.3145 * Tc  # unit: kPa
        Tc_v = Tc / (1 - e / Press * (1 - 0.622))  # virtural temp, ref: atmospheric science an introductory survey, P66
        SHc_v, LEc_v, Fcc_v = cal_flux(uz_DR, Tc_v, qc, co2c, rho_a)
        Tc_v_mean, Tc_v_prime = cal_meanprime(Tc_v)
        Tc_v_var = np.nanmean(Tc_v_prime * Tc_v_prime)
        # calculate flux using potential virtual temperature
        Tc_vp = Tc_v * (Press / 100) ** 0.286
        SHc_vp, LEc_vp, Fcc_vp = cal_flux(uz_DR, Tc_vp, qc, co2c, rho_a)
        Tc_vp_mean, Tc_vp_prime = cal_meanprime(Tc_vp)
        Tc_vp_var = np.nanmean(Tc_vp_prime * Tc_vp_prime)

        # calculate other required parameters
        ux_DR_mean, ux_DR_prime = cal_meanprime(ux)
        ux_DR_var = np.nanmean(ux_DR_prime * ux_DR_prime)
        uy_DR_mean, uy_DR_prime = cal_meanprime(uy)
        uy_DR_var = np.nanmean(uy_DR_prime * uy_DR_prime)
        uz_DR_mean, uz_DR_prime = cal_meanprime(uz)
        uz_DR_var = np.nanmean(uz_DR_prime * uz_DR_prime)
        Ts_mean, Ts_prime = cal_meanprime(Ts)
        Ts_var = np.nanmean(Ts_prime * Ts_prime)
        Tc_mean, Tc_prime = cal_meanprime(Tc)
        Tc_var = np.nanmean(Tc_prime * Tc_prime)
        q_mean, q_prime = cal_meanprime(q)
        q_mean = q_mean * 1000  # convert from kg/m^3 to g/m^3
        q_prime = q_prime * 1000  # convert from kg/m^3 to g/m^3
        q_var = np.nanmean(q_prime * q_prime)
        qc_mean, qc_prime = cal_meanprime(qc)
        qc_mean = qc_mean * 1000  # convert from kg/m^3 to g/m^3
        qc_prime = qc_prime * 1000  # convert from kg/m^3 to g/m^3
        qc_var = np.nanmean(qc_prime * qc_prime)
        co2_mean, co2_prime = cal_meanprime(co2)
        co2_mean = co2_mean * 10 ** 6  # convert from kg/m^3 to mg/m^3
        co2_prime = co2_prime * 10 ** 6  # convert from kg/m^3 to mg/m^3
        co2_var = np.nanmean(co2_prime * co2_prime)
        co2c_mean, co2c_prime = cal_meanprime(co2c)
        co2c_mean = co2c_mean * 10 ** 6  # convert from kg/m^3 to mg/m^3
        co2c_prime = co2c_prime * 10 ** 6  # convert from kg/m^3 to mg/m^3
        co2c_var = np.nanmean(co2c_prime * co2c_prime)
        Press_mean, Press_prime = cal_meanprime(Press)
        Press_var = np.nanmean(Press_prime * Press_prime)

        FW_mean, FW_prime = cal_meanprime(FW)
        FW_var = np.nanmean(FW_prime * FW_prime)
        Ta_mean, Ta_prime = cal_meanprime(Ta)
        Ta_var = np.nanmean(Ta_prime * Ta_prime)
        RH_mean, RH_prime = cal_meanprime(RH)
        RH_var = np.nanmean(RH_prime * RH_prime)

        WS_ave = np.nanmean((ux_DR ** 2 + uy_DR ** 2) ** (1 / 2))
        Winddirc = cal_winddir(ux_DR_mean, uy_DR_mean, angle)
        # print([np.nanmean(Tc),np.nanmean(FW)])
        cor_wu = np.nanmean(uz_DR_prime * ux_DR_prime) / (uz_DR_var ** (1 / 2) * ux_DR_var ** (1 / 2))
        cor_wv = np.nanmean(uz_DR_prime * uy_DR_prime) / (uz_DR_var ** (1 / 2) * uy_DR_var ** (1 / 2))
        cor_wTc = np.nanmean(uz_DR_prime * Tc_prime) / (uz_DR_var ** (1 / 2) * Tc_var ** (1 / 2))
        cor_wqc = np.nanmean(uz_DR_prime * qc_prime) / (uz_DR_var ** (1 / 2) * qc_var ** (1 / 2))
        cor_wco2c = np.nanmean(uz_DR_prime * co2c_prime) / (uz_DR_var ** (1 / 2) * co2c_var ** (1 / 2))
        cor_wFW = np.nanmean(uz_DR_prime * FW_prime) / (uz_DR_var ** (1 / 2) * FW_var ** (1 / 2))
        wwu = np.nanmean(uz_DR_prime * uz_DR_prime * ux_DR_prime)
        wwv = np.nanmean(uz_DR_prime * uz_DR_prime * uy_DR_prime)
        wwTc = np.nanmean(uz_DR_prime * uz_DR_prime * Tc_prime)
        wwqc = np.nanmean(uz_DR_prime * uz_DR_prime * qc_prime)
        wwco2c = np.nanmean(uz_DR_prime * uz_DR_prime * co2c_prime)
        wwFW = np.nanmean(uz_DR_prime * uz_DR_prime * FW_prime)
        rho_a_mean = np.nanmean(rho_a)
        TKE = (ux_DR_var + uy_DR_var + uz_DR_var) / 2  # m^2/s^2
        wu = np.nanmean(uz_DR_prime * ux_DR_prime)
        wv = np.nanmean(uz_DR_prime * uy_DR_prime)
        wTc = np.nanmean(uz_DR_prime * Tc_prime)
        wqc = np.nanmean(uz_DR_prime * qc_prime)
        tau = (wu ** 2 + wv ** 2) ** (1 / 2)
        u_star = tau ** (1 / 2)
        tau = rho_a_mean * tau  # Reynolds' stress
        L = -Tc_mean * (u_star ** 3) / (0.4 * 9.8 * wTc)  # k = 0.4
        zL = (z1 - dz) / L  # stability parameter
        T_star = -wTc / u_star
        q_star = -wqc / u_star
        Bo = SHc / LEc
        try:
            # use FW
            rho_a_FW, rho_d_FW = rho_a_correct(FW, q, Press)
            qc_FW = correct_q_WPL(FW, q, rho_d_FW)
            rho_a_FW, rho_d_FW = rho_a_correct(FW, qc_FW, Press)
            co2c_FW = correct_co2_WPL(FW, qc_FW, co2, rho_d_FW)
            SHc_FW, LEc_FW, Fcc_FW = cal_flux(uz_DR, FW, qc_FW, co2c_FW, rho_a_FW)
        except:
            # no FW at this layer
            SHc_FW = np.nan
            LEc_FW = np.nan
            Fcc_FW = np.nan

        try:
            # calculate low frequency h2o density from 1 hz measurement for furthre flux calculation
            es = 6.11 * 10 ** (7.5 * (Ta - 273.15) / (273.3 + Ta - 273.15))  # satuated vapor pressure, hpa
            ea = es * RH / 100
            Ta_mean = np.nanmean(Ta)
            q_low = ea * 100 * 18 * 10 ** (-3) / 8.314 / Ta
            q_low_mean = np.nanmean(q_low)
            rho_a_low, rho_d_low = rho_a_correct(Ta, q_low, Press)
            Tc_low = correct_T_SND_lowfreq(Ts, q, rho_d_low, np.nanmean(q_low))
            qc_low = correct_q_WPL_lowfreq(Tc_low, q, rho_d_low, np.nanmean(q_low), np.nanmean(Ta))
            # loop correction
            for loop in np.arange(0, 5):
                # Tc_low_old = Tc_low.copy()
                # qc_low_old = qc_low.copy()
                Tc_low = correct_T_SND_lowfreq(Ts, qc_low, rho_d_low, np.nanmean(q_low))
                qc_low = correct_q_WPL_lowfreq(Tc_low, q, rho_d_low, np.nanmean(q_low), np.nanmean(Ta))
            co2c_low = correct_co2_WPL_lowfreq(Tc_low, qc_low, co2, rho_d_low, np.nanmean(q_low), np.nanmean(Ta))
            # calculate flux from low frequency correction
            SHc_low, LEc_low, Fcc_low = cal_flux(uz_DR, Tc_low, qc_low, co2c_low, rho_a_low)
        except:
            # no low frequency measurement at this layer
            SHc_low = np.nan
            LEc_low = np.nan
            Fcc_low = np.nan

        result = [SH, LE, Fc, SHc, LEc, Fcc, SHc_low, LEc_low, Fcc_low, SHc_FW, LEc_FW, Fcc_FW,
                  Bo, TKE, tau, u_star, T_star, q_star, zL, L, WS_ave, Winddirc, rho_a_mean,
                  ux_DR_mean, uy_DR_mean, uz_DR_mean, Tc_mean, Ts_mean, qc_mean, q_mean, co2c_mean, co2_mean,
                  Press_mean, FW_mean, Ta_mean, RH_mean,
                  ux_DR_var, uy_DR_var, uz_DR_var, Tc_var, Ts_var, qc_var, q_var, co2c_var, co2_var, Press_var, FW_var,
                  cor_wu, cor_wv, cor_wTc, cor_wqc, cor_wco2c, cor_wFW,
                  wwu, wwv, wwTc, wwqc, wwco2c, wwFW,
                  SHc_p, LEc_p, SHc_v, LEc_v, SHc_vp, LEc_vp,
                  Tc_p_mean, Tc_v_mean, Tc_vp_mean, Tc_p_var, Tc_v_var, Tc_vp_var
                  ]
        df_flux.iloc[row_num, :] = result
    return df_flux



# ---------------------------------------
# config
# ---------------------------------------
data_in = '/weka/data/lab/heping_liu/haozhang/EC_Hn1/Hn1_data/2019'
dz = 0  # zero plane displacement, unit: m
output_path = '/weka/data/lab/heping_liu/haozhang/EC_Hn1/result/'+ '2019flux' + 'Hn1' + '.xlsx'  # output path
dates = pd.date_range('2019.01.01', '2019.03.25')  # modify based on the date_in
nc_files = glob.glob(os.path.join(data_in, '*.nc'))
df_flux_all = None

for date in dates:
    print('Start Processing... ', date)
    date_formated = date.strftime('%Y.%m.%d')
    matching_files = [f for f in nc_files if date_formated in os.path.basename(f)]

    try:
        if len(matching_files) > 0:
            with nc.Dataset(matching_files[0], 'r', format='NETCDF4') as ds:
                df_day = cal_flux_from1daync(ds,dz, date)

            if df_flux_all is None:
                df_flux_all = df_day
            else:
                df_flux_all = pd.concat([df_flux_all, df_day], axis=0)
        else:
            print('No data for', date_formated)

    except Exception as e:
        print(f'Error on {date_formated}: {e}')
        traceback.print_exc()
        continue

if df_flux_all is not None and not df_flux_all.empty:
    df_flux_all.to_excel(output_path, na_rep='NaN', index=True)
    print('Saved to', output_path)
else:
    print('No valid flux data to save.')

