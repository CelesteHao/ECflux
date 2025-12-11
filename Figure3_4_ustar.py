# This function is for figure3_ in NEE report
# show u* threshold
# just a brief show
# Author: Hao Zhang, Date: 12/10/2025

# import modules
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# setting the path
input_path = 'C:/split/output.txt'
output_path = 'C:/split/figure'

# whether draw figure
drawfigs = [True, True]
figure_name = ['ustar_distribution','ustar_season']
figure_size = (12,7)
data = pd.read_csv(input_path, sep='\t', na_values='-9999', header=0, skiprows=[1], parse_dates=['Date Time'])
data = data.set_index('Date Time')
# this is because python will consider 00:00 as the day's data not the last day, so do a shift
# and make a copy avoid influence the raw data
data1 = data.copy()
data1.index = data1.index - pd.Timedelta(minutes=30)
# test input data
# print(data.head())
# calculate the daily value
# columns1 = ['H_f', 'H_orig', 'LE_f', 'LE_orig', 'NEE_uStar_f', 'NEE_uStar_orig'] # fig1
# columns2 = ['NEE_uStar_f', 'GPP_DT_uStar','Reco_DT_uStar','GPP_uStar_f','Reco_uStar']
# daily = data1[columns1].resample('D')
# daily_mean  = daily.mean()
# daily_count = daily.count()
#
# daily1 = data1[columns2].resample('D')
# daily_sum = daily1.sum()
# daily_sum_count = daily1.count()
# daily_mean[daily_count < 24] = np.nan

night = data['PotRad_U50'] == 0
ustar = data.loc[night, 'Ustar']
print(len(ustar))
ustar05 = data['Ustar_U05_Thres'].iloc[0]
ustar16 = data['Ustar_U16_Thres'].iloc[0]
ustar25 = data['Ustar_U25_Thres'].iloc[0]
ustar50 = data['Ustar_U50_Thres'].iloc[0]
ustar75 = data['Ustar_U75_Thres'].iloc[0]
ustar84 = data['Ustar_U84_Thres'].iloc[0]
ustar95 = data['Ustar_U95_Thres'].iloc[0]

# print(ustar05,
#       ustar16,
#       ustar25,
#       ustar50,
#       ustar75,
#       ustar84,
#       ustar95)

# figure
# a figure show the u* distribution
if drawfigs[0]:
    al = 0.3
    fig1 = plt.figure(dpi=300, figsize=figure_size, clear=True)
    ax1 = fig1.add_subplot(1, 2, 1)
    ax1.hist(ustar, bins=200, density=False, alpha=al, edgecolor='black', facecolor='blue')
    ax1.axvline(ustar05, color='black', linestyle='--', linewidth=1, label='U05')
    ax1.axvline(ustar16, color='black', linestyle='--', linewidth=1, label='U16')
    ax1.axvline(ustar25, color='black', linestyle='--', linewidth=1, label='U25')
    ax1.axvline(ustar50, color='red', linestyle='-',  linewidth=3.0, label='U50')
    ax1.axvline(ustar75, color='black', linestyle='--', linewidth=1, label='U75')
    ax1.axvline(ustar84, color='black', linestyle='--', linewidth=1, label='U84')
    ax1.axvline(ustar95, color='black', linestyle='--', linewidth=1, label='U95')
    ax1.set_xlabel('uStar\n (m/s)',fontsize=16)
    # ax1.set_xlim(0, 0.4)
    ax1.set_ylabel('Number of uStars', fontsize=16)
    # ax1.legend(loc='best', fontsize=8)
    ax1.grid(True, alpha=0.3)


    ax2 = fig1.add_subplot(1, 2, 2)
    ax2.hist(ustar, bins=200, density=False, alpha=al, edgecolor='black', facecolor='blue')
    ax2.axvline(ustar05, color='red', linestyle='--', linewidth=1, label='U05')
    ax2.axvline(ustar16, color='red', linestyle='--', linewidth=1, label='U16')
    ax2.axvline(ustar25, color='red', linestyle='--', linewidth=1, label='U25')
    ax2.axvline(ustar50, color='red', linestyle='-',  linewidth=3.0, label='U50')
    ax2.axvline(ustar75, color='red', linestyle='--', linewidth=1, label='U75')
    ax2.axvline(ustar84, color='red', linestyle='--', linewidth=1, label='U84')
    ax2.axvline(ustar95, color='red', linestyle='--', linewidth=1, label='U95')
    ax2.set_xlabel('uStar\n (m/s)', fontsize=16)
    ax2.set_xlim(0, 0.2)
    # ax2.set_ylabel('Number of uStars',fontsize=14)
    ax2.legend(loc='best', fontsize=14)
    ax2.grid(True, alpha=0.3)
    plt.savefig(output_path + '/' + figure_name[0] + '.jpg', dpi=300,
                format='jpg')




print('NEE before bootstrap:', np.sum(data['NEE_uStar_f']* 1800 * 12e-6))
print('NEE u05:', np.sum(data['NEE_U05_f'] * 1800 * 12e-6))
print('NEE u16:', np.sum(data['NEE_U16_f'] * 1800 * 12e-6))
print('NEE u25:', np.sum(data['NEE_U25_f'] * 1800 * 12e-6))
print('NEE u50:', np.sum(data['NEE_U50_f'] * 1800 * 12e-6))
print('NEE u75:', np.sum(data['NEE_U75_f'] * 1800 * 12e-6))
print('NEE u84:', np.sum(data['NEE_U84_f'] * 1800 * 12e-6))
print('NEE u95:', np.sum(data['NEE_U95_f'] * 1800 * 12e-6))

# draw seasonal u* filtering percent
winter_flag = data1['season'].astype(str).str.endswith('12')
spring_flag = data1['season'].astype(str).str.endswith('03')
summer_flag = data1['season'].astype(str).str.endswith('06')
autumn_flag = data1['season'].astype(str).str.endswith('09')
night_flag = data1['PotRad_U50'] == 0
base_night = night_flag & data1['NEE'].notna()
N_night_all = base_night.sum()

# U05
mask_winter_orig = base_night & winter_flag
mask_winter_kept_u05 = mask_winter_orig & data1['NEE_U05_orig'].notna()
N_winter_removed_u05 = mask_winter_orig.sum() - mask_winter_kept_u05.sum()
winter_removed_frac_u05 = 100.0 * N_winter_removed_u05 / N_night_all

mask_spring_orig = base_night & spring_flag
mask_spring_kept_u05 = mask_spring_orig & data1['NEE_U05_orig'].notna()
N_spring_removed_u05 = mask_spring_orig.sum() - mask_spring_kept_u05.sum()
spring_removed_frac_u05 = 100.0 * N_spring_removed_u05 / N_night_all

mask_summer_orig = base_night & summer_flag
mask_summer_kept_u05 = mask_summer_orig & data1['NEE_U05_orig'].notna()
N_summer_removed_u05 = mask_summer_orig.sum() - mask_summer_kept_u05.sum()
summer_removed_frac_u05 = 100.0 * N_summer_removed_u05 / N_night_all

mask_autumn_orig = base_night & autumn_flag
mask_autumn_kept_u05 = mask_autumn_orig & data1['NEE_U05_orig'].notna()
N_autumn_removed_u05 = mask_autumn_orig.sum() - mask_autumn_kept_u05.sum()
autumn_removed_frac_u05 = 100.0 * N_autumn_removed_u05 / N_night_all

winter_removed_u05 = winter_removed_frac_u05
spring_removed_u05 = spring_removed_frac_u05
summer_removed_u05 = summer_removed_frac_u05
autumn_removed_u05 = autumn_removed_frac_u05

# U16
mask_winter_kept_u16 = mask_winter_orig & data1['NEE_U16_orig'].notna()
N_winter_removed_u16 = mask_winter_orig.sum() - mask_winter_kept_u16.sum()
winter_removed_frac_u16 = 100.0 * N_winter_removed_u16 / N_night_all

mask_spring_kept_u16 = mask_spring_orig & data1['NEE_U16_orig'].notna()
N_spring_removed_u16 = mask_spring_orig.sum() - mask_spring_kept_u16.sum()
spring_removed_frac_u16 = 100.0 * N_spring_removed_u16 / N_night_all

mask_summer_kept_u16 = mask_summer_orig & data1['NEE_U16_orig'].notna()
N_summer_removed_u16 = mask_summer_orig.sum() - mask_summer_kept_u16.sum()
summer_removed_frac_u16 = 100.0 * N_summer_removed_u16 / N_night_all

mask_autumn_kept_u16 = mask_autumn_orig & data1['NEE_U16_orig'].notna()
N_autumn_removed_u16 = mask_autumn_orig.sum() - mask_autumn_kept_u16.sum()
autumn_removed_frac_u16 = 100.0 * N_autumn_removed_u16 / N_night_all

winter_removed_u16 = winter_removed_frac_u16
spring_removed_u16 = spring_removed_frac_u16
summer_removed_u16 = summer_removed_frac_u16
autumn_removed_u16 = autumn_removed_frac_u16

# u25
mask_winter_kept_u25 = mask_winter_orig & data1['NEE_U25_orig'].notna()
N_winter_removed_u25 = mask_winter_orig.sum() - mask_winter_kept_u25.sum()
winter_removed_frac_u25 = 100.0 * N_winter_removed_u25 / N_night_all

mask_spring_kept_u25 = mask_spring_orig & data1['NEE_U25_orig'].notna()
N_spring_removed_u25 = mask_spring_orig.sum() - mask_spring_kept_u25.sum()
spring_removed_frac_u25 = 100.0 * N_spring_removed_u25 / N_night_all

mask_summer_kept_u25 = mask_summer_orig & data1['NEE_U25_orig'].notna()
N_summer_removed_u25 = mask_summer_orig.sum() - mask_summer_kept_u25.sum()
summer_removed_frac_u25 = 100.0 * N_summer_removed_u25 / N_night_all

mask_autumn_kept_u25 = mask_autumn_orig & data1['NEE_U25_orig'].notna()
N_autumn_removed_u25 = mask_autumn_orig.sum() - mask_autumn_kept_u25.sum()
autumn_removed_frac_u25 = 100.0 * N_autumn_removed_u25 / N_night_all

winter_removed_u25 = winter_removed_frac_u25
spring_removed_u25 = spring_removed_frac_u25
summer_removed_u25 = summer_removed_frac_u25
autumn_removed_u25 = autumn_removed_frac_u25

# u50
mask_winter_kept_u50 = mask_winter_orig & data1['NEE_U50_orig'].notna()
N_winter_removed_u50 = mask_winter_orig.sum() - mask_winter_kept_u50.sum()
winter_removed_frac_u50 = 100.0 * N_winter_removed_u50 / N_night_all

mask_spring_kept_u50 = mask_spring_orig & data1['NEE_U50_orig'].notna()
N_spring_removed_u50 = mask_spring_orig.sum() - mask_spring_kept_u50.sum()
spring_removed_frac_u50 = 100.0 * N_spring_removed_u50 / N_night_all

mask_summer_kept_u50 = mask_summer_orig & data1['NEE_U50_orig'].notna()
N_summer_removed_u50 = mask_summer_orig.sum() - mask_summer_kept_u50.sum()
summer_removed_frac_u50 = 100.0 * N_summer_removed_u50 / N_night_all

mask_autumn_kept_u50 = mask_autumn_orig & data1['NEE_U50_orig'].notna()
N_autumn_removed_u50 = mask_autumn_orig.sum() - mask_autumn_kept_u50.sum()
autumn_removed_frac_u50 = 100.0 * N_autumn_removed_u50 / N_night_all

winter_removed_u50 = winter_removed_frac_u50
spring_removed_u50 = spring_removed_frac_u50
summer_removed_u50 = summer_removed_frac_u50
autumn_removed_u50 = autumn_removed_frac_u50

# 75
mask_winter_kept_u75 = mask_winter_orig & data1['NEE_U75_orig'].notna()
N_winter_removed_u75 = mask_winter_orig.sum() - mask_winter_kept_u75.sum()
winter_removed_frac_u75 = 100.0 * N_winter_removed_u75 / N_night_all

mask_spring_kept_u75 = mask_spring_orig & data1['NEE_U75_orig'].notna()
N_spring_removed_u75 = mask_spring_orig.sum() - mask_spring_kept_u75.sum()
spring_removed_frac_u75 = 100.0 * N_spring_removed_u75 / N_night_all

mask_summer_kept_u75 = mask_summer_orig & data1['NEE_U75_orig'].notna()
N_summer_removed_u75 = mask_summer_orig.sum() - mask_summer_kept_u75.sum()
summer_removed_frac_u75 = 100.0 * N_summer_removed_u75 / N_night_all

mask_autumn_kept_u75 = mask_autumn_orig & data1['NEE_U75_orig'].notna()
N_autumn_removed_u75 = mask_autumn_orig.sum() - mask_autumn_kept_u75.sum()
autumn_removed_frac_u75 = 100.0 * N_autumn_removed_u75 / N_night_all

winter_removed_u75 = winter_removed_frac_u75
spring_removed_u75 = spring_removed_frac_u75
summer_removed_u75 = summer_removed_frac_u75
autumn_removed_u75 = autumn_removed_frac_u75

# 84
mask_winter_kept_u84 = mask_winter_orig & data1['NEE_U84_orig'].notna()
N_winter_removed_u84 = mask_winter_orig.sum() - mask_winter_kept_u84.sum()
winter_removed_frac_u84 = 100.0 * N_winter_removed_u84 / N_night_all

mask_spring_kept_u84 = mask_spring_orig & data1['NEE_U84_orig'].notna()
N_spring_removed_u84 = mask_spring_orig.sum() - mask_spring_kept_u84.sum()
spring_removed_frac_u84 = 100.0 * N_spring_removed_u84 / N_night_all

mask_summer_kept_u84 = mask_summer_orig & data1['NEE_U84_orig'].notna()
N_summer_removed_u84 = mask_summer_orig.sum() - mask_summer_kept_u84.sum()
summer_removed_frac_u84 = 100.0 * N_summer_removed_u84 / N_night_all

mask_autumn_kept_u84 = mask_autumn_orig & data1['NEE_U84_orig'].notna()
N_autumn_removed_u84 = mask_autumn_orig.sum() - mask_autumn_kept_u84.sum()
autumn_removed_frac_u84 = 100.0 * N_autumn_removed_u84 / N_night_all

winter_removed_u84 = winter_removed_frac_u84
spring_removed_u84 = spring_removed_frac_u84
summer_removed_u84 = summer_removed_frac_u84
autumn_removed_u84 = autumn_removed_frac_u84

# 95
mask_winter_kept_u95 = mask_winter_orig & data1['NEE_U95_orig'].notna()
N_winter_removed_u95 = mask_winter_orig.sum() - mask_winter_kept_u95.sum()
winter_removed_frac_u95 = 100.0 * N_winter_removed_u95 / N_night_all

mask_spring_kept_u95 = mask_spring_orig & data1['NEE_U95_orig'].notna()
N_spring_removed_u95 = mask_spring_orig.sum() - mask_spring_kept_u95.sum()
spring_removed_frac_u95 = 100.0 * N_spring_removed_u95 / N_night_all

mask_summer_kept_u95 = mask_summer_orig & data1['NEE_U95_orig'].notna()
N_summer_removed_u95 = mask_summer_orig.sum() - mask_summer_kept_u95.sum()
summer_removed_frac_u95 = 100.0 * N_summer_removed_u95 / N_night_all

mask_autumn_kept_u95 = mask_autumn_orig & data1['NEE_U95_orig'].notna()
N_autumn_removed_u95 = mask_autumn_orig.sum() - mask_autumn_kept_u95.sum()
autumn_removed_frac_u95 = 100.0 * N_autumn_removed_u95 / N_night_all

winter_removed_u95 = winter_removed_frac_u95
spring_removed_u95 = spring_removed_frac_u95
summer_removed_u95 = summer_removed_frac_u95
autumn_removed_u95 = autumn_removed_frac_u95

# put in a lst
import numpy as np
import matplotlib.pyplot as plt

x_label = ['U05', 'U16', 'U25', 'U50', 'U75', 'U84', 'U95']
x = np.arange(len(x_label))

winter_lst = np.array([
    winter_removed_u05,
    winter_removed_u16,
    winter_removed_u25,
    winter_removed_u50,
    winter_removed_u75,
    winter_removed_u84,
    winter_removed_u95,
])

spring_lst = np.array([
    spring_removed_u05,
    spring_removed_u16,
    spring_removed_u25,
    spring_removed_u50,
    spring_removed_u75,
    spring_removed_u84,
    spring_removed_u95,
])

summer_lst = np.array([
    summer_removed_u05,
    summer_removed_u16,
    summer_removed_u25,
    summer_removed_u50,
    summer_removed_u75,
    summer_removed_u84,
    summer_removed_u95,
])

autumn_lst = np.array([
    autumn_removed_u05,
    autumn_removed_u16,
    autumn_removed_u25,
    autumn_removed_u50,
    autumn_removed_u75,
    autumn_removed_u84,
    autumn_removed_u95,
])


if drawfigs[1]:
    al = 0.4
    fig0 = plt.figure(dpi=300, figsize=figure_size, clear=True)
    ax = fig0.add_subplot(1, 1, 1)
    bottom = np.zeros(len(x_label))
    ax.bar(x, winter_lst, bottom=bottom, label='Winter(DJF)',alpha=al)
    bottom += winter_lst

    ax.bar(x, spring_lst, bottom=bottom, label='Spring(MAM)',alpha=al)
    bottom += spring_lst

    ax.bar(x, summer_lst, bottom=bottom, label='Summer(JJA)',alpha=al)
    bottom += summer_lst

    ax.bar(x, autumn_lst, bottom=bottom, label='Autumn(SON)',alpha=al)
    bottom += autumn_lst

    ax.set_xticks(x)
    ax.set_xticklabels(x_label)
    ax.set_ylabel('Removed nighttime records of all nighttime records(%)', fontsize=16)
    ax.set_xlabel('uStar thresholds', fontsize=16)
    ax.grid(axis='y', alpha=al)
    ax.legend(fontsize=16)

    plt.savefig(output_path + '/' + figure_name[1] + '.jpg', dpi=300,
                format='jpg')










