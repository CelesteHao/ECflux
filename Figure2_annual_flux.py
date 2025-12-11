# This function is for figure2 in NEE report
# show the annual flux: NEE/H/LE and the gap filled data
# just a brief show
# Author: Hao Zhang, Date: 12/09/2025

# import modules
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates

# setting the path
input_path = 'C:/split/output.txt'
output_path = 'C:/split/figure'

# whether draw figure
drawfigs = [True, True, True, True, True, True]
figure_name = ['daily_flux_filled']
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
columns1 = ['H_f', 'H_orig', 'LE_f', 'LE_orig', 'NEE_uStar_f', 'NEE_uStar_orig'] # fig1
columns2 = ['NEE_uStar_f', 'GPP_DT_uStar','Reco_DT_uStar','GPP_uStar_f','Reco_uStar']
daily = data1[columns1].resample('D')
daily_mean  = daily.mean()
daily_count = daily.count()

daily1 = data1[columns2].resample('D')
daily_sum = daily1.sum()
daily_sum_count = daily1.count()
# daily_mean[daily_count < 24] = np.nan

# figure 0
# a figure with the H/LE/CO2 of gap and filled

if drawfigs[0]:
    # set alpha
    al = 0.3
    ticks = pd.date_range(start='2023-12-08', end='2024-12-06', freq='91D')
    fig0 = plt.figure(dpi=300, figsize=figure_size, clear=True)
    # subplot(3,1,1) for H
    ax1 = fig0.add_subplot(3,1,1)
    ax1.plot(daily_mean.index, daily_mean['H_f'],
             marker='.', linestyle='-', label='after gap-filling',color='red',linewidth=1,markersize=2)
    ax1.plot(daily_mean.index, daily_mean['H_orig'],
             marker='.', linestyle='--', label='before gap-filling',color='blue',linewidth=1,markersize=2)
    ax1.set_ylim(bottom=-50, top=250)
    ax1.set_ylabel('Sensible Heat Flux\n(W/m^2)', fontsize = 14)
    ax1.grid(True, alpha=0.5)
    ax1.set_xticks(ticks)
    ax1.tick_params(axis='x', labelbottom=False)

    # subplot(3,1,2) for LE
    ax2 = fig0.add_subplot(3,1,2)
    ax2.plot(daily_mean.index, daily_mean['LE_f'],
             marker='.', linestyle='-', label='after gap-filling',color='red',linewidth=1,markersize=2)
    ax2.plot(daily_mean.index, daily_mean['LE_orig'],
             marker='.', linestyle='--', label='before gap-filling',color='blue',linewidth=1,markersize=2)
    ax2.set_ylim(bottom=-25, top=75)
    ax2.set_ylabel('Latent Heat Flux\n(W/m^2)', fontsize = 14)
    ax2.grid(True, alpha=0.5)
    ax2.legend(fontsize=14)
    ax2.set_xticks(ticks)
    ax2.tick_params(axis='x', labelbottom=False)

    # subplot(3,1,3) for CO2
    ax3 = fig0.add_subplot(3,1,3)
    ax3.plot(daily_mean.index, daily_mean['NEE_uStar_f'],
             marker='.', linestyle='-', label='after gap-filling',color='red',linewidth=1,markersize=2)
    ax3.plot(daily_mean.index, daily_mean['NEE_uStar_orig'],
             marker='.', linestyle='--', label='before gap-filling',color='blue',linewidth=1,markersize=2)
    ax3.set_ylim(bottom=-10, top=10)
    ax3.set_ylabel('NEE\n(μmol CO2/m^2 s)', fontsize = 14)
    ax3.grid(True, alpha=0.5)
    ax3.set_xlabel('Date', fontsize = 14)
    ax3.set_xticks(ticks)
    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%y'))
    plt.savefig(output_path+'/'+figure_name[0]+'.jpg' , dpi=300,
                format='jpg')

    # show the value of filling
    print('H before gap_filling:', np.nanmean(data['H_orig']))
    print('H after gap_filling:', np.nanmean(data['H_f']))
    print('LE before gap_filling:', np.nanmean(data['LE_orig']))
    print('LE after gap_filling:', np.nanmean(data['LE_f']))
    print('NEE before gap_filling:', np.nanmean(data['NEE_uStar_orig']))
    print('NEE after gap_filling:', np.nanmean(data['NEE_uStar_f']))
    print('NEE before gap_filling_yr:', np.nansum(data['NEE_uStar_orig'] * 1800 * 12e-6))
    print('NEE after gap_filling_yr:', np.sum(data['NEE_uStar_f']* 1800 * 12e-6))

    print('NEE before gap_filling:', data['NEE_uStar_orig'].isna().sum())
    print('NEE gap:',  data['NEE_uStar_orig'].isna().sum()/len(data['NEE_uStar_orig']))






