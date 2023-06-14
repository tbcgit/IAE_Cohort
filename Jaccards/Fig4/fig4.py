#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
plt.style.use('seaborn-v0_8-paper')

# Origen de los datos
cohort = ['I', 'A', 'E']
colors = {'I':'#28d528', 'A':'#ff0808', 'E':'#4141ff'}

df = pd.read_csv('jaccard_list_TaxV2.tsv', sep='\t')


# Sample label
labels = list()
for idx, row in df.iterrows():
    labels.append(row[0][0])
#endfor
df['group'] = labels

xx = np.array(df.time_dif)
yy = np.array(df.Jaccard)
group = list(df.group)
distances = pd.DataFrame({'x':xx, 'y':yy, 'group':group})

# Binned Means and Standard Deviations

## Infants
infants = distances[distances.group == 'I'].copy()

# Bins
xmin = infants.x.min()
xmax = infants.x.max()
steps = int((xmax - xmin)/60)
bins = np.linspace(xmin, xmax, steps, dtype=np.float64)
#
# Group by bin
infants['bins'] = pd.cut(infants['x'], bins)
binned = infants.groupby('bins')

i_x_mean = list()
i_x_std = list()
i_y_mean = list()
i_y_std = list()
for b, g in binned:
    if len(g) != 0:
        i_x_mean.append(g.x.mean())
        i_y_mean.append(g.y.mean())
        if len(g) != 1:
            i_x_std.append(g.x.std())
            i_y_std.append(g.y.std())
        #endif
    #endif
#endfor
i_xmean = np.array(i_x_mean)
i_ymean = np.array(i_y_mean)
i_xstd = np.array(i_x_std)
i_ystd = np.array(i_y_std)


## Adults
adults = distances[distances.group == 'A'].copy()

# Bins
xmin = adults.x.min()
xmax = adults.x.max()
steps = int((xmax - xmin)/60)
bins = np.linspace(xmin, xmax, steps, dtype=np.float64)
#
# Group by bin
adults['bins'] = pd.cut(adults['x'], bins)
binned = adults.groupby('bins')

a_x_mean = list()
a_x_std = list()
a_y_mean = list()
a_y_std = list()
for b, g in binned:
    if len(g) != 0:
        a_x_mean.append(g.x.mean())
        a_y_mean.append(g.y.mean())
        if len(g) != 1:
            a_x_std.append(g.x.std())
            a_y_std.append(g.y.std())
        #endif
    #endif
#endfor
a_xmean = np.array(a_x_mean)
a_ymean = np.array(a_y_mean)
a_xstd = np.array(a_x_std)
a_ystd = np.array(a_y_std)

## Elders
elders = distances[distances.group == 'E'].copy()

# Bins
xmin = elders.x.min()
xmax = elders.x.max()
steps = int((xmax - xmin)/60)
bins = np.linspace(xmin, xmax, steps, dtype=np.float64)
#
# Group by bin
elders['bins'] = pd.cut(elders['x'], bins)
binned = elders.groupby('bins')

e_x_mean = list()
e_x_std = list()
e_y_mean = list()
e_y_std = list()
for b, g in binned:
    if len(g) != 0:
        e_x_mean.append(g.x.mean())
        e_y_mean.append(g.y.mean())
        if len(g) != 1:
            e_x_std.append(g.x.std())
            e_y_std.append(g.y.std())
        #endif
    #endif
#endfor
e_xmean = np.array(e_x_mean)
e_ymean = np.array(e_y_mean)
e_xstd = np.array(e_x_std)
e_ystd = np.array(e_y_std)



# Figure 4
fig, axs = plt.subplots(3, 1, sharex=True, 
                        gridspec_kw={'height_ratios': [1, 1, 1]}, 
                        figsize=(11.0, 9.0))
# Some horizontal space between axes
fig.subplots_adjust(hspace=0.03)

##########################################################
# Infants
##########################################################
ax = axs[0]
x = infants.x
y = infants.y
color = colors['I']
lab = 'I'
ax.plot(x, y, 'o', color=color, alpha=0.60, markersize=4, 
        markeredgecolor='black', markeredgewidth=0.1, label=lab)

i_ylow = i_ymean - i_ystd
i_yupp = i_ymean + i_ystd
ax.plot(i_xmean, i_ymean, '-', color=colors['I'], lw= 2, alpha=1.0, 
        label='Average')
ax.plot(i_xmean, i_ylow, '--', color=colors['I'], lw=1, alpha=0.70,
       label='Standard deviation')
ax.plot(i_xmean, i_yupp, '--', color=colors['I'], lw=1, alpha=0.70)
ax.fill_between(i_xmean, i_ylow, i_yupp, facecolor=colors['I'], 
        alpha=0.10, interpolate=True)
ax.legend(loc='best', shadow=False, facecolor='wheat',
              frameon=1, fontsize=11)

# Major and minor ticks
ax.grid(visible=True, which='major', linewidth=2.0, linestyle='-',
            color='gray', alpha=0.25)
ax.grid(visible=True, which='minor', linewidth=1.0, linestyle='-',
            color='gray', alpha=0.25)
ax.tick_params(direction='in', which='minor', length=3)
ax.tick_params(direction='in', which='major', length=6)

ax.tick_params(axis='both', which='major', labelsize=14)
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax.set_ylabel('Jaccard index', fontsize=14, style='italic', fontweight="bold")

###########################
# Current y-limits
ylimits = ax.get_ylim()

##########################################################
# Adults
##########################################################
ax = axs[1]
ax.set_ylim(ylimits)

x = adults.x
y = adults.y
color = colors['A']
lab = 'A'
ax.plot(x, y, 'o', color=color, alpha=0.60, markersize=4, 
        markeredgecolor='black', markeredgewidth=0.1, label=lab)

a_ylow = a_ymean - a_ystd
a_yupp = a_ymean + a_ystd
ax.plot(a_xmean, a_ymean, '-', color=colors['A'], lw= 2, alpha=1.0, 
        label='Average')
ax.plot(a_xmean, a_ylow, '--', color=colors['A'], lw=1, alpha=0.70,
       label='Standard deviation')
ax.plot(a_xmean, a_yupp, '--', color=colors['A'], lw=1, alpha=0.70)
ax.fill_between(a_xmean, a_ylow, a_yupp, facecolor=colors['A'], 
        alpha=0.10, interpolate=True)
ax.legend(loc='best', shadow=False, facecolor='wheat',
              frameon=1, fontsize=11)

# Major and minor ticks
ax.grid(visible=True, which='major', linewidth=2.0, linestyle='-',
            color='gray', alpha=0.25)
ax.grid(visible=True, which='minor', linewidth=1.0, linestyle='-',
            color='gray', alpha=0.25)
ax.tick_params(direction='in', which='minor', length=3)
ax.tick_params(direction='in', which='major', length=6)

ax.tick_params(axis='both', which='major', labelsize=14)
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax.set_ylabel('Jaccard index', fontsize=14, style='italic', fontweight="bold")

##########################################################
# Elders
##########################################################
ax = axs[2]
ax.set_ylim(ylimits)

x = elders.x
y = elders.y
color = colors['E']
lab = 'E'
ax.plot(x, y, 'o', color=color, alpha=0.60, markersize=4, 
        markeredgecolor='black', markeredgewidth=0.1, label=lab)

e_ylow = e_ymean - e_ystd
e_yupp = e_ymean + e_ystd
ax.plot(e_xmean, e_ymean, '-', color=colors['E'], lw= 2, alpha=1.0, 
        label='Average')
ax.plot(e_xmean, e_ylow, '--', color=colors['E'], lw=1, alpha=0.70,
       label='Standard deviation')
ax.plot(e_xmean, e_yupp, '--', color=colors['E'], lw=1, alpha=0.70)
ax.fill_between(e_xmean, e_ylow, e_yupp, facecolor=colors['E'], 
        alpha=0.10, interpolate=True)
ax.legend(loc='best', shadow=False, facecolor='wheat',
              frameon=1, fontsize=11)

# Major and minor ticks
ax.grid(visible=True, which='major', linewidth=2.0, linestyle='-',
            color='gray', alpha=0.25)
ax.grid(visible=True, which='minor', linewidth=1.0, linestyle='-',
            color='gray', alpha=0.25)
ax.tick_params(direction='in', which='minor', length=3)
ax.tick_params(direction='in', which='major', length=6)

ax.tick_params(axis='both', which='major', labelsize=14)
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax.set_xlabel('Time lapse (days)', fontsize=14, style='italic', fontweight="bold" )
ax.set_ylabel('Jaccard index', fontsize=14, style='italic', fontweight="bold")
#
fig.tight_layout()
fig.savefig('Fig4.pdf')





