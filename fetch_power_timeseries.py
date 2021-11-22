#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 17:15:07 2018

@author: qinggang
"""

from nilearn import image
from nilearn import datasets
from nilearn import input_data
from nilearn.connectome import ConnectivityMeasure
import numpy as np
sub1_fname = "/home/qinggang/research/MIDUS_refresher/resting_conn/conn_MR_new/results/preprocessing/niftiDATA_Subject001_Condition000.nii"
sub1_data = image.load_img(sub1_fname)
power = datasets.fetch_coords_power_2011()
coords = np.vstack((power.rois['x'], power.rois['y'], power.rois['z'])).T
spheres_masker = input_data.NiftiSpheresMasker(
        seeds = coords, radius = 5., allow_overlap=False)
timeseries2 = spheres_masker.fit_transform(sub1_fname)

all_sub = "/home/qinggang/research/MIDUS_refresher/resting_conn/conn_MR_new/results/preprocessing/niftiDATA_Subject00*_Condition000.nii"
all_sub_data = image.load_img(all_sub)


all_subid_num = list(range(120))[1:]
all_subid_str = []
for i in all_subid_num:
    item = str(i)
    all_subid_str.append(item)
all_subid = []
for i in all_subid_str:
    if len(i) == 1:
        all_subid.append('00'+ i)
    elif len(i) == 2:
        all_subid.append('0' + i)
    else:
        all_subid.append(i)

timeseries_list2 = []
for sub in all_subid:
    sub_fname = "/home/qinggang/research/MIDUS_refresher/resting_conn/" + \
    "conn_MR_new/results/preprocessing/niftiDATA_Subject" + str(sub) + \
    "_Condition000.nii"
    timeseries = spheres_masker.fit_transform(sub_fname)
    timeseries_list2.append(timeseries)

corr_measure = ConnectivityMeasure(kind = 'correlation')
corr_matrix = corr_measure.fit_transform(timeseries_list2)