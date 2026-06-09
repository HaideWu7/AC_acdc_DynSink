#!/bin/bash

matlab -batch "run_acdc_file('preprocessed_all_data.csv','DynamicCS', true, 'SMPSFile', 'SMPS_example.csv', 'SMPSLogBase', 'ln')"
