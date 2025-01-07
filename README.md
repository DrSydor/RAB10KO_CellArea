# RAB10KO_CellArea
A Python script to segment GFP-transfected cells using a Cellpose model and measure the cell area of cells with a mean intensity value above a given threshold.

This script automatically takes a series of images and segments the cytosol of cells transfected with GFP or stained in Green with CellTracker, and then measures the cell area for that cell. Note that the cell area output is in pixels^2 and needs to be converted to um^2 based on the objective/microscope used.

The input of this script is a folder of exported multi-channel images from Volocity. Each exported image is inside a folder named after the image and contains each channel as a separate file. This script ignores the channels other than green and assumes the green channel is named 'T00001C02Z001.tif'.

This script was created for the paper "Salmonella exploits membrane reservoirs for invasion of host cells": 
https://www.nature.com/articles/s41467-024-47183-x
