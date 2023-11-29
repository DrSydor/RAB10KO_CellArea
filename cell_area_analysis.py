'''
This script automatically takes a series of images and segments the cytosol
of cells transfected with GFP or stained in Green with CellTracker, and then
measures the cell area for that cell. Note that the cell area output is in
pixels^2 and needs to be converted to um^2 based on the objective/microscope
used.

The input of this script is a folder of exported multi-channel images from
Volocity. Each exported image is inside a folder named after the image and
contains each channel as a separate file. This script ignores the channels
other than green and assumes the green channel is named 'T00001C02Z001.tif'.
'''

from cellpose import models, plot
from skimage import io
from skimage.measure import regionprops_table
from skimage.segmentation import clear_border
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set up Cellpose model for cytosol segmentation using pre-trained model
model = models.CellposeModel(gpu=False, pretrained_model = '20231128_cyto_2')
channels=[0,0] # Grayscale cytosol and no DNA channel

# File handling
path = os.getcwd() + "/Images/"
directory = os.listdir(path)

# Create a dataframe to hold the data
df = pd.DataFrame()

for folder in directory:
    image_488 = io.imread(path + folder + "\T00001C01Z001.tif") # Cytosol- GFP

    masks, flows, styles = model.eval(image_488,
                                      diameter=146, #determined during training
                                      resample = True, # Results in more accurate boundaries
                                      channels=channels,
                                      cellprob_threshold = 0,
                                      flow_threshold = 0.4,
                                      stitch_threshold = 0,
                                      min_size = 50)

    masks = clear_border(masks) # remove any masks that touch the edges

    # There must be a nicer, more Pythonic way of determining column name, but this works too:
    if '_mycRAB10_' in folder:
        column_name = 'myc-RAB10'
    elif '_mycPMRAB10_' in folder:
        column_name = 'myc-RAB10-PM'
    elif 'RAB10KOHenle_GFP_DAPI' in folder:
        column_name = 'RAB10 KO'
    elif 'WTHenle' in folder:
        column_name = 'WT Henle'
    else: # It should never get here, but if it does, atleast we'll know the offending folder name
        column_name = 'None'
        print(folder)

    # Measure the mask properties, specifically area, mean intensity, and centroid
    df_488 = pd.DataFrame(regionprops_table(masks, intensity_image = image_488, properties = ('centroid', 'area', 'intensity_mean')))
    df_488['Condition'] = column_name
    df_488['Image'] = str(folder)

    # Drop any rows where the mean intensity is below 7000 (determined based on microscope settings for this dataset, so not universal)
    # This allows for screening of the cells for those that are transfected
    df_488 = df_488.drop(df_488[df_488.intensity_mean < 7000].index)

    df_488 = df_488[df_488.area > 500] # Remove any cell areas less than 500 (no cell would be that small; cleans up any mis-masked cells)

    # Create an annotated image of the cells to assist with identifying the masks
    plt2, ax = plt.subplots()
    ax.imshow(image_488)
    ax.scatter(df_488['centroid-1'], df_488['centroid-0'])
    for index, row in df_488.iterrows():
        ax.annotate(index, (row['centroid-1'], row['centroid-0']))
    ax.set_axis_off()
    ax.set_title('Cell mask labels')
    filename2 = 'Masked Images/' + str(folder) + '-ref.png'
    plt2.savefig(filename2, dpi = 300)

    # I am going to plot all masks and save as a way to keep records of what is going on
    fig = plt.figure(figsize=(12,5))
    plot.show_segmentation(fig, image_488, masks, flows[0], channels=channels)
    plt.title(column_name)
    plt.tight_layout()
    filename = 'Masked Images/' + str(folder) + '.png'
    plt.savefig(filename, dpi = 300)

    df = pd.concat([df, df_488], ignore_index=True) # Append dataframe to master dataframe

df.to_csv('cell_area_measurements.csv', index = False)

# Plot the data

plt.close('all')

order = ['WT Henle', 'RAB10 KO', 'myc-RAB10', 'myc-RAB10-PM']

for condition in order:
    print(str(condition) + ": " + str(df['Condition'].value_counts()[condition]))

ax = sns.boxplot(x = 'Condition',y = 'area', data = df, order=order)
sns.stripplot(data=df, x = 'Condition',y = 'area', dodge=True, ax=ax)
fig = ax.figure

fig.savefig('Cell Area Plot.pdf', dpi = 300)
