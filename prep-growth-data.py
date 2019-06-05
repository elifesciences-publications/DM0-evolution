'''
prep-growth-data.py by Rohan Maddamsetti.

This script loads the well labels, and the raw time series data
from the plate reader.

It then puts these data into 'tidy' format in a data frame.
'tidy' is the format that Hadley Wickham uses.

Then write the growth curve experiment data into one table
in the results directory. Data in this table will be used by make-figures.R.
'''

from os.path import join
import pandas as pd

def read_plate_labels(f):
    ''' return a tidy pandas dataframe with the well layout.'''
    fh = open(f)
    cols = { 'Name':[], 'Well':[] }
    for i,l in enumerate(fh):
        l = l.strip()
        if i == 0:
            continue
        data = l.split(',')
        row = data[0]
        names = data[1:]
        column = [str(i) for i in [1,2,3,4,5,6,7,8,9,10,11,12]]
        well = [row+x for x in column]
        assert len(well) == len(names)
        for i in range(len(well)):
            cols['Well'].append(well[i])
            cols['Name'].append(names[i])
    return pd.DataFrame(cols)

def read_growth_data(f):
    ''' return a tidy pandas dataframe of the growth curve time series.'''
    cols = {'Time':[], 'Temperature':[],'Well':[],'OD420':[]}
    fh = open(f)
    well = []
    for i,l in enumerate(fh):
        l = l.strip()
        if i == 0:
            well = l.split(',')[2:]
            continue
        data = l.split(',')
        time = data[0]
        temp = data[1]
        OD420 = data[2:]
        assert len(well) == len(OD420)
        for i in range(len(well)):
            cols['Time'].append(time)
            cols['Temperature'].append(temp)
            cols['Well'].append(well[i])
            cols['OD420'].append(float(OD420[i]))
    return pd.DataFrame(cols)

def prep_growth_data(plate_layout_file, raw_DM0_growth_file, raw_DM25_growth_file, outf):
    proj_dir = "/Users/Rohandinho/BoxSync/active-projects/DM0-evolution/"
    data_dir = "data/rohan-formatted/"

    wellfile = join(proj_dir,data_dir, plate_layout_file)
    plate_labels = read_plate_labels(wellfile)

    raw_DM0_data = join(proj_dir,data_dir, raw_DM0_growth_file)
    tidy_DM0_data = read_growth_data(raw_DM0_data)
    ## label the samples using plate_labels.
    tidy_DM0_data = pd.merge(tidy_DM0_data,plate_labels,how='outer',on='Well')

    raw_DM25_data = join(proj_dir,data_dir, raw_DM25_growth_file)
    tidy_DM25_data = read_growth_data(raw_DM25_data)
    ## label the samples. using plate_labels.
    tidy_DM25_data = pd.merge(tidy_DM25_data,plate_labels,how='outer',on='Well')

    ## label by experiment and merge the two data frames.
    tidy_DM0_data['Experiment'] = pd.Series(['DM0-growth' for x in range(len(tidy_DM0_data.index))], index=tidy_DM0_data.index)
    tidy_DM25_data['Experiment'] = pd.Series(['DM25-growth' for x in range(len(tidy_DM25_data.index))], index=tidy_DM25_data.index)

    full_table = pd.concat([tidy_DM0_data,tidy_DM25_data])
    ## now write to results directory.
    outfile = join(proj_dir,"results/", outf)
    full_table.to_csv(outfile)


def main():

    ## DM0-evolved population data
    prep_growth_data("DM0-evolved-pop-growth-plate-layout.csv",
                     "DM0-evolved-pop-DM0-growth-4-18-13.csv",
                     "DM0-evolved-pop-DM25-growth-4-15-13.csv",
                     "DM0-evolved-pop-growth.csv")

    ## DM0-evolved clone data
    prep_growth_data("DM0-evolved-clone-growth-plate-layout.csv",
                     "DM0-evolved-clone-DM0-growth-5-17-19.csv",
                     "DM0-evolved-clone-DM25-growth-5-15-19.csv",
                     "DM0-evolved-clone-growth.csv")

main()
