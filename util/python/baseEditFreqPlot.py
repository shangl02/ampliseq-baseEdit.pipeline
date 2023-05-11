import sys
import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot  as plt


def plot(infile, outfile):
    ## load freq table, subset, transpose and *100
    td=pd.read_table(infile, delimiter=',')
    data=td[['%A','%C','%G','%T']].T
    data=data.mul(100)  
    pos=td['start'].T

    ## rename columns
    refNts=td.originally.values
    for index, value in enumerate(refNts):
        refNts[index] = str(pos[index]) +' ('+value[0]+')'
    data.columns=refNts
    #print(data)

    ## heatmap
    fig, ax = plt.subplots(figsize=(20,3))         # Sample fig size in inches
    #labels=data.applymap(lambda v:v if v!=0 else '')
    plot = sns.heatmap(data, 
                       annot=True,fmt=".2f",
                       linewidth=.5,
                       #cmap=sns.diverging_palette(145, 300, s=60, as_cmap=True),
                       cmap=sns.color_palette("light:b", as_cmap=True),
                       square=True)
    fig = plot.get_figure()
    fig.savefig(outfile)

if __name__ == "__main__":
    try:
        #dir=r'X:\users\shangl02\Test\ampliCan.test2\AS15\output\vs1\edit_percent'
        dir=sys.argv[1]
        print(dir)
        for file in os.listdir(dir):
            if file.endswith(".csv"):
                print(file)
                infile = os.path.join(dir,file)
                plotFile = os.path.join(dir, file+'.freqHeatmap.png')
                plot(infile, plotFile)

        print("The output file generation is complete.")  
    
    except:
        print("The output file generation has failed. Kindly check.")
