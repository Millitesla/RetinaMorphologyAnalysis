# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 22:40:23 2016

Neuro Morphology analysis tool

@author: ruff
"""
import seaborn as sns
import neurom as nm
import pandas as pd
import glob
import math
import matplotlib.pyplot as plt
import os
import numpy as np
import pylab as p
from scipy import misc
from __future__ import print_function

from neurom.core.dataformat import COLS
import neurom as nm
from neurom import geom
from neurom.fst import sectionfunc
from neurom.core import Tree
from neurom.core.types import tree_type_checker, NEURITES
from neurom import morphmath as mm
import numpy as np
#Don't forget to convert swc files so it can be read by neurom

os.chdir('/Users/ruff/OneDrive/Retina Python Tools/DataV2') 

# Determine Stratification of neurons
def StratifyAnalysis (data, filename):
    '''
    Takes Series of images with stratification pattern (binary) and Sums each row, bins into 10, Sums with all other images 
    of same series (=same cell) and plots an image using heatmap
    '''
    binnum = 5
    allallstratify = pd.DataFrame()
    for cell in range(len(data)):
        
        images = glob.glob(os.curdir + '/InFigures/' + '*' + data['SQL_ID1'][cell] + ' ' + '*')
        allstratify=pd.DataFrame(columns={'stratify'}, index=np.linspace(0,binnum-1,binnum)).fillna(0)
        for image in images:
            stratify = pd.DataFrame()
            depth = pd.DataFrame()
            stratify = misc.imread(image, flatten = True)
            stratify = pd.DataFrame(np.sum(stratify, axis=1)).iloc[::-1]
            stratify = stratify.reset_index(drop=True)
            stratify = stratify.rename(columns = {0:'stratify'})
            depth = pd.DataFrame(np.linspace(0, 1, len(stratify)))
            depth = depth.rename(columns={0:'depth'})
            stratify['depth']=depth.depth # add depth column to stratify
            stratify.depth=stratify.depth*binnum
            stratify.depth = stratify.depth.astype(int) #convert to integer for subsequent sumation
            stratify['depth'].replace(to_replace=binnum, value=binnum-1, inplace=True, method='pad', axis=None)
            stratify = pd.pivot_table(data=stratify, columns=stratify.depth , aggfunc='sum').transpose()
            allstratify['stratify'] = stratify.stratify + allstratify.stratify
        allstratify = allstratify.stratify/allstratify.stratify.max()#normalize
        allstratify = pd.DataFrame(allstratify) #convert to dataframe
        allstratify = allstratify.transpose()
        #allstratify['SQL_ID1'] = image[12:30]
        allallstratify = pd.concat([allallstratify, allstratify])
    allallstratify = allallstratify.rename(columns = {0.0:'Layer1',1.0:'Layer2',2.0:'Layer3',3.0:'Layer4',4.0:'Layer5'})
    #allallstratify = allallstratify.rename(columns = {0.0:'Layer1',1.0:'Layer2',2.0:'Layer3',3.0:'Layer4',4.0:'Layer5',5.0:'Layer6',6.0:'Layer7',7.0:'Layer8',8.0:'Layer9',9.0:'Layer10'})#Use with 10 layers
    allallstratify = allallstratify.reset_index(drop=True)    
    #return(allallstratify)

#Filter for celltype (choose either RGC or Amacrine) #############################
data = pd.read_excel(os.curdir + '/Alldata_Project_Retina.xlsx', sheetname='CellMorphology')
data = data[data['Year'].isin([2017])]
data = data[data['Experiment'].isin([14,22,26])]
data = data[data['Sub_Type'].isin(['RGC'])]
data = data.reset_index(drop=True)

stratification = StratifyAnalysis (data, ' RGCs')   #Stratification properties

data = pd.read_excel(os.curdir + '/Alldata_Project_Retina.xlsx', sheetname='CellMorphology')
data = data[data['Year'].isin([2017])]
data = data[data['Experiment'].isin([14,22,26])]
data = data[data['Sub_Type'].isin(['Amacrine'])]
data = data.reset_index(drop=True)

stratification = StratifyAnalysis (data, ' Amacrines')   #Stratification properties

##############################################

PcaData = pd.DataFrame()                    #Make empty DataFrame for Principle component analysis
PcaData['SQL_ID1'] = data.loc[:,'SQL_ID1']                                  #Add Filename SQL_ID in first column
PcaData['SomaSize'] = data.loc[:,'Area']    #Add Soma Size
PcaData = pd.concat([PcaData, stratification], axis=1)#Add Stratification
PcaData = PcaData.set_index(PcaData.SQL_ID1, drop=False)

####SORT DATA according to stratification pattern##########
#PcaData10 = PcaData[PcaData.Layer10==1]
#PcaData9 = PcaData[PcaData.Layer9==1]
#PcaData8 = PcaData[PcaData.Layer8==1]
#PcaData7 = PcaData[PcaData.Layer7==1]
#PcaData6 = PcaData[PcaData.Layer6==1]    

PcaData5 = PcaData[PcaData.Layer5==1]
PcaData4 = PcaData[PcaData.Layer4==1]
PcaData3 = PcaData[PcaData.Layer3==1]
PcaData2 = PcaData[PcaData.Layer2==1]
PcaData1 = PcaData[PcaData.Layer1==1]    
PcaDataSorted = PcaData1.append([PcaData2, PcaData3, PcaData4, PcaData5])
#PcaDataSorted = PcaData1.append([PcaData2, PcaData3, PcaData4, PcaData5, PcaData6, PcaData7, PcaData8, PcaData9])#Use for 10 bins
###################################################

#####Plotting data#################################
plt.figure()
ax1=plt.subplot(122)
plt.title('Soma')
plt.imshow(PcaDataSorted[['SomaSize']], cmap='Reds' ,vmin=0, vmax=350, aspect=0.5)
ax2=plt.subplot(121)
plt.title('Strat')
plt.imshow(PcaDataSorted[['Layer1','Layer2','Layer3','Layer4','Layer5']], cmap='viridis', interpolation='nearest', aspect=0.5)
#plt.imshow(PcaDataSorted[['Layer1','Layer2','Layer3','Layer4','Layer5', 'Layer6','Layer7','Layer8','Layer9','Layer10']], cmap='viridis', interpolation='nearest', aspect=1.8)
#Use for 10 bins     
p.savefig(os.curdir + '/Figures/'  +'Amacrine Stratification'+'.png', frameon= False, transparent=False)
###################################################

#####Frequency distribution of SomaSize Plotted with cntrl data############
##RGC control data
ctrl = pd.read_excel(os.curdir + '/Alldata_Project_Retina.xlsx', sheetname='CellSize')
ctrl = ctrl[ctrl['Year'].isin([2016])]
ctrl = ctrl[ctrl['Experiment'].isin([89])]
ctrl = ctrl[ctrl['Type'].isin([2])]
ctrl = ctrl.reset_index(drop=True)

fig=plt.figure()
fig=plt.style.use('seaborn-white')
ax=sns.distplot(ctrl.Area, bins=5, hist=False, rug=False, kde=True, norm_hist=False, hist_kws=dict(alpha=0.5), kde_kws={"shade": True})
#SomaSize of either RGC or Amacrines
ax=sns.distplot(PcaDataSorted.SomaSize, bins=5, hist=False, rug=False, kde=True, norm_hist=False, hist_kws=dict(alpha=0.5), kde_kws={"shade": True}, color='r')

##Amacrine control data
ctrl = pd.read_excel(os.curdir + '/Alldata_Project_Retina.xlsx', sheetname='CellSize')
ctrl = ctrl[ctrl['Year'].isin([2016])]
ctrl = ctrl[ctrl['Experiment'].isin([89])]
ctrl = ctrl[ctrl['Type'].isin([1])]
ctrl = ctrl.reset_index(drop=True)

fig=plt.figure()
fig=plt.style.use('seaborn-white')
ax=sns.distplot(ctrl.Area, bins=5, hist=False, rug=False, kde=True, norm_hist=False, hist_kws=dict(alpha=0.5), kde_kws={"shade": True})
#SomaSize of either RGC or Amacrines
ax=sns.distplot(PcaDataSorted.SomaSize, bins=5, hist=False, rug=False, kde=True, norm_hist=False, hist_kws=dict(alpha=0.5), kde_kws={"shade": True}, color='r')

###################################################

#Plot Single Values like SomaSize
PcaDataSorted[['SQL_ID1','SomaSize']].plot(kind='bar')
plt.ylim(0,300)

plt.figure()
plt.scatter(x=PcaData.SomaSize, y=PcaData.Layer2, color = 'r')
plt.scatter(x=PcaData.SomaSize, y=PcaData.Layer3, color = 'g')
plt.scatter(x=PcaData.SomaSize, y=PcaData.Layer4, color = 'b')
plt.scatter(x=PcaData.SomaSize, y=PcaData.Layer5, color = 'y')


#Check for Tomato intensity

#Do PCA here
def doPCA(data):
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    pca.fit(data)
    return pca

pca = doPCA(PcaData)







   
files2 = glob.glob(os.curdir + '/*edit.swc') #load all Filenames of SWC files

for item in files2:
    nrn =  nm.load_neuron(item)
    
    
#nrn = nm.load_neuron(os.curdir + '/Image001-007-01.CNG.swc')

from neurom import viewer

fig, ax = viewer.draw(nrn)
fig.show()
fig, ax = viewer.draw(nrn, mode='3d') # valid modes '2d', '3d', 'dendrogram'
fig.show()



def sec_len(sec):
    '''Return the length of a section'''
    return mm.section_length(sec.points)

print('Total neurite length (sections):',
      sum(sec_len(s) for s in nm.iter_sections(nrn)))

# Get length of all neurites in cell by iterating over segments,
# and summing the segment lengths.
# This should yield the same result as iterating over sections.
print('Total neurite length (segments):',
      sum(mm.segment_length(s) for s in nm.iter_segments(nrn)))

# get volume of all neurites in cell by summing over segment
# volumes
print('Total neurite volume:',
      sum(mm.segment_volume(s) for s in nm.iter_segments(nrn)))

# get area of all neurites in cell by summing over segment
# areas
print('Total neurite surface area:',
      sum(mm.segment_area(s) for s in nm.iter_segments(nrn)))

# get total number of neurite points in cell.
def n_points(sec):
    '''number of points in a section'''
    n = len(sec.points)
    # Non-root sections have duplicate first point
    return n if sec.parent is None else n - 1

print('Total number of points:',
      sum(n_points(s) for s in nm.iter_sections(nrn)))

# get mean radius of neurite points in cell.
# p[COLS.R] yields the radius for point p.
# Note: this includes duplicated points at beginning of
# non-trunk sections
print('Mean radius of points:',
      np.mean([s.points[:, COLS.R] for s in nm.iter_sections(nrn)]))

# get mean radius of neurite points in cell.
# p[COLS.R] yields the radius for point p.
# Note: this includes duplicated points at beginning of
# non-trunk sections
pts = [p[COLS.R] for s in nrn.sections[1:] for p in s.points]
print('Mean radius of points:',
      np.mean(pts))

# get mean radius of segments
print('Mean radius of segments:',
      np.mean(list(mm.segment_radius(s) for s in nm.iter_segments(nrn))))

# get stats for the segment taper rate, for different types of neurite
for ttype in NEURITES:
    ttt = ttype
    seg_taper_rate = [mm.segment_taper_rate(s)
                      for s in nm.iter_segments(nrn, neurite_filter=tree_type_checker(ttt))]

    print('Segment taper rate (', ttype,
          '):\n  mean=', np.mean(seg_taper_rate),
          ', std=', np.std(seg_taper_rate),
          ', min=', np.min(seg_taper_rate),
          ', max=', np.max(seg_taper_rate),
          sep='')

# Number of bifurcation points.
print('Number of bifurcation points:',
      sum(1 for _ in nm.iter_sections(nrn,
                                      iterator_type=Tree.ibifurcation_point)))

# Number of bifurcation points for apical dendrites
print('Number of bifurcation points (apical dendrites):',
      sum(1 for _ in nm.iter_sections(nrn,
                                      iterator_type=Tree.ibifurcation_point,
                                      neurite_filter=tree_type_checker(nm.APICAL_DENDRITE))))

# Maximum branch order
print('Maximum branch order:',
      max(sectionfunc.branch_order(s) for s in nm.iter_sections(nrn)))

# Neuron's bounding box
# Note: does not account for soma radius
print('Bounding box ((min x, y, z), (max x, y, z))', geom.bounding_box(nrn))


#with open(os.curdir + '/exp14_17 Pos_R1_7 M-exported-000.swc', 'r') as f, open(os.curdir + '/new_scw.scw', 'w') as ff:
#    lines = f.readlines()
#    
#    for i, line in enumerate(lines):
#        row = line.split(' ')
#        
#        if i > 5:
#            row[-2] = '1.0'
#        
#        ff.write(' '.join(row))