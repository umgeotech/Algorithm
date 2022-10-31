import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
from mpl_toolkits.mplot3d import Axes3D
warnings.filterwarnings("ignore")
import matplotlib.colors as colors

mycolor=['#4A708B','#7EC0EE','#A2CD5A','#FFD700','#CD3333','#8B1A1A']
cmap_color = colors.LinearSegmentedColormap.from_list('my_list', mycolor)

def correlation(dx, dy, dd, IX, IY, ID):
    R = math.exp(-((math.pi*pow(dx,2))/pow(IX,2))-((math.pi*pow(dy,2))/pow(IY,2))-((math.pi*pow(dd,2))/pow(ID,2)))
    return R

def Probability(C):
    P = C/np.sum(C)
    return P

def Entropy(P):
    e = P*np.log(P)
    E = -np.sum(e)
    return E

def Cor(A,IX,IY,ID):
    if len(A) != 0:
        for k in range(0,len(A)):
            dx = np.array(EXP['X'  ])[0]-np.array(A['X'])[k]
            dy = np.array(EXP['Y'  ])[0]-np.array(A['Y'])[k]
            dd = np.array(EXP['Dep'])[0]-np.array(A['Dep'])[k]
            C = correlation(dx,dy,dd, IX,IY,ID)

            if k > 0:
                C_A = C+C_A
            else:
                C_A = C
    else:
        C_A = 0

    return C_A

def params(Borehole):
    XMAX = Borehole['X'].max()
    XMIN = Borehole['X'].min()
    YMAX = Borehole['Y'].max()
    YMIN = Borehole['Y'].min()
    DF=6 #depth factor

    M  = Borehole[(Borehole['Code']=='M' )]
    S  = Borehole[(Borehole['Code']=='S' )]
    R  = Borehole[(Borehole['Code']=='R' )]

    M_diff = M.groupby(['X','Y']).count()['Dep'].max()-M.groupby(['X','Y']).count()['Dep'].min()
    S_diff = S.groupby(['X','Y']).count()['Dep'].max()-S.groupby(['X','Y']).count()['Dep'].min()
    R_diff = R.groupby(['X','Y']).count()['Dep'].max()-R.groupby(['X','Y']).count()['Dep'].min()

    # SoF
    IX_M = (XMAX-XMIN)*(1-(M_diff/Borehole.groupby(['X','Y']).count()['Dep'].max()))
    IY_M = (YMAX-YMIN)*(1-(M_diff/Borehole.groupby(['X','Y']).count()['Dep'].max()))
    ID_M = (DF)*(1-(M_diff/Borehole.groupby(['X','Y']).count()['Dep'].max()))

    IX_S = (XMAX-XMIN)*(1-(S_diff/Borehole.groupby(['X','Y']).count()['Dep'].max()))
    IY_S = (YMAX-YMIN)*(1-(S_diff/Borehole.groupby(['X','Y']).count()['Dep'].max()))
    ID_S = (DF)*(1-(S_diff/Borehole.groupby(['X','Y']).count()['Dep'].max()))

    IX_R = (XMAX-XMIN)*(1-(R_diff/Borehole.groupby(['X','Y']).count()['Dep'].max()))
    IY_R = (YMAX-YMIN)*(1-(R_diff/Borehole.groupby(['X','Y']).count()['Dep'].max()))
    ID_R = (DF)*(1-(R_diff/Borehole.groupby(['X','Y']).count()['Dep'].max()))

    IX = (IX_M,IX_S,IX_R)
    IY = (IY_M,IY_S,IY_R)
    ID = (ID_M,ID_S,ID_R)

    return ID,IX,IY,M,S,R

def visualization(NEW):
    RAW = NEW.query("Y==5 or X==7 or Y==16 or X==23 or Y==25")

    fig = plt.figure(figsize=(7,7),dpi=100)
    ax=plt.subplot(111,projection='3d')

    ax.scatter(RAW['X'], RAW['Y'],RAW['Dep'],s=80,c=RAW['Entropy'],cmap=cmap_color, vmin=0.2,vmax=0.95)

    ax.set_zlabel('Z (m)', fontsize=20,labelpad=16, family='Arial')
    ax.set_ylabel('Y (m)', fontsize=20,labelpad=16, family='Arial')
    ax.set_xlabel('X (m)', fontsize=20,labelpad=16, family='Arial')
    ax.set_xticks(np.linspace(0,30,4))
    ax.set_yticks(np.linspace(0,30,4))
    ax.set_zticks(np.linspace(-18,0,3))
    ax.tick_params(labelsize=20)
    ax.view_init(elev=70., azim=-60)

    plt.show()

########### read emptry voxel data
data = pd.read_csv(r"3D_model.csv", header=0, index_col=0, encoding='utf-8')

########### define emptry voxel and borehole locations
N = data[(data['Code']=='N')]
loc = data[(data['Code']!='N')].drop_duplicates(['X','Y']).reset_index()
EXPL = N.reset_index()
k = 4

########### loop
for a in range(0,len(EXPL)):

    ########### the ath empty voxle
    EXP = EXPL[EXPL.index==a]

    ########### identify the nearest k boreholes
    X_dis = [((EXP['X']-n)**2) for n in loc['X']]
    Y_dis = [((EXP['Y']-n)**2) for n in loc['Y']]
    Dis = np.sqrt([a+b for a,b in zip(X_dis,Y_dis)])
    loc['Dis'] = Dis
    NN = np.argsort(loc['Dis'])
    NN = NN[:k]
    LOC = loc.loc[NN][['X','Y']]
    Borehole = pd.merge(LOC,data,how='left')

    ########### obtain SoF
    ID,IX,IY,M,S,R = params(Borehole)

    ########### obtain correlation
    Sum = (Cor(M,IX[0],IY[0],ID[0])+Cor(S,IX[1],IY[1],ID[1])+Cor(R,IX[2],IY[2],ID[2]))
    C_M = Cor(M,IX[0],IY[0],ID[0])+1/(Sum)
    C_S = Cor(S,IX[1],IY[1],ID[1])+1/(Sum)
    C_R = Cor(R,IX[2],IY[2],ID[2])+1/(Sum)
    Corr = (C_M, C_S, C_R)

    ########### obtain probability
    Pro = Probability(Corr)

    ########### obtain entropy
    EXP['Entropy'] = Entropy(Pro)

    ########### store dataframe
    if a == 0:
        NEW = EXP
    else:
        NEW = pd.concat([EXP,NEW])


visualization(NEW)