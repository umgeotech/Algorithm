import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

# ========= define grid model
x_range = np.arange(-20, 520, 4)
y_range = np.arange(-20, 120, 4)
depth_range = np.arange(0, 40, 1)

X, Y, Depth = np.meshgrid(x_range, y_range, depth_range, indexing='ij')
X_flat = X.flatten()
Y_flat = Y.flatten()
Depth_flat = Depth.flatten()

grid_df = pd.DataFrame({
    'X': X_flat,
    'Y': Y_flat,
    'Dep': -Depth_flat
})

# ========= input virtual BH
VB = pd.read_csv(r"VB-1.csv", header=0, index_col=0, encoding='utf-8')

# ========= define SI-RFB model
def entropy(P):
    e = P * np.log(P)
    return -np.sum(e)

def correlation(dx, dy, dd, IX, IY, ID):
    return np.exp(-((np.pi * dx**2) / IX**2) - ((np.pi * dy**2) / IY**2) - ((np.pi * dd**2) / ID**2))

def params(Borehole):
    ID_ = 3
    
    codes = ['MUD', 'SAND', 'CLAY', 'GRANITE', 'CDG']
    results = {}
    Bmax=Borehole.groupby(['X', 'Y']).count()['Dep'].max()
    
    for code in codes:
        subset = Borehole[Borehole['Code'] == code]
        code_diff = subset.groupby(['X', 'Y']).count()['Dep'].max() - subset.groupby(['X', 'Y']).count()['Dep'].min()
        
        IX = 540 * (1 - (code_diff / Bmax))
        IY = 120 * (1 - (code_diff / Bmax))
        ID = ID_*  (1-  (code_diff / Bmax))

        results[code] = (IX, IY, ID)
    
    return results

def correlation_sum(data, unsampled, params_dict):
    correlations = {}
    
    for code, (IX, IY, ID) in params_dict.items():
        subset = data[data['Code'] == code]
        dx = unsampled['X'].iloc[0] - subset['X'].values
        dy = unsampled['Y'].iloc[0] - subset['Y'].values
        dd = unsampled['Dep'].iloc[0] - subset['Dep'].values

        corrs = correlation(dx, dy, dd, IX, IY, ID)
        C_sum = np.sum(corrs)
        correlations[code] = C_sum
    
    return correlations

# ========= define model inputs
bore_=VB.copy()
loc = bore_.drop_duplicates(['X', 'Y']).reset_index()
EXPL = grid_df.drop_duplicates(['X', 'Y','Dep']).reset_index()
print(len(EXPL))

# ========= run SI-RFB
NEW = pd.DataFrame()
Bore = bore_.copy()
params_dict = params(Bore)
print(params_dict)
for a in range(len(EXPL)):
    
    EXP = EXPL[EXPL.index==a]    

    C_sum = correlation_sum(Bore, EXP, params_dict)
    total_sum = sum(C_sum.values())

    corrs = {code: C_sum[code] + 1/total_sum for code in C_sum}
    Pro = {code: corrs[code] / sum(corrs.values()) for code in corrs}
    small_value = 1e-30
    probs = [max(Pro.get(code, small_value), small_value) for code in C_sum]

    EXP[[f'P_{code}' for code in C_sum]] = probs
    EXP['Entropy'] = entropy(probs)

    NEW = pd.concat([EXP, NEW])

# ========= processing and outputing
RAW = NEW.copy()
RAW['pred']=RAW.loc[:,['P_MUD','P_SAND','P_CLAY','P_GRANITE','P_CDG']].idxmax(1)
RAW=RAW.replace({'P_MUD':'MUD','P_SAND':'SAND','P_CLAY':'CLAY','P_GRANITE':'GRANITE','P_CDG':'CDG'})
RAW[['X','Dep',	'P_MUD','P_SAND','P_CLAY','P_GRANITE','P_CDG','Entropy','pred']]

RAW.to_csv('uncertinty_3DRF.csv')