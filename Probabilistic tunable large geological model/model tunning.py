import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import warnings
warnings.filterwarnings("ignore")

# ========= input data
DATA = pd.read_csv(r"uncertinty_3DRF.csv", header=0, index_col=0, encoding='utf-8')
Borehole = pd.read_csv(r"XIWAN_bore.csv", header=0, index_col=0, encoding='utf-8')

# ========= define SSBU
def HTM(data):
    code_order = ['MUD', 'SAND', 'CLAY', 'CDG', 'GRANITE']
    transition_matrix = pd.DataFrame(np.zeros((len(code_order), len(code_order))), index=code_order, columns=code_order)
    X=data.Number.unique()
    for number in X:
        sorted_df = data[data['Number'] == number].sort_values('X')

        for i in range(len(sorted_df) - 1):
            current_code = sorted_df.iloc[i]['pred']
            next_code = sorted_df.iloc[i + 1]['pred']
            transition_matrix.at[current_code, next_code] += 1
            
    for number in X:
        sorted_df = data[data['Number'] == number].sort_values('Y')

        for i in range(len(sorted_df) - 1):
            current_code = sorted_df.iloc[i]['pred']
            next_code = sorted_df.iloc[i + 1]['pred']
            transition_matrix.at[current_code, next_code] += 1
    
    transition_matrix = transition_matrix.div(transition_matrix.sum(axis=1), axis=0)
    transition_matrix = transition_matrix.fillna(0)

    V = np.mean(np.diagonal(transition_matrix))

    return V

def zone(HTM_):
    LIKE=[]

    for i in range(1000):
        
        likelihood = np.array([[0, 0, 1, 0, 0]])
        v=HTM_
        u=(1-v)/4
        TPT = np.array([[v, u, u, u,u],
                        [u, v, u, u,u],
                        [u, u, v, u,u],
                        [u, u, u, v ,u],
                        [u, u, u, u ,v]
                        ])
        
        TPT_m = np.linalg.matrix_power(TPT, i)
        likelihood = np.dot(likelihood, TPT_m)
        
        LIKE.append(likelihood.max())
        
        if i>2:
            if LIKE[-2]-LIKE[-1]<0.001:
                break

    return len(LIKE)

def search_TPT(target):
    for r in np.arange(0.8,0.99,0.001):
        if zone(r)>=target:
            break
    return r.round(4)

def entropy(probabilities):
    return -np.sum(prob * np.log(prob) if prob != 0 else 0 for prob in probabilities)

def find_nearest_points(empty, BORE1,nn):
    empty_values = empty[['X', 'Y', 'Dep']].astype(float).values.reshape(1, -1)
    distances = cdist(empty_values, BORE1[['X', 'Y','Dep']].values, metric='euclidean')
    nearest_indices = distances.argsort()[:, :nn]
    nearest_points = BORE1.iloc[nearest_indices.flatten()].reset_index(drop=True)

    nearest_distances = distances[:, nearest_indices.flatten()]#[:, :nn]
    nearest_points['Distance'] = nearest_distances.flatten()
    nearest_points = nearest_points.sort_values(by='Distance', ascending=False)

    code_weights = {'MUD': 0, 'SAND': 0, 'CLAY': 0, 'CDG': 0, 'GRANITE': 0}

    for i in range(len(nearest_points)):
        code = nearest_points['pred'][i]
        weight = 1 / nearest_points['Distance'][i] if nearest_points['Distance'][i] != 0 else np.inf
        code_weights[code] += weight
        
    if np.inf in code_weights.values():
        for code in code_weights:
            if code_weights[code] == np.inf:
                code_weights[code] = 1
            else:
                code_weights[code] = 0
                
    weights_sum = sum(code_weights.values())
    normalized_weights = {code: weight / weights_sum for code, weight in code_weights.items()}

    return list(normalized_weights.values())

def Bayes(empty, BORE1,TM):
    Prior = np.array(empty[["P_MUD","P_SAND","P_CLAY","P_CDG","P_GRANITE"]])
    likelihood = find_nearest_points(empty, BORE1,1)
    i=int(empty['Dis'])
    v=TM
    u=(1-v)/4
    TPT = np.array([[v, u, u, u,u],
                    [u, v, u, u,u],
                    [u, u, v, u,u],
                    [u, u, u, v ,u],
                    [u, u, u, u ,v]
                    ])
    
    TPT_m = np.linalg.matrix_power(TPT, i)
    likelihood = np.dot(likelihood, TPT_m)
    
    unnormalized_posterior = np.multiply(likelihood, Prior)
    Posterior = unnormalized_posterior / np.sum(unnormalized_posterior)

    return Posterior

# ========= run SSBU
import time
Borehole['pred']=Borehole['Code']
LOCATION=Borehole[['X','Y']].drop_duplicates()
N=DATA.drop_duplicates(['X','Y','Dep'])

for a in range(0,len(LOCATION)):#
    print('======',a)
    stime=time.time()
    XX,YY=LOCATION.iloc[a][0],LOCATION.iloc[a][1]
    Results=[]
    HTM_=0.974
    resolution=4
    LEN=zone(HTM_)
    print(LEN/resolution)

    for n in range(0,int(LEN/resolution)):
        
        BORE1=Borehole[(Borehole['X'] == XX) & (Borehole['Y'] == YY)]

        if n == 0:
            empty=N[(N['X'] == XX) & (N['Y'] == YY)]
            # BORE1=Borehole[(Borehole['X'] == XX) & (Borehole['Y'] == YY)]
        
        if n == 1:
            EE=empty.copy()
            empty = N[
                ((N['X'] >=EE['X'].iloc[0] - resolution) & (N['X'] <= EE['X'].iloc[0]+ resolution) &
                (N['Y'] >= EE['Y'].iloc[0] - resolution) & (N['Y'] <=EE['Y'].iloc[0] + resolution)) &
                ((N['X'] ==EE['X'].iloc[0] - resolution) | (N['X'] == EE['X'].iloc[0]+ resolution) |
                (N['Y'] == EE['Y'].iloc[0] - resolution) | (N['Y'] ==EE['Y'].iloc[0] + resolution))
            ]
        
        if n > 1:
            EE=empty.copy()
            empty = N[
                ((N['X'] >= EE['X'].min() - resolution) & (N['X'] <= EE['X'].max() + resolution) &
                (N['Y'] >=  EE['Y'].min() - resolution) & (N['Y'] <= EE['Y'].max() + resolution)) &
                ((N['X'] == EE['X'].min() - resolution) | (N['X'] == EE['X'].max() + resolution) |
                (N['Y'] ==  EE['Y'].min() - resolution) | (N['Y'] == EE['Y'].max() + resolution))
            ]

        empty['Dis'] = np.sqrt(((empty.X-XX)**2)+((empty.Y-YY)**2))
        Pro=empty.apply(lambda row: Bayes(row, BORE1 ,HTM_), axis=1)
        empty[['Posterior_MUD','Posterior_SAND','Posterior_CLAY','Posterior_CDG','Posterior_GRANITE']]=Pro.apply(pd.Series)
        empty['Posterior_Entropy'] = empty[['Posterior_MUD','Posterior_SAND','Posterior_CLAY','Posterior_CDG','Posterior_GRANITE']].apply(entropy, axis=1)
        empty['Posterior_Code']=empty.loc[:,['Posterior_MUD','Posterior_SAND','Posterior_CLAY','Posterior_CDG','Posterior_GRANITE']].idxmax(1)
        empty['Posterior_Code']=empty['Posterior_Code'].replace({'Posterior_MUD':'MUD','Posterior_SAND':'SAND','Posterior_CLAY':'CLAY','Posterior_CDG':'CDG','Posterior_GRANITE':'GRANITE'})
        # print(n)
        
        Results.append(empty)
        
    Results_=pd.concat(Results)
    data_updated=Results_[['X','Y','Dep','Posterior_MUD','Posterior_SAND','Posterior_CLAY','Posterior_CDG','Posterior_GRANITE','Posterior_Code','Posterior_Entropy']]
    data_updated.columns=['X','Y','Dep',"P_MUD","P_SAND","P_CLAY","P_CDG","P_GRANITE",'pred','Entropy']
    data_updated=data_updated.sort_values('X').drop_duplicates(['X','Y','Dep'])
    N_=N[['X','Y','Dep',"P_MUD","P_SAND","P_CLAY","P_CDG","P_GRANITE",'pred','Entropy']]

    N_.set_index(['X', 'Y', 'Dep'], inplace=True)
    data_updated.set_index(['X', 'Y', 'Dep'], inplace=True)
    N_.update(data_updated)
    N_.reset_index(inplace=True)
    
    etime=time.time()
    T=(etime-stime)
    print('computational cost:',T)

    N=N_.copy()

# ========= output
N_.to_csv('Region_3D_updated_DEM.csv')