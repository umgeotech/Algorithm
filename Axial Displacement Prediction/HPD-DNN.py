import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from tensorflow.keras.initializers import glorot_normal
from tensorflow import keras
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.utils.np_utils import to_categorical
from keras.layers import Dense
import keras.backend as K

E1 = pd.read_csv(r"Segment-1.csv", header=0, index_col=0, encoding='utf-8').reset_index().dropna().reset_index()
E2 = pd.read_csv(r"Segment-2.csv", header=0, index_col=0, encoding='utf-8').reset_index().dropna().reset_index()
E4 = pd.read_csv(r"Segment-3.csv", header=0, index_col=0, encoding='utf-8').reset_index().dropna().reset_index()
E5 = pd.read_csv(r"Segment-4.csv", header=0, index_col=0, encoding='utf-8').reset_index().dropna().reset_index()
E13 = pd.read_csv(r"Segment-5.csv", header=0, index_col=0, encoding='utf-8').reset_index().dropna().reset_index()
E14 = pd.read_csv(r"Segment-6.csv", header=0, index_col=0, encoding='utf-8').reset_index().dropna().reset_index()
E17 = pd.read_csv(r"Segment-7.csv", header=0, index_col=0, encoding='utf-8').reset_index().dropna().reset_index()
E18 = pd.read_csv(r"Segment-8.csv", header=0, index_col=0, encoding='utf-8').reset_index().dropna().reset_index()
E24 = pd.read_csv(r"Segment-9.csv", header=0, index_col=0, encoding='utf-8').reset_index().dropna().reset_index()
E25 = pd.read_csv(r"Segment-10.csv", header=0, index_col=0, encoding='utf-8').reset_index().dropna().reset_index()
E29 = pd.read_csv(r"Segment-11.csv", header=0, index_col=0, encoding='utf-8').reset_index().dropna().reset_index()
E30 = pd.read_csv(r"Segment-12.csv", header=0, index_col=0, encoding='utf-8').reset_index().dropna().reset_index()
E33 = pd.read_csv(r"Segment-13.csv", header=0, index_col=0, encoding='utf-8').reset_index().dropna().reset_index()
Ea = pd.read_csv(r"Segment-14.csv", header=0, index_col=0, encoding='utf-8').reset_index().dropna().reset_index()

EEE=[E1,E2,E4,E5,E13,E14,E17,E18,E24,E25,E29,E30,E33,Ea]

def r2(y_true, y_pred):
    a = K.square(y_pred - y_true)
    b = K.sum(a)
    c = K.mean(y_true)
    d = K.square(y_true - c)
    e = K.sum(d)
    f = 1 - b/e
    return f

def model(Train,Test,layers,node):

    ####################### Data processing

    m=to_categorical(Train['Month'], dtype ="uint8")

    Train['M_1']=m[:,1]
    Train['M_2']=m[:,2]
    Train['M_3']=m[:,3]
    Train['M_4']=m[:,4]
    Train['M_5']=m[:,5]
    Train['M_6']=m[:,6]
    Train['M_7']=m[:,7]
    Train['M_8']=m[:,8]
    Train['M_9']=m[:,9]
    Train['M_10']=m[:,10]
    Train['M_11']=m[:,11]
    Train['M_12']=m[:,12]

    n=to_categorical(Test['Month'], dtype ="uint8")

    Test['M_1']= n[:,1]
    Test['M_2']= n[:,2]
    Test['M_3']= n[:,3]
    Test['M_4']= n[:,4]
    Test['M_5']= n[:,5]
    Test['M_6']= n[:,6]
    Test['M_7']= n[:,7]
    Test['M_8']= n[:,8]
    Test['M_9']= n[:,9]
    Test['M_10']=n[:,10]
    Test['M_11']=n[:,11]
    Test['M_12']=n[:,12]

    b=to_categorical(Train['State'], dtype ="uint8")

    Train['S_1']=b[:,1]
    Train['S_2']=b[:,2]
    Train['S_3']=b[:,3]

    a=to_categorical(Test['State'], dtype ="uint8")

    Test['S_1']=a[:,1]
    Test['S_2']=a[:,2]
    Test['S_3']=a[:,3]

    ####################### Training

    X=Train[['M_1','M_2','M_3','M_4','M_5','M_6','M_7','M_8','M_9','M_10','M_11','M_12','pre_DIS','FLX','Z','S_1','S_2','S_3']] # HPD=DNN
    y=Train[['DIS']]

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=10)

    seed=5
    n_layers = layers # Number of hidden layers
    n_nodes = node # Number of nodes per hidden layer
    model = Sequential()
    for layer in np.arange(n_layers):
        if layer == 0:
            model.add(Dense(n_nodes, kernel_initializer=glorot_normal(seed=seed),activation='relu', input_shape=(np.shape(X_train)[1],)))
        else:
            model.add(Dense(n_nodes, kernel_initializer=glorot_normal(seed=seed),activation='relu'))
        model.add(Dropout(0.015,seed=seed))
    model.add(Dense(1, kernel_initializer=glorot_normal(seed=seed),activation='linear'))

    model.compile(
        optimizer=keras.optimizers.Adadelta(lr=1, rho=0.95, epsilon=None, decay=0.0),
        loss='mae',
        metrics=['mae',r2])

    model.fit(x=X_train, y=y_train, epochs=50,batch_size=200,validation_data=(X_test, y_test),shuffle=False)

    model.predict(X_test)

    X_test=Test[['M_1','M_2','M_3','M_4','M_5','M_6','M_7','M_8','M_9','M_10','M_11','M_12','pre_DIS','FLX','Z','S_1','S_2','S_3']] # HPD-DNN
    y_test=Test[['DIS']]

    ####################### Results

    print('r2:',r2_score(y_test, model.predict(X_test)))
    print('MAE:',mean_absolute_error(y_test, model.predict(X_test)))
    print('MSE:',mean_squared_error(y_test, model.predict(X_test)))

################################ Test program

layers=3
node=10
Train=[EEE[0],EEE[1],EEE[2],EEE[3],EEE[4],EEE[5],EEE[6],EEE[7],EEE[8],EEE[9],EEE[10],EEE[11],EEE[12]] #Test generalization using Segment-14
Train=pd.concat(Train)
Test=EEE[13].copy()

model(Train,Test,layers,node)