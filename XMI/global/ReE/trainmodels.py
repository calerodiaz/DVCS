import sys
import ROOT
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from tensorflow_addons.metrics import RSquare
import matplotlib.pyplot as plt

def get_input_data(irep):
    replica_name = f'./replicas/replica_{irep}_ReE_xmi_data_5pct_sets_21to25.root'
    df = ROOT.RDataFrame("sampled_ReE", replica_name) 
    dfcols = df.AsNumpy()
    kin_array = np.array(dfcols['t'])
    # kin_array = np.dstack((dfcols['xB'], dfcols['t']))
    # kin_array = kin_array.reshape(kin_array.shape[1:])
    ReE_array = np.array(dfcols['ReE'])
    #print("kinematics array: \n",kin_array)
    print("kinematics array shape: \n", kin_array.shape)
    #print("ReE array: \n",ReE_array)
    print("ReE array shape: \n", ReE_array.shape)
    return kin_array, ReE_array

def build_models():
    model = tf.keras.Sequential([
    tf.keras.layers.Flatten(input_shape=(1,)),
    tf.keras.layers.Dense(140, activation='tanh'),    
    tf.keras.layers.Dense(140, activation='tanh'),
    tf.keras.layers.Dense(70, activation='relu'),
    tf.keras.layers.Dense(200, activation='sigmoid'),
    tf.keras.layers.Dense(120, activation='tanh'),
    tf.keras.layers.Dense(170, activation='relu'),
    tf.keras.layers.Dense(1, activation='linear')
    ])  

    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate = 0.0001),             
                  loss=tf.keras.losses.Huber(reduction='none', delta=41),
                  metrics=['MeanSquaredError', RSquare()]) 
    return model

def fit_data(replica_num, kin_array, ReE_array):
    X_train, X_test, y_train, y_test = train_test_split(kin_array, ReE_array, test_size=0.2, random_state=42)
    # Create a callback to stop training early after reaching a certain value for the validation loss.
    stop_early = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=50, verbose=1, mode='auto') # used it for models 6 with patience = 10
    # model_names = ['model']
    batch_size = 6 # model 7, 8
    # history = [None] * len(model_names)    
    models = build_models()
    # for idx, model_name in enumerate(model_names): 
    models.summary()
    history = models.fit(X_train, y_train, epochs=1000, batch_size=batch_size, validation_data=(X_test, y_test), callbacks=[stop_early]) 
    models.save('models/gfit_replica_'+str(replica_num)+'.keras') 
    np.save('models/history_gfit_replica_'+str(replica_num)+'.npy',history.history)

for irep in range(1,101):
    kin_array, ReE_array = get_input_data(irep)
    fit_data(irep, kin_array, ReE_array)


# Hyperparameters:
# layers: 6
# units_0: 140
# act_0: tanh
# learning_rate: 0.0001
# reduction: none
# delta: 41
# loss: Huber
# optimizer: Adam
# units_1: 140
# act_1: tanh
# units_2: 70
# act_2: relu
# units_3: 200
# act_3: sigmoid
# batch_size: 6
# units_4: 120
# act_4: tanh
# units_5: 170
# act_5: relu
# Score: 2.5084686967602465e-06

# Trial 04 summary
# Hyperparameters:
# layers: 6
# units_0: 140
# act_0: tanh
# learning_rate: 0.01
# reduction: none
# delta: 19
# loss: Huber
# optimizer: Adam
# units_1: 170
# act_1: tanh
# units_2: 170
# act_2: tanh
# units_3: 190
# act_3: relu
# batch_size: 3
# units_4: 200
# act_4: tanh
# units_5: 120
# act_5: relu
# Score: 1.4604714124288876e-05
