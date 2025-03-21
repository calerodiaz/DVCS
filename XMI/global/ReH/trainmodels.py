import sys
import ROOT
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from tensorflow_addons.metrics import RSquare
import matplotlib.pyplot as plt

def get_input_data(irep):
    replica_name = f'./replicas/replica_{irep}_ReH_xmi_data_5pct_sets_21to25.root'
    df = ROOT.RDataFrame("sampled_ReH", replica_name) 
    dfcols = df.AsNumpy()
    kin_array = np.array(dfcols['t'])
    # kin_array = np.dstack((dfcols['xB'], dfcols['t']))
    # kin_array = kin_array.reshape(kin_array.shape[1:])
    ReH_array = np.array(dfcols['ReH'])
    #print("kinematics array: \n",kin_array)
    print("kinematics array shape: \n", kin_array.shape)
    #print("ReH array: \n",ReH_array)
    print("ReH array shape: \n", ReH_array.shape)
    return kin_array, ReH_array

def build_models():
    model = tf.keras.Sequential([
    tf.keras.layers.Flatten(input_shape=(1,)),
    tf.keras.layers.Dense(20, activation='relu'),    
    tf.keras.layers.Dense(40, activation='tanh'),
    tf.keras.layers.Dense(180, activation='sigmoid'),
    tf.keras.layers.Dense(130, activation='relu'),
    tf.keras.layers.Dense(20, activation='relu'),
    tf.keras.layers.Dense(1, activation='linear')
    ])  

    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate = 0.0001),             
                  loss=tf.keras.losses.MeanSquaredError(reduction='sum'),
                  metrics=['MeanSquaredError', RSquare()]) 
    return model

def fit_data(replica_num, kin_array, ReH_array):
    X_train, X_test, y_train, y_test = train_test_split(kin_array, ReH_array, test_size=0.2, random_state=42)
    # Create a callback to stop training early after reaching a certain value for the validation loss.
    stop_early = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=50, verbose=1, mode='auto') # used it for models 6 with patience = 10
    # model_names = ['model']
    batch_size = 5 # model 7, 8
    # history = [None] * len(model_names)    
    models = build_models()
    # for idx, model_name in enumerate(model_names): 
    models.summary()
    history = models.fit(X_train, y_train, epochs=1000, batch_size=batch_size, validation_data=(X_test, y_test), callbacks=[stop_early]) 
    models.save('models/gfit_replica_'+str(replica_num)+'.keras') 
    np.save('models/history_gfit_replica_'+str(replica_num)+'.npy',history.history)

for irep in range(1,101):
    kin_array, ReH_array = get_input_data(irep)
    fit_data(irep, kin_array, ReH_array)

# Trial 02 summary
# Hyperparameters:
# layers: 5
# units_0: 20
# act_0: relu
# learning_rate: 0.0001
# reduction: sum
# delta: 14
# loss: MSE
# optimizer: Adam
# units_1: 40
# act_1: tanh
# units_2: 180
# act_2: sigmoid
# units_3: 130
# act_3: relu
# units_4: 20
# act_4: relu
# units_5: 20
# act_5: sigmoid
# batch_size: 5
# Score: 5.791136118205031e-06

# Trial 03 summary
# Hyperparameters:
# layers: 4
# units_0: 90
# act_0: relu
# learning_rate: 0.01
# reduction: sum_over_batch_size
# delta: 24
# loss: Huber
# optimizer: Adam
# units_1: 190
# act_1: relu
# units_2: 200
# act_2: tanh
# units_3: 130
# act_3: tanh
# units_4: 180
# act_4: tanh
# units_5: 70
# act_5: tanh
# batch_size: 9
# Score: 7.033996553218458e-06