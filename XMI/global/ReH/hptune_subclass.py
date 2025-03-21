import ROOT
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D  
import keras_tuner as kt
from sklearn.model_selection import train_test_split
from tensorflow_addons.metrics import RSquare
from sklearn.preprocessing import MinMaxScaler, StandardScaler

class MyHyperModel(kt.HyperModel):
    def build(self, hp):
        model = tf.keras.Sequential()
        model.add(tf.keras.layers.Flatten(input_shape=(1,)))
        # Tune the number of units (neurons (nodes)) and the number of dense hidden layer
        # Activation is going to be a list of choices
        for i in range(hp.Int('layers', 1, 6)):
            model.add(tf.keras.layers.Dense(units=hp.Int('units_' + str(i), 20, 200, step=10),
                                            activation=hp.Choice('act_' + str(i), ['relu', 'tanh', 'sigmoid'])))        
        
        # Output layer, 1 output neurons corresponding to 1 output classes
        model.add(tf.keras.layers.Dense(1, activation='linear'))
        
        # Tune model compilation hyperparameters
        hp_learning_rate = hp.Choice('learning_rate', values=[1e-2, 1e-3, 1e-4])   
        hp_reduction = hp.Choice('reduction', values=['none', 'sum_over_batch_size', 'sum'])  
        hp_delta = hp.Int('delta', 0, 50, step=1)  
        hp_loss = hp.Choice('loss', values=['MSE', 'Huber'])  
        # hp_loss = hp.Choice('loss', values=['MSE'])  
        loss_dict = {"MSE": tf.keras.losses.MeanSquaredError(reduction=hp_reduction),
                     "Huber": tf.keras.losses.Huber(reduction=hp_reduction, delta=hp_delta)}
        # hp_optimizers = hp.Choice('optimizer', values=["Adam", "SGD", "Nadam"])
        hp_optimizers = hp.Choice('optimizer', values=["Adam"])
        optimizers_dict = {"Adam":    tf.keras.optimizers.Adam(learning_rate=hp_learning_rate),
                           "SGD":     tf.keras.optimizers.SGD(learning_rate=hp_learning_rate),
                           "Nadam":     tf.keras.optimizers.Nadam(learning_rate=hp_learning_rate),
                           "Adagrad":     tf.keras.optimizers.Adagrad(learning_rate=hp_learning_rate)
                           }
    
        model.compile(optimizer=optimizers_dict[hp_optimizers], loss=loss_dict[hp_loss],
                      metrics=['MeanSquaredError', RSquare()])
        return model

    def fit(self, hp, model, *args, **kwargs):
        return model.fit(
            *args,
            batch_size=hp.Int('batch_size', 1, 10, step=1),
            **kwargs,
        )
    def evaluate(self, model, x, y):
        return model.evaluate(x,y)
    def predict(self, model, x):
        return model.predict(x)   

#df = ROOT.RDataFrame("sampled_ReH", "/media/lily/Data/GPDs/TMVA_ANN-global/CFFs_Model_FromKM15/PseudoMapFit/sampled_ReH_XMapCF_KM15_5pct.root") #with 425000 sets (sampling ReH)
#df = ROOT.RDataFrame("XMapCF_results", "/media/lily/Data/GPDs/tf_global/TMVA_tfKeras_85sets/XMapCF_results_KM15_5pct.root") # resutlts from XMapCF with 85 sets (no sampling)
# df = ROOT.RDataFrame("cffs_km15", "cffs_km15_HallA_85sets.root") # KM15 cffs at HallA kinematics 

# df = ROOT.RDataFrame("sampled_ReH", "./replicas/sampled_ReH_xmi_data_5pct.root") # 100 replicas per set
# replica_no = 0
# # print(f"Hello, My name is {name} and I'm {age} years old.")

# # replica = f'<= {replica_no}*100 <= rdfentry_ < {replica_no+1}*100'
# df = df1.Filter(ROOT.Form("0 <= rdfentry_ < 100"))

 
# `Count` action
# The `Count` allows to retrieve the number of the entries that passed the
# filters. Here we show how the automatic selection of the column kicks
# in in case the user specifies none.
# entries1 = d.Filter(cutb1) \
#             .Filter(cutb1b2) \
#             .Count();

df = ROOT.RDF.FromCSV('./replicas/results_HallA_sets_21to25_data_5pct.csv')

df.Display().Print()
dfcols = df.AsNumpy()
# kin_array = np.dstack((dfcols['xB'], dfcols['t']))
kin_array = np.array(dfcols['t'])
# kin_array = kin_array.reshape(kin_array.shape[1:])
ReH_array = np.array(dfcols['xmi_ReH'])
print("kinematics array: \n",kin_array)
print("kinematics array shape: \n", kin_array.shape)
print("ReH array: \n",ReH_array)
print("ReH array shape: \n", ReH_array.shape)





# scaler = MinMaxScaler()
X_train, X_test, y_train, y_test = train_test_split(kin_array, ReH_array, test_size=0.2, random_state=42)
# normalize input data
# scaler.fit(X_train)
# transform training dataset
# X_train = scaler.transform(X_train)
# transform test dataset
# X_test = scaler.transform(X_test)

# Instantiate the tuner to perform the hypertuning. The Keras Tuner has four tuners available - RandomSearch, Hyperband, BayesianOptimization, and Sklearn.
tuner = kt.RandomSearch(
    MyHyperModel(),
    objective=kt.Objective("val_mean_squared_error", "min"),
    # objective=kt.Objective("val_r_square", "max"),
    max_trials=10,
    executions_per_trial = 1,
    overwrite=True,
    directory="tuner",
    project_name="tune_hypermodel_n1",
)

tuner.search_space_summary()

# Create a callback to stop training early after reaching a certain value for the validation loss.
stop_early = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=50, verbose=1, mode='auto')

# Will stop training if the "val_loss" hasn't improved in 3 epochs.
tuner.search(X_train, y_train, epochs=1000, validation_data=(X_test, y_test), callbacks=[stop_early])

# Build the model with the optimal hyperparameters and train it on the data for 50 epochs
hypermodel = MyHyperModel()
best_hps = tuner.get_best_hyperparameters()[0]
best_model = hypermodel.build(best_hps)
history = hypermodel.fit(best_hps, best_model, X_train, y_train, epochs=1000, validation_data=(X_test, y_test), callbacks=[stop_early])

val_acc_per_epoch = history.history['val_mean_squared_error']
# val_acc_per_epoch = history.history['val_r_square']
# best_epoch = val_acc_per_epoch.index(max(val_acc_per_epoch)) + 1
best_epoch = val_acc_per_epoch.index(min(val_acc_per_epoch)) + 1
print('Best epoch: %d' % (best_epoch,))

# Re-instantiate the hypermodel and train it with the optimal number of epochs from above.
best_hypermodel2 = hypermodel.build(best_hps)

# Retrain the model
train_history = hypermodel.fit(best_hps, best_hypermodel2, X_train, y_train, epochs=best_epoch, validation_data=(X_test, y_test))

# Print best model results
best_model.summary()
tuner.results_summary()

# # #evalute the model
# test_loss, test_acc = hypermodel.evaluate(best_hypermodel2, X_test, y_test)
# print('Test accuracy:', test_acc)

loss = train_history.history['loss']
val_loss = train_history.history['val_loss']
plt.plot(loss)
plt.plot(val_loss)
plt.legend(['model_loss', 'val_loss'])
plt.show()

# -------------------------------------------------------------------------------------

# # Kinematic input to predict ReH
# # Input file with 100K kinematic points
# filename = '/media/lily/Data/GPDs/tf_global/tf/genkinKM15cffs.root'
# treename = 'genkinKM15cffs'
# df_in = ROOT.RDataFrame(treename, filename)
# df_in_cols = df_in.AsNumpy()
# kin_array_in = np.dstack((df_in_cols['xB'], df_in_cols['t']))
# kin_array_in = kin_array_in.reshape(kin_array_in.shape[1:])

# ReH_tf = hypermodel.predict(best_hypermodel2, kin_array_in)

# print("kin in array shape: \n", kin_array_in.shape)
# print("ReH tf array shape: \n", ReH_tf.shape)

# df_plots = ROOT.RDF.MakeNumpyDataFrame({"QQ":df_in_cols['QQ'],"xB":df_in_cols['xB'],"t":df_in_cols['t'],"ReH_km15":df_in_cols['ReH_km15'],"ReH_tf":ReH_tf}) 
# df_plots.Display().Print()
# # print(ReH_tf)

# # 1D graphs ReH vs t 
# gr_ReH_vs_t_tf = df_plots.Filter("xB > 0.34 && xB < 0.35").Graph("t", "ReH_tf")
# gr_ReH_vs_t_km15 = df_plots.Filter("xB > 0.34 && xB < 0.35").Graph("t", "ReH_km15")
# gr_ReH_vs_t_tf.Sort()
# gr_ReH_vs_t_km15.Sort()

# # 2D graphs ReH vs t vs xB
# gr2D_km15 = ROOT.TGraph2D(100, df_in_cols['t'], df_in_cols['xB'], df_in_cols['ReH_km15'])
# gr2D_tf = ROOT.TGraph2D(100, df_in_cols['t'], df_in_cols['xB'], ReH_tf)
# gr2D_km15.SetTitle("KM15 ; t; xB; ReH")
# gr2D_tf.SetTitle("Keras ; t; xB; ReH")

# With mathplot
# 1D
# plt1D = plt.scatter(df_in_cols['t'][(df_in_cols['xB'] > 0.34) & (df_in_cols['xB'] < 0.35)], df_in_cols['ReH_km15'][(df_in_cols['xB'] > 0.34) & (df_in_cols['xB'] < 0.35)], color ="red")
# plt1D = plt.scatter(df_in_cols['t'][(df_in_cols['xB'] > 0.34) & (df_in_cols['xB'] < 0.35)], ReH_tf[(df_in_cols['xB'] > 0.34) & (df_in_cols['xB'] < 0.35)], color ="green")
# plt1D = plt.xlabel("$t[GeV^{2}]$")
# plt1D = plt.ylabel("ReH")
# plt1D = plt.legend(["KM15" , "Keras"])
# plt1D = plt.show()

# plt1D = plt.scatter(df_in_cols['xB'][(df_in_cols['t'] > -0.3) & (df_in_cols['t'] < -0.1)], df_in_cols['ReH_km15'][(df_in_cols['t'] > -0.3) & (df_in_cols['t'] < -0.1)], color ="red")
# plt1D = plt.scatter(df_in_cols['xB'][(df_in_cols['t'] > -0.3) & (df_in_cols['t'] < -0.1)], ReH_tf[(df_in_cols['t'] > -0.3) & (df_in_cols['t'] < -0.1)], color ="green")
# plt1D = plt.xlabel("xB")
# plt1D = plt.ylabel("ReH")
# plt1D = plt.legend(["KM15" , "Keras"])
# plt1D = plt.show()

# # Draw grphs
# mg = ROOT.TMultiGraph('mg', 'TMultiGraph')
# mg.SetTitle("xB[0.34, 0.35] ; t; ReH")  
# mg.Add(gr_ReH_vs_t_km15.GetPtr())
# mg.Add(gr_ReH_vs_t_tf.GetPtr())

# gr_ReH_vs_t_km15.SetTitle("KM15; t[GeV^{2}];ReH")
# gr_ReH_vs_t_km15.SetLineColor(ROOT.kRed)
# gr_ReH_vs_t_km15.SetLineWidth(2)
# gr_ReH_vs_t_km15.SetMarkerStyle(21)
# gr_ReH_vs_t_km15.SetMarkerSize(1)
# gr_ReH_vs_t_km15.SetMarkerColor(ROOT.kRed)

# gr_ReH_vs_t_tf.SetTitle("Keras; t[GeV^{2}];ReH")
# gr_ReH_vs_t_tf.SetLineColor(ROOT.kGreen)
# gr_ReH_vs_t_tf.SetLineWidth(2)
# gr_ReH_vs_t_tf.SetMarkerStyle(21)
# gr_ReH_vs_t_tf.SetMarkerSize(1)
# gr_ReH_vs_t_tf.SetMarkerColor(ROOT.kGreen)

# c1 = ROOT.TCanvas()
# mg.Draw("AP")
# ROOT.gPad.BuildLegend(0.78, 0.8, 0.88, 0.9)
# c1.SaveAs("ReH_vs_t_tuner_2.png")

# c2 = ROOT.TCanvas()
# gr2D_km15.Draw('tri2')
# c2.SaveAs("km15_2D_tuner_2.png")
 
# c3 = ROOT.TCanvas()
# gr2D_tf.Draw('tri2')
# c3.SaveAs("tf_2D_tuner_2.png")






# array_examples = np.load('/media/lily/Data/GPDs/tf_global/HitMatrix.npy',allow_pickle=True)
# array_labels = np.load('/media/lily/Data/GPDs/tf_global/XandYSlope.npy',allow_pickle=True)
# print("input 1: \n", kincols)
# print("input 2: \n", ReHcol)
# print("input array_examples 1 shape: \n", array_examples.shape)
# print("input array_labels 2 shape: \n", array_labels.shape)

# kincols = df.AsNumpy(exclude=["k","QQ","ReH"])
# ReHcol = df.AsNumpy(["ReH"])
# xBcol = df.AsNumpy(["xB"])
# print("Read-out of the kinematics:\n{}\n".format(kincols))
# # print("Read-out of the ReH:\n{}\n".format(ReHcol))
# print("set 0: \n", kincols['xB'][0], kincols['t'][0])

# xBarr = np.asarray(kincols['xB'][:])
# tarr = np.asarray(kincols['t'][:])
# kinarr = np.dstack((xBarr, tarr))
# ReHarr = np.asarray(ReHcol['ReH'][:])
# kinarr = kinarr.reshape(kinarr.shape[1:])
# print("\n",kinarr)
# print("input kin array 1 shape: \n", kinarr.shape)
# print("input ReH array 2 shape: \n", ReHarr.shape)
# kincols=np.array(kincols)
# print(kincols.shape)

# print("input array_examples 1: \n", array_examples)
# print("input array_labels 2: \n", array_labels)
#cols = df.AsNumpy(["xB", "t", "ReH"]) # retrieve columns "x" and "y" as NumPy arrays
# print(cols["xB"], cols["t"], cols["ReH"]) # the values of the cols dictionary are NumPy arrays
#kincols.shape
#print("input 1 shape: \n", kincols.shape())

#X_train, X_test, y_train, y_test = train_test_split(xBcol, ReHcol, test_size=0.3, random_state=42)


# 2D
# Plot the surface.
# fig_km15, ax_km15 = plt.subplots(subplot_kw={"projection": "3d"})
# surf_km15 = ax_km15.plot_surface(df_in_cols['t'], df_in_cols['xB'], df_in_cols['ReH_km15'], cmap=cm.coolwarm, linewidth=0, antialiased=False)
# fig_km15.colorbar(surf_km15, shrink=0.5, aspect=5) # Add a color bar which maps values to colors.
# plt.show()
# fig = plt.figure()
# ax = Axes3D(fig)
# surf = ax.plot_trisurf(df_in_cols['t'], df_in_cols['xB'], df_in_cols['ReH_km15'], linewidth=0.1)
# fig.colorbar(surf, shrink=0.5, aspect=5)
# plt.show()

# # Creating figure
# fig = plt.figure(figsize =(14, 9))
# ax = plt.axes(projection ='3d')
 
# # Creating plot
# ax.plot_surface(df_in_cols['t'], df_in_cols['xB'], df_in_cols['ReH_km15'])
 
# # show plot
# plt.show()

