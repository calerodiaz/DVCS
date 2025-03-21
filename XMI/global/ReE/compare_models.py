import sys
import ROOT
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from tensorflow_addons.metrics import RSquare

def band_limits(a):
    amax = np.max(a, 0) # from each column (epoch)
    amin = np.min(a, 0) # from each column (epoch)
    amean = np.mean(a, 0) # from each column (epoch)
    astd = np.std(a, 0)
    return amean, amin, amax, astd

# def get_kinematics(dfcols):
#     kin_array = np.dstack((dfcols['xB'], dfcols['t']))
#     kin_array = kin_array.reshape(kin_array.shape[1:])
#     return kin_array

def load_iterations(kin_array, iterations):
    loss = []
    val_loss = []
    r_square = []
    val_r_square = []
    ReE_tf = []
    # kin_array = get_kinematics(dfcols)
    for i in range(iterations):
        model = tf.keras.models.load_model('./models/gfit_replica_'+str(i+1)+'.keras') # Load model
        history = np.load('./models/history_gfit_replica_'+str(i+1)+'.npy',allow_pickle='TRUE').item() # Load history
        loss_itr = history['loss']
        val_loss_itr = history['val_loss']
        r_square_itr = history['r_square']
        val_r_square_itr = history['val_r_square']
        loss.append(loss_itr)
        val_loss.append(val_loss_itr)
        r_square.append(r_square_itr)
        val_r_square.append(val_r_square_itr)    
        ReE_tf_itr = model(kin_array)
        ReE_tf.append(ReE_tf_itr) 
        # print("kinematics: ", kin_array)
        # print("ReE: ", ReE_tf)

    return np.array(loss, dtype=object),np.array(val_loss, dtype=object),np.array(r_square, dtype=object),np.array(val_r_square, dtype=object),np.array(ReE_tf)

def models_analysis(kinarray, iterations):
    loss = [] 
    val_loss = [] 
    r_square = [] 
    val_r_square = []
    ReE_tf = []
    # ReE_tf_mean = []
    # ReE_tf_min = []
    # ReE_tf_max = []
    # ReE_tf_std = []

    loss_m, val_loss_m, r_square_m, val_r_square_m, ReE_tf_m = load_iterations(kinarray, iterations)
    loss.append(loss_m)
    val_loss.append(val_loss_m)
    r_square.append(r_square_m)
    val_r_square.append(val_r_square_m)
    ReE_tf.append(ReE_tf_m)

    ReE_tf_mean, ReE_tf_min, ReE_tf_max, ReE_tf_std = band_limits(ReE_tf_m)
    # ReE_tf_mean_m, ReE_tf_min_m, ReE_tf_max_m, ReE_tf_std_m = band_limits(ReE_tf_m)
    # ReE_tf_mean.append(ReE_tf_mean_m)
    # ReE_tf_min.append(ReE_tf_min_m)
    # ReE_tf_max.append(ReE_tf_max_m)
    # ReE_tf_std.append(ReE_tf_std_m)
    
    models_params = {'loss': loss, 'val_loss': val_loss, 'r_square': r_square, 'val_r_square': val_r_square}
    models_ReE = {'ReE_tf': ReE_tf, 'ReE_tf_mean': ReE_tf_mean, 'ReE_tf_min': ReE_tf_min, 'ReE_tf_max': ReE_tf_max, 'ReE_tf_std': ReE_tf_std}

    return models_params, models_ReE 

# plot first 5 sets as a function of t
df_xmi = ROOT.RDF.FromCSV("./replicas/results_HallA_sets_21to25_data_5pct.csv")
dfcols_xmi = df_xmi.AsNumpy()

# plt.scatter(dfcols['t'], dfcols['xmi_ReE'])
# plt.errorbar(dfcols_xmi['t'], dfcols_xmi['xmi_ReE'], yerr=dfcols_xmi['exmi_ReE'], fmt='o')
# plt.show()

# Input kinematics for predictions, also contains km15 model value of ReE
# df = ROOT.RDataFrame("genkin_36_1", "./genkin_36_1.root") # kinematic sets and KM15 cffs 
# dfcols = df.AsNumpy()

# predict model n = sys.argv[1] times (number of iterations)
iterations = int(sys.argv[1])

# # Available models
# # models = [1, 2, 3, 4, 5, 6, 7, 8]
# models = [6, 7]
# N = len(models)
kin_array = np.arange(-0.6, -0.2, 0.001)

# Get ReE for all the models listed in 'models'
models_params, models_ReE = models_analysis(kin_array, iterations)

# Make a dataframe with all the values
df_dict = {"t":kin_array}
mean_dict = {"ReE_tf_mean": models_ReE['ReE_tf_mean']}
min_dict = {"ReE_tf_min": models_ReE['ReE_tf_min']}
max_dict = {"ReE_tf_max": models_ReE['ReE_tf_max']}
std_dict = {"ReE_tf_std": models_ReE['ReE_tf_std']}
# delta_dict = {"deltaReE": dfcols['ReE_km15'] - list(np.concatenate(models_ReE['ReE_tf_mean'][c])) for c in range(N)}

df_dict.update(mean_dict)
df_dict.update(min_dict)
df_dict.update(max_dict)
df_dict.update(std_dict)

df2 = ROOT.RDF.FromNumpy(df_dict)
df2.Display().Print()
df_np =df2.AsNumpy()


# print('t:', df_np['t'])
# print('ReE_tf_mean:', mean_dict['ReE_tf_mean'])
# t-dep analysis plots
# plt.plot(df_np['t'], df_np['ReE_km15'], color ="red", alpha=0.5)
plt.errorbar(dfcols_xmi['t'], dfcols_xmi['xmi_ReE'], yerr=dfcols_xmi['exmi_ReE'], fmt='o', color='red', label ='$\chi$MI')
plt.plot(df_dict['t'], df_dict['ReE_tf_mean'], color ="teal", alpha=0.5, label ='DNN model')
plt.fill_between(df_dict['t'], df_np['ReE_tf_mean']-df_np['ReE_tf_std'], df_np['ReE_tf_mean']+df_np['ReE_tf_std'], color ="teal", alpha=0.3)
plt.xlabel("$t[GeV^{2}]$", fontsize = 18)
plt.ylabel("$\mathfrak{Re}\mathcal{E}$", fontsize = 18)
plt.legend(fontsize = 18, loc ="upper left")
# plt.set_title('ReE vs t xB['+ str(xB_bin)+']')

plt.savefig("tdep2_vlosscut_std.jpg", dpi = 300)
plt.show()
# xB_bins = [0.36, 0.39, 0.48, 0.61]
# fig, ax = plt.subplots(2, len(xB_bins),figsize=(25.20,15.80))
# for idx_xB, xB_bin in enumerate(xB_bins): 
#     # For ReH vs t
#     ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_km15'][df_np['xB'] == xB_bin], color ="red", alpha=0.5)
#     # ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_1'][df_np['xB'] == xB_bin], color ="teal", alpha=0.5)
#     # ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_2'][df_np['xB'] == xB_bin], color ="blue", alpha=0.5)
#     # ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_3'][df_np['xB'] == xB_bin], color ="green", alpha=0.5)
#     # ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_4'][df_np['xB'] == xB_bin], color ="blueviolet", alpha=0.5)
#     # ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_5'][df_np['xB'] == xB_bin], color ="deeppink", alpha=0.5)
#     ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_6'][df_np['xB'] == xB_bin], color ="skyblue", alpha=0.5)
#     ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_7'][df_np['xB'] == xB_bin], color ="gold", alpha=0.5)
#     # ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_8'][df_np['xB'] == xB_bin], color ="darkorange", alpha=0.5)
#     # ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_1'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_1'][df_np['xB'] == xB_bin], color ="teal", alpha=0.3)
#     # ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_2'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_2'][df_np['xB'] == xB_bin], color ="blue", alpha=0.3)
#     # ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_3'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_3'][df_np['xB'] == xB_bin], color ="green", alpha=0.3)
#     # ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_4'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_4'][df_np['xB'] == xB_bin], color ="blueviolet", alpha=0.3)
#     # ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_5'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_5'][df_np['xB'] == xB_bin], color ="deeppink", alpha=0.3)
#     ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_6'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_6'][df_np['xB'] == xB_bin], color ="skyblue", alpha=0.3)
#     ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_7'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_7'][df_np['xB'] == xB_bin], color ="gold", alpha=0.3)
#     # ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_8'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_8'][df_np['xB'] == xB_bin], color ="darkorange", alpha=0.3)
#     ax[1, idx_xB].set_xlabel("$t[GeV^{2}]$")
#     ax[1, idx_xB].set_ylabel("ReH")
#     ax[1, idx_xB].legend(["KM15" , "Keras"])
#     ax[1, idx_xB].set_title('ReH vs t xB['+ str(xB_bin)+']')

# plt.savefig("tdep67.jpg", dpi = 300)



# # Draw model loss and accuracy (coeff. of determination) for model 1 iteration 1
# x_epoch = np.arange(0, len(models_params['loss'][0][0]), 1)

# fig, axis = plt.subplots(2, 2)
# axis[0, 0].plot(x_epoch,models_params['loss'][0][0],label='loss_mean')
# axis[0, 0].legend(['model_loss'])
# axis[0, 0].set_title("model 1, iteration 1")
# axis[0, 0].set_xlabel("epoch")
# axis[0, 0].set_ylabel("loss")

# axis[0, 1].plot(x_epoch,models_params['val_loss'][0][0],label='val_loss_mean')
# axis[0, 1].legend(['val_loss'])
# axis[0, 1].set_xlabel("epoch")
# axis[0, 1].set_ylabel("loss")

# axis[1, 0].plot(x_epoch,models_params['r_square'][0][0],label='r_square_mean')
# axis[1, 0].legend(['model_$R^2$'])
# axis[1, 0].set_xlabel("epoch")
# axis[1, 0].set_ylabel('$R^2$')

# axis[1, 1].plot(x_epoch,models_params['val_r_square'][0][0],label='val_r_square_mean')
# axis[1, 1].legend(['val_$R^2$'])
# axis[1, 1].set_xlabel("epoch")
# axis[1, 1].set_ylabel('$R^2$')

# # # Draw ReH distribution at some bins for model 1
# # fig, ax = plt.subplots(2, 4,figsize=(20.20,10.80))
# # ax = ax.flatten()
# # for i in range(8): 
# #     ax[i].hist(models_ReH['ReH_tf'][0][:,i*220]) #model 1
# #     ax[i].axvline(x = dfcols['ReH_km15'][i*220], color = 'r', label = 'km15')
# #     ax[i].axvline(x = df_np['ReH_tf_mean_1'][i*220], color = 'b', label = 'model1')
# #     ax[i].legend(['KM15', 'model1'])
# #     ax[i].set_title('model 1, xB = '+ str(df_np['xB'][i*220])+', t = '+str(df_np['t'][i*220]))
# #     ax[i].set_xlabel("ReH")
# #     ax[i].set_ylabel("dN/dReH")

# # t-dep analysis plots
# xB_bins = [0.36, 0.39, 0.48, 0.61]
# fig, ax = plt.subplots(2, len(xB_bins),figsize=(25.20,15.80))
# for idx_xB, xB_bin in enumerate(xB_bins): 
#     # ReH_km15 - ReH_mean_model1
#     # ax[0, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['deltaReH_1'][df_np['xB'] == xB_bin], color ="teal", alpha=0.5)
#     # ax[0, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['deltaReH_2'][df_np['xB'] == xB_bin], color ="blue", alpha=0.5)
#     # ax[0, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['deltaReH_3'][df_np['xB'] == xB_bin], color ="green", alpha=0.5)
#     # ax[0, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['deltaReH_4'][df_np['xB'] == xB_bin], color ="blueviolet", alpha=0.5)
#     # ax[0, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['deltaReH_5'][df_np['xB'] == xB_bin], color ="deeppink", alpha=0.5)
#     ax[0, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['deltaReH_6'][df_np['xB'] == xB_bin], color ="skyblue", alpha=0.5)
#     ax[0, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['deltaReH_7'][df_np['xB'] == xB_bin], color ="gold", alpha=0.5)
#     # ax[0, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['deltaReH_8'][df_np['xB'] == xB_bin], color ="darkorange", alpha=0.5)
#     ax[0, idx_xB].axhline(y = 0, color = 'black', linestyle = '-', alpha=0.5)
#     ax[0, idx_xB].set_xlabel("$t[GeV^{2}]$")
#     ax[0, idx_xB].set_ylabel("$\Delta ReH$")
#     ax[0, idx_xB].legend(["model1"])
#     ax[0, idx_xB].set_title('xB = '+ str(xB_bin))
#     # For ReH vs t
#     ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_km15'][df_np['xB'] == xB_bin], color ="red", alpha=0.5)
#     # ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_1'][df_np['xB'] == xB_bin], color ="teal", alpha=0.5)
#     # ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_2'][df_np['xB'] == xB_bin], color ="blue", alpha=0.5)
#     # ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_3'][df_np['xB'] == xB_bin], color ="green", alpha=0.5)
#     # ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_4'][df_np['xB'] == xB_bin], color ="blueviolet", alpha=0.5)
#     # ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_5'][df_np['xB'] == xB_bin], color ="deeppink", alpha=0.5)
#     ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_6'][df_np['xB'] == xB_bin], color ="skyblue", alpha=0.5)
#     ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_7'][df_np['xB'] == xB_bin], color ="gold", alpha=0.5)
#     # ax[1, idx_xB].plot(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_mean_8'][df_np['xB'] == xB_bin], color ="darkorange", alpha=0.5)
#     # ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_1'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_1'][df_np['xB'] == xB_bin], color ="teal", alpha=0.3)
#     # ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_2'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_2'][df_np['xB'] == xB_bin], color ="blue", alpha=0.3)
#     # ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_3'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_3'][df_np['xB'] == xB_bin], color ="green", alpha=0.3)
#     # ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_4'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_4'][df_np['xB'] == xB_bin], color ="blueviolet", alpha=0.3)
#     # ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_5'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_5'][df_np['xB'] == xB_bin], color ="deeppink", alpha=0.3)
#     ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_6'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_6'][df_np['xB'] == xB_bin], color ="skyblue", alpha=0.3)
#     ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_7'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_7'][df_np['xB'] == xB_bin], color ="gold", alpha=0.3)
#     # ax[1, idx_xB].fill_between(df_np['t'][df_np['xB'] == xB_bin], df_np['ReH_tf_min_8'][df_np['xB'] == xB_bin], df_np['ReH_tf_max_8'][df_np['xB'] == xB_bin], color ="darkorange", alpha=0.3)
#     ax[1, idx_xB].set_xlabel("$t[GeV^{2}]$")
#     ax[1, idx_xB].set_ylabel("ReH")
#     ax[1, idx_xB].legend(["KM15" , "Keras"])
#     ax[1, idx_xB].set_title('ReH vs t xB['+ str(xB_bin)+']')

# plt.savefig("tdep67.jpg", dpi = 300)

# # xB-dep analysis plots
# t_bins = [df_np['t'][1], df_np['t'][60], df_np['t'][150], df_np['t'][240]]
# fig, ax = plt.subplots(2, len(t_bins), figsize=(25.20,15.80))
# for idx_t, t_bin in enumerate(t_bins): 
#     # ReH_km15 - ReH_mean_model1
#     # ax[0, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['deltaReH_1'][df_np['t'] == t_bin], color ="teal", alpha=0.5)
#     # ax[0, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['deltaReH_2'][df_np['t'] == t_bin], color ="blue", alpha=0.5)
#     #ax[0, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['deltaReH_3'][df_np['t'] == t_bin], color ="green", alpha=0.5)
#     #ax[0, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['deltaReH_4'][df_np['t'] == t_bin], color ="blueviolet", alpha=0.5)
#     # ax[0, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['deltaReH_5'][df_np['t'] == t_bin], color ="deeppink", alpha=0.5)
#     ax[0, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['deltaReH_6'][df_np['t'] == t_bin], color ="skyblue", alpha=0.5)
#     ax[0, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['deltaReH_7'][df_np['t'] == t_bin], color ="gold", alpha=0.5)
#     # ax[0, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['deltaReH_8'][df_np['t'] == t_bin], color ="darkorange", alpha=0.5)
#     ax[0, idx_t].axhline(y = 0, color = 'black', linestyle = '-', alpha=0.5)
#     ax[0, idx_t].set_xlabel("$x_B$")
#     ax[0, idx_t].set_ylabel("$\Delta ReH$")
#     ax[0, idx_t].legend(["model1"])
#     ax[0, idx_t].set_title('t = '+ str(t_bin))
#     # For ReH vs xB
#     ax[1, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_km15'][df_np['t'] == t_bin], color ="red", alpha=0.5)
#     # ax[1, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_mean_1'][df_np['t'] == t_bin], color ="teal", alpha=0.5)
#     # ax[1, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_mean_2'][df_np['t'] == t_bin], color ="blue", alpha=0.5)
#     # ax[1, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_mean_3'][df_np['t'] == t_bin], color ="green", alpha=0.5)
#     # ax[1, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_mean_4'][df_np['t'] == t_bin], color ="blueviolet", alpha=0.5)
#     # ax[1, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_mean_5'][df_np['t'] == t_bin], color ="deeppink", alpha=0.5)
#     ax[1, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_mean_6'][df_np['t'] == t_bin], color ="skyblue", alpha=0.5)
#     ax[1, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_mean_7'][df_np['t'] == t_bin], color ="gold", alpha=0.5)
#     # ax[1, idx_t].plot(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_mean_8'][df_np['t'] == t_bin], color ="darkorange", alpha=0.5)
#     # ax[1, idx_t].fill_between(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_min_1'][df_np['t'] == t_bin], df_np['ReH_tf_max_1'][df_np['t'] == t_bin], color ="teal", alpha=0.3)
#     # ax[1, idx_t].fill_between(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_min_2'][df_np['t'] == t_bin], df_np['ReH_tf_max_2'][df_np['t'] == t_bin], color ="blue", alpha=0.3)
#     # ax[1, idx_t].fill_between(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_min_3'][df_np['t'] == t_bin], df_np['ReH_tf_max_3'][df_np['t'] == t_bin], color ="green", alpha=0.3)
#     # ax[1, idx_t].fill_between(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_min_4'][df_np['t'] == t_bin], df_np['ReH_tf_max_4'][df_np['t'] == t_bin], color ="blueviolet", alpha=0.3)
#     # ax[1, idx_t].fill_between(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_min_5'][df_np['t'] == t_bin], df_np['ReH_tf_max_5'][df_np['t'] == t_bin], color ="deeppink", alpha=0.3)
#     ax[1, idx_t].fill_between(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_min_6'][df_np['t'] == t_bin], df_np['ReH_tf_max_6'][df_np['t'] == t_bin], color ="skyblue", alpha=0.3)
#     ax[1, idx_t].fill_between(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_min_7'][df_np['t'] == t_bin], df_np['ReH_tf_max_7'][df_np['t'] == t_bin], color ="gold", alpha=0.3)
#     # ax[1, idx_t].fill_between(df_np['xB'][df_np['t'] == t_bin], df_np['ReH_tf_min_8'][df_np['t'] == t_bin], df_np['ReH_tf_max_8'][df_np['t'] == t_bin], color ="darkorange", alpha=0.3)
#     ax[1, idx_t].set_xlabel("$x_B$")
#     ax[1, idx_t].set_ylabel("ReH")
#     ax[1, idx_t].legend(["KM15" , "Keras"])
#     ax[1, idx_t].set_title('t['+ str(t_bin)+']')

# plt.savefig("xBdep67.jpg", dpi = 300)

# plt.show()


# tf.keras.backend.clear_session()
