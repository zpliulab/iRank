# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 19:35:50 2019

@author: sea_sunshine
"""
import os
#os.chdir("D:\\shx_bioinformatics\\TCGA_Project\\used_for_plot_ROC_and_static")
os.chdir("D:\\git_data\\procedure\\HCC_procedure\\information_integration\\ORIrms")
#os.chdir("D:\\shx_bioinformatics\\TCGA_Project\\Copy_of_reg_mi_s")
import numpy as np
from sklearn import metrics
from scipy import interp
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
#from matplotlib.backends.backend_pdf import PdfPages
import scipy.io as scio
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
 
dataFile = 'RMS_3_33_new.mat'
data = scio.loadmat(dataFile)
Normal =data['Normal_pr']
Normal.shape[1]
Disease =data['D']
d_used =Disease[:,1]#只有一个
##
##贴标签
lable=np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
mean_fpr = np.linspace(0, 1, 200)
tprs=[]
aucs=[]
j=0
#pdf = PdfPages('my_figure.pdf')
f=plt.figure()
for i in range(0,Normal.shape[1]):#
    #i=0
    N_used=Normal[:,i]
    score = np.hstack((d_used,N_used))
    fpr, tpr, thresholds = metrics.roc_curve(lable, score)
    tprs.append(interp(mean_fpr, fpr, tpr))
    #此处是关键
    tprs[-1][0] = 0.0
    roc_auc = auc(fpr, tpr)
    roc_auc
    aucs.append(roc_auc)
#    plt.plot(fpr, tpr, lw=1, alpha=0.3,
#             label='ROC fold %d (AUC = %0.2f)' % (j, roc_auc))
    j += 1
plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
         label='Chance', alpha=.8)
mean_tpr = np.mean(tprs, axis=0)#列求和
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
itp = interp1d(mean_fpr, mean_tpr, kind='linear')
window_size, poly_order = 99, 3
list_x_new = np.linspace(min(mean_fpr), max(mean_tpr), 1000)
yy_sg = savgol_filter(itp(list_x_new), window_size, poly_order)
plt.plot(list_x_new, yy_sg, 'b', label= r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
         lw=2, alpha=.8)
#plt.plot(mean_fpr, mean_tpr, color='b',
#        label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
#         lw=2.5, alpha=.8)
std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
itp1 = interp1d(mean_fpr, tprs_upper, kind='linear')
yy_sg1 = savgol_filter(itp1(list_x_new), window_size, poly_order)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
itp2 = interp1d(mean_fpr, tprs_lower, kind='linear')
yy_sg2 = savgol_filter(itp2(list_x_new), window_size, poly_order)
plt.fill_between(list_x_new, yy_sg2, yy_sg1, color='grey', alpha=.2,
                 label=r'$\pm$ 1 std. dev.')
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('RMS_new')
plt.legend(loc="lower right")
plt.show()
f.savefig("RMS_new.pdf", bbox_inches='tight')
#pdf.savefig()  
#plt.close()
#pdf.close()
'''
import scipy.interpolate as interpolate
t, c, k = interpolate.splrep(mean_fpr, mean_tpr, s=0, k=4)
list_x_new = np.linspace(min(mean_fpr), max(mean_tpr), 1000)
spline = interpolate.BSpline(t, c, k, extrapolate=False)
#list_y_smooth = spline(mean_fpr, mean_tpr, list_x_new)
plt.plot(mean_fpr, mean_tpr, 'bo', label='Original points')
plt.plot(list_x_new, spline(list_x_new), 'r', label='BSpline')
plt.grid()
plt.legend(loc='best')
plt.show()

from scipy import interpolate
f = interpolate.interp1d(mean_fpr, mean_tpr, kind="linear")
y_int = f(list_x_new)
plt.plot(list_x_new, y_int, 'r', label='BSpline')



import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
# interpolate + smooth
itp = interp1d(mean_fpr, mean_tpr, kind='linear')
window_size, poly_order = 99, 3
yy_sg = savgol_filter(itp(list_x_new), window_size, poly_order)
# or fit to a global function
#def func(x, A, B, x0, sigma):
#    return A+B*np.tanh((x-x0)/sigma)

#fit, _ = curve_fit(func,mean_fpr, mean_tpr)
#yy_fit = func(list_x_new, *fit)
fig, ax = plt.subplots(figsize=(7, 4))
#plt.plot(mean_fpr, mean_tpr, 'r.', label= 'Unsmoothed curve')
#plt.plot(list_x_new, yy_fit, 'b--', label=r"$f(x) = A + B \tanh\left(\frac{x-x_0}{\sigma}\right)$")
plt.plot(list_x_new, yy_sg, 'k', label= "Smoothed curve")
plt.legend(loc='best')
'''