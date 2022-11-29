#!/usr/bin/env python3# -*- coding: utf-8 -*-"""Created on Mon Aug  8 16:12:38 2022@author: karan94"""import numpy as npimport matplotlib.pyplot as pltplt.style.use('PhD_IMFT')# %% Read the dataellp_part1 = np.loadtxt('/Users/karan94/Documents/Ellipsoidal/Data_EllBoite/test_AR_1/ell_part.stat')ellp_part2 = np.loadtxt('/Users/karan94/Documents/Ellipsoidal/Data_EllBoite/test_AR_2_0/ell_part.stat')ellp_part1_5 = np.loadtxt('/Users/karan94/Documents/Ellipsoidal/Data_EllBoite/test_AR_1_5/ell_part.stat')ellp_part2_5 = np.loadtxt('/Users/karan94/Documents/Ellipsoidal/Data_EllBoite/test_AR_2_5/ell_part.stat')ellp_part3 = np.loadtxt('/Users/karan94/Documents/Ellipsoidal/Data_EllBoite/test_AR_3_0/ell_part.stat')# %% Absolute Statisticsfig = plt.figure(1, figsize=(12,12))plt.plot(ellp_part1[:,0], ellp_part1[:,19], 'k', linestyle='-', label=r'$\lambda$=1')plt.plot(ellp_part1[:,0], ellp_part1[:,20], 'k', linestyle='-')plt.plot(ellp_part1[:,0], ellp_part1[:,19]+ellp_part1[:,20], 'k', linestyle='--')plt.plot(ellp_part1_5[:,0], ellp_part1_5[:,19], 'xkcd:grey blue', linestyle='-', label=r'$\lambda$=1.5')plt.plot(ellp_part1_5[:,0], ellp_part1_5[:,20], 'xkcd:grey blue', linestyle='-')plt.plot(ellp_part1_5[:,0], ellp_part1_5[:,19]+ellp_part1_5[:,20], 'k', linestyle='--')plt.plot(ellp_part2[:,0], ellp_part2[:,19], 'xkcd:cement', linestyle='-', label=r'$\lambda$=2')plt.plot(ellp_part2[:,0], ellp_part2[:,20], 'xkcd:cement', linestyle='-')plt.plot(ellp_part2[:,0], ellp_part2[:,19]+ellp_part2[:,20], 'k', linestyle='--')plt.plot(ellp_part2_5[:,0], ellp_part2_5[:,19], 'xkcd:slate', linestyle='-', label=r'$\lambda$=2.5')plt.plot(ellp_part2_5[:,0], ellp_part2_5[:,20], 'xkcd:slate', linestyle='-')plt.plot(ellp_part2_5[:,0], ellp_part2_5[:,19]+ellp_part2_5[:,20], 'k', linestyle='--')plt.plot(ellp_part3[:,0], ellp_part3[:,19], 'xkcd:cobalt blue', linestyle='-', label=r'$\lambda$=3')plt.plot(ellp_part3[:,0], ellp_part3[:,20], 'xkcd:cobalt blue', linestyle='-')plt.plot(ellp_part3[:,0], ellp_part3[:,19]+ellp_part3[:,20], 'k', linestyle='--', label='Total')plt.xlabel('time (s)')plt.ylabel(r'K.E. $m^2s^{-2}$')# plt.ylim([0.0, 0.125])plt.legend(loc='lower right',ncol=2)plt.grid()plt.title(r'$N_p$=20000, $L_{(x=y=z)}$=0.314')plt.savefig('Absolute_KE_vs_time.pdf',dpi=500)plt.show()plt.close()# %% # fig = plt.figure(2, figsize=(12,12))# plt.plot(ellp_part2[:,0], ellp_part2[:,1], 'k', linestyle='-', label=r'$u_p$')# plt.plot(ellp_part2[:,0], ellp_part2[:,2], 'k', linestyle=':', label=r'$v_p$')# plt.plot(ellp_part2[:,0], ellp_part2[:,3], 'k', linestyle='--',label=r'$w_p$')# plt.xlabel('time (s)')# plt.ylabel(r'Velocity $ms^{-1}$')# #plt.ylim([0.0, 0.125])# plt.legend(loc='lower right')# plt.grid()# plt.title(r'$N_p$=20000, $L_{(x=y=z)}$=0.314')# plt.savefig('Absolute_TransVel_vs_time.pdf',dpi=500)# plt.show()# plt.close()# fig = plt.figure(3, figsize=(12,12))# plt.plot(ellp_part2[:,0], ellp_part2[:,4], 'k', linestyle='-', label=r'${\omega}_{p,x}$')# plt.plot(ellp_part2[:,0], ellp_part2[:,5], 'k', linestyle=':', label=r'${\omega}_{p,y}$')# plt.plot(ellp_part2[:,0], ellp_part2[:,6], 'k', linestyle='--',label=r'${\omega}_{p,z}$')# plt.xlabel('time (s)')# plt.ylabel(r' Angular Velocity $s^{-1}$')# #plt.ylim([0.0, 0.125])# plt.legend(loc='lower right')# plt.grid()# plt.title(r'$N_p$=20000, $L_{(x=y=z)}$=0.314')# plt.savefig('Absolute_AngVel_vs_time.pdf',dpi=500)# plt.show()# plt.close()# %% ellp_part2_tran = np.loadtxt('/Users/karan94/Documents/Ellipsoidal/Data_EllBoite/test_AR_2_0/ell_part_translation.stat')ellp_part2_rotation = np.loadtxt('/Users/karan94/Documents/Ellipsoidal/Data_EllBoite/test_AR_2_0/ell_part_angular.stat')a = 500.0e-06; Lambda = 2.0b = a/Lambda; c = bI_xx = (b*b + c*c)/5.0I_yy = (a*a + c*c)/5.0I_zz = (a*a + b*b)/5.0fig = plt.figure(9, figsize=(12,12))plt.plot(ellp_part2_tran[:,0], ellp_part2_tran[:,4], 'k', label=r'<$u_p u_p$>')plt.plot(ellp_part2_tran[:,0], ellp_part2_tran[:,5], 'xkcd:grey blue', label=r'<$v_p v_p$>')plt.plot(ellp_part2_tran[:,0], ellp_part2_tran[:,6], 'xkcd:cement', label=r'<$w_p w_p$>')plt.plot(ellp_part2_rotation[:,0], I_xx*ellp_part2_rotation[:,4], 'k', label=r'<${\omega}_{p,x} {\omega}_{p,x}$>')plt.plot(ellp_part2_rotation[:,0], I_yy*ellp_part2_rotation[:,5], 'xkcd:grey blue', label=r'<${\omega}_{p,y} {\omega}_{p,y}$>')plt.plot(ellp_part2_rotation[:,0], I_zz*ellp_part2_rotation[:,6], 'xkcd:cement', label=r'<${\omega}_{p,z} {\omega}_{p,z}$>')plt.xlabel('time (s)')plt.ylabel(r'$m^{2}s^{-1}$')#plt.ylim([0.0, 0.125])plt.legend(loc='lower right', ncol=2)plt.grid()plt.title(r'$N_p$=20000, $L_{(x=y=z)}$=0.314')#plt.savefig('Fluct_Trans_Stress_vs_time.pdf',dpi=500)plt.show()plt.close()# qp_tran = (ellp_part2_tran[:,4]+ellp_part2_tran[:,5]+ellp_part2_tran[:,6])/3.0# fig = plt.figure(10, figsize=(12,12))# plt.plot(ellp_part2_rotation[:,0], I_xx*ellp_part2_rotation[:,4], 'k', label=r'<${\omega}_{p,x} {\omega}_{p,x}$>')# plt.plot(ellp_part2_rotation[:,0], I_yy*ellp_part2_rotation[:,5], 'xkcd:grey blue', label=r'<${\omega}_{p,y} {\omega}_{p,y}$>')# plt.plot(ellp_part2_rotation[:,0], I_zz*ellp_part2_rotation[:,6], 'xkcd:cement', label=r'<${\omega}_{p,z} {\omega}_{p,z}$>')# plt.xlabel('time (s)')# plt.ylabel(r'$m^{2}s^{-1}$')# #plt.ylim([0.0, 0.125])# plt.legend(loc='upper right', ncol=3)# plt.grid()# plt.title(r'$N_p$=20000, $L_{(x=y=z)}$=0.314')# #plt.savefig('Fluct_Ang_Stress_vs_time.pdf',dpi=500)# plt.show()# plt.close()# qp_rot = (I_xx*ellp_part2_rotation[:,4]+ I_yy*ellp_part2_rotation[:,5]+I_zz*ellp_part2_rotation[:,6])/3.0# #%%# qp_tran = (ellp_part2_tran[:,4]+ellp_part2_tran[:,5]+ellp_part2_tran[:,6])# qp_rot = (I_xx*ellp_part2_rotation[:,4]+ I_yy*ellp_part2_rotation[:,5]+I_zz*ellp_part2_rotation[:,6])# fig = plt.figure(11, figsize=(12,12))# plt.plot(ellp_part2_tran[:,0], qp_tran, 'k', label=r'$q_p$-Translation')# plt.plot(ellp_part2_rotation[:,0], qp_rot, 'r', label=r'$q_p$-Rotation')# plt.plot(ellp_part2_tran[:,0], qp_tran+qp_rot, 'k', label=r'Total q_p')# plt.xlabel('time (s)')# plt.ylabel(r'$m^{2}s^{-1}$')# #plt.ylim([0.0, 0.125])# plt.legend(loc='upper right')# plt.grid()# plt.title(r'$N_p$=50000, $L_{(x=y=z)}$=0.314')# plt.show()# plt.close()# # %%# fig = plt.figure(10, figsize=(12,12))# # plt.plot(ellp_part2_rotation[:,0], ellp_part2_rotation[:,4], 'k', label=r'<$I_{xx}{\omega}_{p,x} {\omega}_{p,x}$>')# # plt.plot(ellp_part2_rotation[:,0], ellp_part2_rotation[:,5], 'k', label=r'<$I_{yy}{\omega}_{p,y} {\omega}_{p,y}$>')# # plt.plot(ellp_part2_rotation[:,0], ellp_part2_rotation[:,6], 'k', label=r'<$I_{zz}{\omega}_{p,z} {\omega}_{p,z}$>')# plt.plot(ellp_part2_rotation[:,0], ellp_part2_rotation[:,7], label=r'<$I_{xy}{\omega}_{p,x} {\omega}_{p,y}$>')# plt.plot(ellp_part2_rotation[:,0], ellp_part2_rotation[:,8], label=r'<$I_{xz}{\omega}_{p,x} {\omega}_{p,z}$>')# plt.plot(ellp_part2_rotation[:,0], ellp_part2_rotation[:,9], label=r'<$I_{yz}{\omega}_{p,y} {\omega}_{p,z}$>')# plt.xlabel('time (s)')# plt.ylabel(r'$m^{2}s^{-1}$')# #plt.ylim([0.0, 0.125])# plt.legend(loc='upper right')# plt.grid()# plt.title(r'$N_p$=50000, $L_{(x=y=z)}$=0.314')# plt.savefig('Fluct_Ang_Stress_vs_time.pdf',dpi=500)# plt.show()# plt.close()