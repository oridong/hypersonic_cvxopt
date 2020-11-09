% Skye Mceowen
% Qualifying Exam Plots
% Nov4, 2020

clear all, close all, clc
% Select directory (later use if statement to pick quals vs wang)
addpath('matfiles/')

ev = vehicle;

h = figure;
ev.plot_bank_sweep(h)


%savefig('matfiles/BankAngle.fig')
%save('matfiles/bankangle_figs.mat')