clear all; close all; clc;

Ar = load (['waist_31/sPPLT_N_0.5/As_temporal_r.dat']);
Ai = load (['waist_31/sPPLT_N_0.5/As_temporal_i.dat']);

A1 = Ar + 1i*Ai;

hold on
plot(abs(A1).^2)
% zlim   ([0,6])
% shading interp;

ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];

