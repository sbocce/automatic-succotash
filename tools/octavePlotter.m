close all
clear
clc

% This script plots a VDF given as a series of X, Y and Z points.

% Set parameters (find them on the heading)
Nx = 122;
Nv = 60;

dd = load('./output/file_00000021.dat');

xx = dd(:,1);
vv = dd(:,2);
ff = dd(:,3);

% Reshape the stuff
XX = repmat(xx(1:Nv:end), 1, Nv);
VV = repmat(vv(1:Nv)', Nx, 1);
FF = reshape(ff, Nv, Nx)';

% PLOT data
figure
surf(XX, VV, FF)
xlabel('x [m]')
ylabel('v [m/s]')
zlabel('f [s/ m6]')
