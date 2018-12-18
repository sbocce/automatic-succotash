close all
clear
clc

page_screen_output(0);

% This script plots a VDF given as a series of X, Y and Z points.

% Set parameters (find them on the heading)
Nx = 100;
Nv = 100;

figure
for kk = 0:4:120
%for kk = 0:2

  fprintf('Step %d\n', kk);

  filename = sprintf('../output/file_%08d.dat', kk);
  dd = load(filename);
  
  xx = dd(:,1);
  vv = dd(:,2);
  ff = dd(:,3);
  
  % Reshape the stuff
  XX = repmat(xx(1:Nv:end), 1, Nv);
  VV = repmat(vv(1:Nv)', Nx, 1);
  FF = reshape(ff, Nv, Nx)';
  
  % PLOT data
  surf(XX, VV, FF)
  shading flat
  xlabel('x [m]')
  ylabel('v [m/s]')
  zlabel('f [s/m^6]')
%  zlim([0,1.2]);
%  view(70,60)
  view(140,40)

  pause(0.001)
  if(kk == 0)
    pause()
  end

end

% PLOT VDF LINE BY LINE
figure
plot3(VV', XX', FF','k')
xlabel('v [m/s]')
ylabel('x [m]')
zlabel('f [s^3/m^6]')
view(140,40)
