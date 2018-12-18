close all
clear
clc

% Computes moments of the solution VDF for sake of comparing 
% with the Pearson-IV distribution

% #######  PHYSICAL PARAMETERS  #######
m = 2.17e-25;
q = 1.602e-19;
kB = 1.38e-23; % [J/K]

Tt = 500; % [K] (Neutrals) Temperature in tangential direction
Tr = 500; % [K] (Neutrals) Temperature in radius direction

% #######  LOAD THE SOLUTION FILE  #######
% Set parameters (find them on the heading)
Nx = 100;
Nv = 100;

file_ID = 80;
filename = sprintf('../output/file_%08d.dat', file_ID);
dd = load(filename);

xx = dd(:,1);
vv = dd(:,2);
ff = dd(:,3);

% Reshape the stuff
XX = repmat(xx(1:Nv:end), 1, Nv);
VV = repmat(vv(1:Nv)', Nx, 1);
FF = reshape(ff, Nv, Nx)';

% ######  Compute moments at each location  #######

% Commodity values
intExpTr = 1; % sqrt(2*pi*kB*Tr / m);
intExpTt = 1; % sqrt(2*pi*kB*Tt / m);

v_n    = [];
v_rho  = [];
v_rhou = [];
v_u    = [];
v_P    = [];
v_Pxx  = [];
v_Prr  = [];

v_mu2  = [];
v_mu3  = [];
v_mu4  = [];

x_vec = XX(:,1);
v_vec = VV(1,:);

for(i = 1:Nx)

  % Extract variables
  x_now = XX(i,1);
  f_vec = FF(i,:); 

  % (Mass) Densities
  n   = trapz(v_vec, f_vec); % [1/m3] number density
  rho = m*n; % [kg/m3]

  % Momentum
  rhou = m*trapz(v_vec, v_vec.*f_vec);

  % Average velocity
%  u = rhou./(rho + max(rho/10000, 1e-25));
  u = rhou./(rho);

  % Pressure tensor
  c_x = v_vec - u;
  P_xx = m*trapz(v_vec, c_x.*c_x.*f_vec) * intExpTr * intExpTt;
  P_xr = 0;
  P_rr = n*kB*Tr;
  P_tt = n*kB*Tt;

  P = (P_xx + P_rr + P_tt)/3.0;

  % Second, third and fourth order moments for the Pearson-IV distribution
  v_mu2(i) = trapz(v_vec, c_x.^2.*f_vec);
  v_mu3(i) = trapz(v_vec, c_x.^3.*f_vec);
  v_mu4(i) = trapz(v_vec, c_x.^4.*f_vec);

  v_PEARS_rho(i)   = m*trapz(v_vec, f_vec);
  v_PEARS_rhov(i)  = m*trapz(v_vec, v_vec.*f_vec);
  v_PEARS_rhoth(i) = m*trapz(v_vec, c_x.^2.*f_vec);
  v_PEARS_q(i)     = m*trapz(v_vec, c_x.^3.*f_vec);
  v_PEARS_delta(i) = m*trapz(v_vec, c_x.^4.*f_vec);

  % Save everything into vectors
  v_n(i)    = n;
  v_rho(i)  = rho;
  v_rhou(i) = rhou;
  v_u(i)    = u;
  v_P(i)    = P;
  v_Pxx(i)  = P_xx;
  v_Prr(i)  = P_rr;

end

% SOME PLOTS
v_Q = v_mu3./(v_mu2.^(3/2));
v_D = v_mu4./(v_mu2.^2);


figure
loglog(v_Q(2:end-1), v_D(2:end-1), 'o', 'linewidth', 2);

%figure
%for ii = 2:numel(v_Q)-1
%  hold on
%  loglog(v_Q(ii), v_D(ii), 'o', 'linewidth', 2);
%  pause(0.01)
%end


% MORE PLOTS
v_PEARS_th = v_PEARS_rhoth./v_PEARS_rho;

xxx = v_PEARS_q./(v_PEARS_rho.*v_PEARS_th.^(3.0/2.0));
yyy = v_PEARS_delta./(v_PEARS_rho.*v_PEARS_th.^2);

figure
loglog(xxx, yyy, 'o', 'linewidth', 2)

%  
%  figure
%  subplot(3,1,1)
%  plot(x_vec', v_n, 'b', 'linewidth', 2)
%  xlabel('position [m]')
%  ylabel('number density [1/m3]')
%  
%  subplot(3,1,2)
%  plot(x_vec', v_u, 'b', 'linewidth', 2)
%  ylabel('Velocity [m/s]')
%  
%  subplot(3,1,3)
%  plot(x_vec', v_P, 'b', 'linewidth', 2)
%  ylabel('Pressure [Pa]')
%  
%  % Plot pressure tensor
%  figure
%  semilogy(x_vec', v_Pxx, '-r', 'linewidth', 2)
%  hold on
%  semilogy(x_vec', v_Prr, '-g', 'linewidth', 2)
%  semilogy(x_vec', v_P, '-b', 'linewidth', 2)
%  legend('Pxx','Prr = Ptt','P')
%  ylim([1e-10, Inf])
%  
%  figure
%  plot(x_vec, v_Pxx, 'r', 'linewidth', 2)
%  hold on
%  plot(x_vec, v_Prr, 'g', 'linewidth', 2)
%  plot(x_vec, v_P, 'b', 'linewidth', 2)
%  xlabel('Position [m]')
%  ylabel('Pressures [Pa]')
%  legend('Pxx','Prr = Ptt','P')
%  
%  %%% Plot comparison with Paper data
%  dd_paper = load('../data/n_ions.dat');
%  figure
%  plot(x_vec', v_n/1e6, 'b', 'linewidth', 2)
%  hold on
%  plot(dd_paper(:,1), dd_paper(:,2), 'or', 'linewidth', 2)
%  xlabel('Position [m]')
%  ylabel('Number density [cm^{-3}]')
