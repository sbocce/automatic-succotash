close all
clear
clc


page_screen_output(0);

% Plots solution in time

files_list = dir('../output/mom_*');

figure
for ii = 1:2:numel(files_list)

  % Load file
  dd = load(['../output/',files_list(ii).name]);
  dd = dd(2:end-1, :); % Exclude ghost cells

  fprintf('Data from: %s\n', files_list(ii).name);

  t_now = dd(1,1);
  xx    = dd(:,2);
  n     = dd(:,3);
  u     = dd(:,4);
  P     = dd(:,5);
  q     = dd(:,6);
  r     = dd(:,7);

  % Print stuff
  subplot(2,3,1)
  plot(xx, n, 'linewidth', 2);
  hold off
  title(['t = ', num2str(t_now), ' s']);
  grid on
  xlabel('Position [m]')
  ylabel('Number density [1/m3]')

  subplot(2,3,2)
  plot(xx, u, 'linewidth', 2);
  grid on
  ylabel('Velocity [m/s]')

  subplot(2,3,3)
  plot(xx, P, 'linewidth', 2);
  grid on
  ylabel('P [Pa]')

  subplot(2,3,4)
  plot(xx, q, 'linewidth', 2);
  grid on
  ylabel('q [W/m2]')

  subplot(2,3,5)
  plot(xx, r, 'linewidth', 2);
  grid on
  ylabel('r [..]')

  pause(0.001)
end
