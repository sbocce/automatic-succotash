close all
clear
clc

% Cooks up some smooth fields
dd = load('../data/E_field_boeuf_1998.dat');
xx = dd(:,1);
EE = dd(:,2);

% Plot the field
figure
plot(xx, EE, 'or', 'linewidth', 2)

% Cook up the field
EE(1) = 0.0;
EE(2) = 80;

xadd = linspace(0.04,  0.05, 8)';
Eadd = linspace(20375, 20000, 8)';
xx = [xx; xadd(2:end)];
EE = [EE; Eadd(2:end)];

% Mirror the dataset (for polynomial regression reasons)
EE2 = [-EE(end:-1:2); EE];
xx2 = [-xx(end:-1:2); xx];

% Plot the field
hold on
plot(xx2, EE2, 'xg', 'linewidth', 2)
grid on


% Polyfit
%N = 11;
N = 20;
coefs = polyfit(xx2, EE2, N);
xval = linspace(-0.05,0.05,1000);
Eval = polyval(coefs, xval);

hold on 
plot(xval, Eval, '-m', 'linewidth', 2)

% Extract what is at x > 0.02
%pos_start = find(xval > 0.0245);
pos_start = find(xval > 0.026);
pos_start = pos_start(1);
pos_end   = find(xval > 0.04);
pos_end   = pos_end(1);

x_part2 = xval(pos_start:pos_end);
E_part2 = Eval(pos_start:pos_end);

plot(x_part2, E_part2, '--k', 'linewidth', 2);

%% Now attach with a straight line..
%x_part3 = [linspace(0, x_part2(1), 300), x_part2(2:end)];
%E_part3 = [linspace(0, E_part2(1), 300), E_part2(2:end)];
%
%plot(x_part3, E_part3, '--c', 'linewidth',2)

%%% Now attach a parabola that reaches the same second derivative
%%dDATA = ( E_part2(2) - E_part2(1) ) / (x_part2(2) - x_part2(1));
%%xxx3 = linspace(0, x_part2(1), 300);
%%EEE3 = dDATA/(2*x_part2(1))*xxx3.^2;

%xxx3 = linspace(0, x_part2(1), 300);
%EEE3 = E_part2(1)/(x_part2(1)^2) * xxx3.^2;
%x_part3 = [xxx3, x_part2(2:end)];
%E_part3 = [EEE3, E_part2(2:end)];
%
%plot(x_part3, E_part3, '--c', 'linewidth',2)

%Attach a parabola that reaches the same final point

NN3 = 16.97
xxx3 = linspace(0, x_part2(1), 300);
EEE3 = E_part2(1)/(x_part2(1)^NN3) * xxx3.^NN3;
x_part3 = [xxx3, x_part2(2:end)];
E_part3 = [EEE3, E_part2(2:end)];

plot(x_part3, E_part3, '--c', 'linewidth',2)

plot(x_part2(2), E_part2(2), 'xk')


% PLOT the first and second derivative as a diagnostic tool

figure
plot(x_part3(2:end), diff(E_part3));

figure
plot(x_part3(2:end-1), diff(diff(E_part3)));

