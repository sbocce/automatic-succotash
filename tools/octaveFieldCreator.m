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


% Plot the field
hold on
plot(xx, EE, 'ob', 'linewidth', 2)
grid on

% Augment the dataset (for polynomial regression reasons)
EE2 = [-EE(end:-1:2); EE];
xx2 = [-xx(end:-1:2); xx];

% Plot the field
hold on
plot(xx2, EE2, 'xg', 'linewidth', 2)
grid on

% Polyfit
%N = 11;
N = 13;
coefs = polyfit(xx2, EE2, N);
xval = linspace(-0.05,0.05,1000);
Eval = polyval(coefs, xval);

hold on 
plot(xval, Eval, '-m', 'linewidth', 2)

% Extract what is at x > 0.02
pos_start = find(xval > 0.02);
pos_start = pos_start(1);
pos_end   = find(xval > 0.038);
pos_end   = pos_end(1);

x_part2 = xval(pos_start:pos_end);
E_part2 = Eval(pos_start:pos_end);

plot(x_part2, E_part2, '--k', 'linewidth', 2);





