close all
clear
clc

% Cooks up some smooth fields
%dd = load('../data/E_field_boeuf_1998.dat');
dd = load('../data/S_ionization_boeuf_1998_data.dat');
xx = dd(:,1);
SS = dd(:,2);

% Plot the field
figure
plot(xx, SS, 'or', 'linewidth', 2)

% Cook up the field
xadd = linspace(0.04,  0.05, 8)';
Sadd = linspace(SS(end), 8*(SS(end)-SS(end-1)), 8)';
xx = [xx; xadd(2:end)];
SS = [SS; Sadd(2:end)];


% Mirror the dataset (for polynomial regression reasons)
SS2 = [-SS(end:-1:2); SS];
xx2 = [-xx(end:-1:2); xx];

% Plot the field
hold on
plot(xx2, SS2, 'xg', 'linewidth', 2)
grid on

% Polyfit
%N = 12;
N = 25;
coefs = polyfit(xx2, SS2, N);
xval = linspace(-0.05,0.05,1000);
Sval = polyval(coefs, xval);

hold on 
plot(xval, Sval, '-m', 'linewidth', 2)

% Extract what is at x > 0
pos_start = find(xval > 0.0);
pos_start = pos_start(1);
pos_end   = find(xval >= 0.04);
pos_end   = pos_end(1);

x_part2 = xval(pos_start:pos_end);
S_part2 = Sval(pos_start:pos_end);

plot(x_part2, S_part2, '--k', 'linewidth', 2);

% Take only what is positive  and attach some zeros
neg_pos = find(S_part2 < 0);
neg_pos = neg_pos(end);

x_part3 = x_part2(neg_pos+1:end);
S_part3 = S_part2(neg_pos+1:end);


plot(x_part3, S_part3, '--c', 'linewidth',2)

% PLOT the first and second derivative as a diagnostic tool

figure
plot(x_part3(2:end), diff(S_part3));

figure
plot(x_part3(2:end-1), diff(diff(S_part3)));

