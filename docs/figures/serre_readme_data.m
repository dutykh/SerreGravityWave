% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% Data of the figure of the README: the head-on collision of two solitary waves computed by
% script/SerreGravityWave.m, run as it stands. The function returns nothing and draws three
% figures; this script reads back the free surface that its space-time plot holds (51 profiles,
% every 0.8 time units) and the energy that its last figure holds, checks the grid and the times
% it reads, and writes them to eta.csv (t, x, eta) and energy.csv (t, E(t) - E(0)).
%
% Run it from this folder with GNU Octave (the reading follows the surface drawn by the waterfall
% function of Octave):
%   octave --no-gui serre_readme_data.m
% plot_readme_figure.py then draws serre-collision.png from the two files.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'script'));
set(0, 'defaultfigurevisible', 'off');

SerreGravityWave();

% The space-time plot: the surface of the second figure, surrounded by the two rows and
% columns of the curtain that the waterfall function adds on each side.
s = findobj(2, 'type', 'surface');
X = get(s, 'xdata'); Y = get(s, 'ydata'); Z = get(s, 'zdata');
x = X(3:end-2); t = Y(3:end-2); eta = Z(3:end-2, 3:end-2);
N = 2048; l = 75.0;
if numel(x) ~= N || max(abs(x(:) - (1-N/2:N/2)'*2*l/N)) > 1e-12 || max(abs(t(:) - (0:0.8:40)')) > 1e-12
  error('serre_readme_data: the space-time plot does not hold the expected grid and times');
end
[T, Xg] = ndgrid(t(:), x(:));
dlmwrite(fullfile(here, 'eta.csv'), [T(:) Xg(:) eta(:)], 'precision', '%.16e');

% The energy: the third figure shows E(t) minus its mean at the 101 stored times.
e = findobj(3, 'type', 'line');
te = get(e, 'xdata'); E = get(e, 'ydata');
dlmwrite(fullfile(here, 'energy.csv'), [te(:) E(:) - E(1)], 'precision', '%.16e');

% The energy of the initial condition, with the formulas of InitCond and Hamilt of the
% function (g = d = 1, solitary waves of amplitudes 0.2 and 0.3 centred at x = -25 and 25).
xx = x(:); dx = 2*l/N;
a1 = 0.2; a2 = 0.3;
c1 = sqrt(1 + a1); k1 = sqrt(3*a1/(a1 + 1));
c2 = sqrt(1 + a2); k2 = sqrt(3*a2/(a2 + 1));
eta1 = a1*sech(0.5*k1*(xx + 25)).^2; eta2 = a2*sech(0.5*k2*(xx - 25)).^2;
h = 1 + eta1 + eta2; u = (c1*eta1 - c2*eta2)./h;
k = [0:N/2 1-N/2:-1]'*pi/l; k(N/4+2:N/4*3) = 0;
ux = real(ifft(1i*k.*fft(u)));
E0 = 0.5*dx*sum(h.*u.^2 + h.^3.*ux.^2/3 + (eta1 + eta2).^2);
[~, m] = max(max(eta, [], 2));
fprintf('E(0) = %.10f, largest |E(t) - E(0)| = %.3e at t = %.1f\n', E0, max(abs(E - E(1))), te(find(abs(E - E(1)) == max(abs(E - E(1))), 1)));
fprintf('largest elevation %.4f at t = %.1f\n', max(eta(:)), t(m));
