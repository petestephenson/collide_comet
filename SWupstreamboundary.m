%--------------------------------------------------------------------------
% Solar Wind Weighting
%--------------------------------------------------------------------------
lef = [0.08, 0.59];
bot = ones(2,1) * 0.17;
wid = [0.39,0.39];
heig = ones(2,1) * 0.75;

m = 9.11e-31; %kg
kBT = 10 * 1.6e-19;
% kB = 1.38e-23;
vsw = 4e5; %ms^-1
v_range = linspace(-vsw * 20, vsw * 20, 1e5);

vdfx = (m/(2 * pi * kBT))^0.5 * exp(- m .* (v_range - vsw).^2 ./ ( 2 * kBT) );% m^-1 s (if nSW added)
vdfy = (m/(2 * pi * kBT))^0.5 * exp(- m .* (v_range).^2 ./ ( 2 * kBT) );% m^-1 s
vdfz = (m/(2 * pi * kBT))^0.5 * exp(- m .* (v_range).^2 ./ ( 2 * kBT) );% m^-1 s

i_pos = v_range >= 0;
xFlux1 = trapz(v_range(i_pos), v_range(i_pos) .* vdfx(i_pos));% ms^-1
xFlux2 = trapz(v_range, vdfy);% No units
xFlux3 = trapz(v_range, vdfz);% No units
xFlux = xFlux1 * xFlux2 * xFlux3; % x nSW -> m^-2s^-1

yFlux1 = trapz(v_range, vdfx);
yFlux2 = trapz(v_range(i_pos), v_range(i_pos) .* vdfy(i_pos));
yFlux3 = trapz(v_range, vdfz);
yFlux = yFlux1 * yFlux2 * yFlux3;
%--------------------------------------------------------------------------
% Theoretical Vals
%--------------------------------------------------------------------------
vTh = sqrt(2 * kBT / m);
alpha = vsw / vTh;

xFluxTheory = (vTh / (2 * sqrt(pi))) *...
    ( exp(- alpha^2) + sqrt(pi) * alpha * ...
    (1 + erf(alpha)));
yFluxTheory = sqrt(kBT / (2 * pi * m));


%--------------------------------------------------------------------------
% Theoretical Vals
%--------------------------------------------------------------------------
m1 = 9.11e-31;
m2 = 1.67e-29;
vdfx1 = (m1/(2 * pi * kBT))^0.5 * exp(- m1 .* (v_range - vsw).^2 ./ ( 2 * kBT) );% m^-1 s (if nSW added)
vdfx2 = (m2/(2 * pi * kBT))^0.5 * exp(- m2 .* (v_range - vsw).^2 ./ ( 2 * kBT) );% m^-1 s (if nSW added)

f = figure('units', 'centimeters','Position', [1,1,22,11]);

ax1 = axes(f, 'Position', [lef(1), bot(1), wid(1), heig(1)]);
hold(ax1, 'on');
xfill = [v_range(v_range>0), fliplr(v_range(v_range>0))];
yfill1 = [vdfx1(v_range>0), fliplr(zeros(1,sum(v_range>0)))];

fill(ax1, xfill, yfill1, 'b', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
yfill2 = [vdfx2(v_range>0), fliplr(zeros(1,sum(v_range>0)))];
fill(ax1, xfill, yfill2, 'g', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
plot(ax1, [0,0],[0,1.2 * max(vdfx2)], 'k', 'LineWidth', 1.3, 'HandleVisibility', 'off');
plot(ax1, v_range, vdfx1, 'DisplayName', 'm_e = 9.11e-31 kg', 'LineWidth', 1.3);
plot(ax1, v_range, vdfx2, 'DisplayName', 'm_e = m_p /100', 'LineWidth', 1.3);

xlabel(ax1,'X Velocity [ms^{-1}]')
legend(ax1);
xlim(ax1,[-vsw,Inf]);
ylabel(ax1,'vdf [m^{-2}s^1]');
ylim(ax1,[0,1.1 * max(vdfx2)]);
set(ax1, 'FontSize', 14);
% hold on

%--------------------------------------------------------------------------
% Electron Mass Tests
%--------------------------------------------------------------------------
% Looking at how the electron mass effects the energy distribution of pls
% entering the domain from the solar wind. 


m = 9.11e-31; %kg
kBT = 10 * 1.6e-19;
% kB = 1.38e-23;
vsw = 4e5; %ms^-1
v_range = linspace(-vsw * 20, vsw * 20, 1e4);

vdfx = (m/(2 * pi * kBT))^0.5 * exp(- m .* (v_range - vsw).^2 ./ ( 2 * kBT) );% m^-1 s (if nSW added)
vdfy = (m/(2 * pi * kBT))^0.5 * exp(- m .* (v_range).^2 ./ ( 2 * kBT) );% m^-1 s
vdfz = (m/(2 * pi * kBT))^0.5 * exp(- m .* (v_range).^2 ./ ( 2 * kBT) );% m^-1 s
xDens1 = trapz(v_range, vdfx);% No units
xDens2 = trapz(v_range, vdfy);% No units
xDens3 = trapz(v_range, vdfz);% No units
xDens = xDens1 * xDens2 * xDens3; % x nSW -> m^-2s^-1

%--------------------------------------------------------------------------
% Monte Carlo int in cartesian
%--------------------------------------------------------------------------
nSig = 30;
nSamples = 1000000;
v_rand = (rand(1000000,3) - 0.5) * vsw * nSig;
vdf_rand = (m/(2 * pi * kBT))^1.5 * exp(- m .* sum((v_rand - [vsw,0,0]).^2,2) ./ ( 2 * kBT) );% m^-1 s (if nSW added)
I1 = sum(vdf_rand) * (vsw * nSig)^3 / nSamples;

% vdfx = (m/(2 * pi * kBT))^0.5 * exp(- m .* (v_range - vsw).^2 ./ ( 2 * kBT) );% m^-1 s (if nSW added)
% vdfy = (m/(2 * pi * kBT))^0.5 * exp(- m .* (v_range).^2 ./ ( 2 * kBT) );% m^-1 s
% vdfz = (m/(2 * pi * kBT))^0.5 * exp(- m .* (v_range).^2 ./ ( 2 * kBT) );% m^-1 s
%--------------------------------------------------------------------------
% -------------- Conversion to Energy space ----------------------------
%--------------------------------------------------------------------------
m = 1.67e-29;
% m = 9.11e-31;
E_SW = 0.5 * m * vsw ^ 2;
Edf_const = ((m/(2 * pi * kBT)) ^ 1.5) * sqrt(2/m^3) * exp(- E_SW / kBT);

% Monte carlo integration
nSamples = 1000000;
E_max = 200 * 1.6e-19;
E_range = rand(nSamples,1) * E_max;
theta_rand = rand(nSamples,1) * pi;
phi_rand = (rand(nSamples,1) - 0.5) * pi;
Vol_box = E_max * pi ^2;

edf_rand = sin(theta_rand) .* sqrt(E_range) .* exp( - E_range ./ kBT) .* ...
    exp( (m * vsw / kBT) .* sqrt(2 * E_range ./ m) .* sin(theta_rand) .* cos(phi_rand));

I2 = Edf_const .* Vol_box .* sum(edf_rand) / nSamples;

%--------------------------------------------------------------------------
% Integrating over angles for each energy
%--------------------------------------------------------------------------
m = 1.67e-29;
% m = 9.11e-31;


% Monte carlo integration
nEng = 300;
nSamples = 100000;
E_max = 200 * 1.6e-19;
E_range = linspace(0,E_max, nEng);
theta_rand = rand(nSamples,1) * pi;
phi_rand = (rand(nSamples,1) - 0.5) * pi;
Vol_box = pi ^2;
Edf = zeros(2, nEng);
% Low electron
m = [9.11e-31,1.67e-29];
for iMass = 1:2
    E_SW = 0.5 * m(iMass) * vsw ^ 2;
    Edf_const = ((m(iMass)/(2 * pi * kBT)) ^ 1.5) * sqrt(2/m(iMass)^3) * exp(- E_SW / kBT);
    edf_rand = sin(theta_rand) .* sqrt(E_range) .* exp( - E_range ./ kBT) .* ...
        exp( (m(iMass) * vsw / kBT) .* sqrt(2 * E_range ./ m(iMass)) .* sin(theta_rand) .* cos(phi_rand));
    Edf(iMass,:) = Edf_const .* Vol_box .* sum(edf_rand,1) / nSamples;
end
m = 9.11e-31;


% f = figure;
ax2 = axes(f, 'Position', [lef(2), bot(2), wid(2), heig(2)]);
hold(ax2,  'on');
plot(ax2,E_range/1.6e-19, Edf(1,:)*1.6e-19, 'DisplayName', 'm_e = 9.11e-31 kg', 'LineWidth', 1.3);
plot(ax2,E_range/1.6e-19, Edf(2,:)*1.6e-19, 'DisplayName', 'm_e = m_p /100', 'LineWidth', 1.3);
ylabel(ax2, 'EDF [cm^{-3} eV^{-1}]');
xlabel(ax2, 'Energy [eV]');
legend(ax2);
set(ax2, 'FontSize', 14);
% title(ax2, 'Distribution of SW elecs entering at upstream boundary')