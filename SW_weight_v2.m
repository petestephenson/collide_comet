%% SW weighting Test

%% SW Properties
n_SW = 1; %cm^-3
T_SW = 10; %eV
m = 9.11e-31;
sig_v = sqrt(T_SW * 1.6e-19 / 9.11e-31);

v_sw = 400e3;
kBT = 1.6e-19 * T_SW;

%% 
nPls = 1000000;
r = rand(nPls, 4);
r_gaus = zeros(nPls,3);
r_gaus(:,1) = sqrt(-2 .* log(r(:,1))) .* cos( 2 .* pi .* r(:,2));
r_gaus(:,2) = sqrt(-2 .* log(r(:,1))) .* sin( 2 .* pi .* r(:,2));
r_gaus(:,3) = sqrt(-2 .* log(r(:,3))) .* cos( 2 .* pi .* r(:,4));

v_rand = sig_v * r_gaus;
v_rand(:,1) = v_rand(:,1) + v_sw;

f = figure;
v_grid = linspace(- 4 * sig_v, 4 * sig_v, 100);
h = cell(3,1);
ax = axes(f);
dim_name = {'x', 'y', 'z'};
hold(ax, 'on');
for i_dim = 1:3
    h{i_dim} = histogram(ax, v_rand(:, i_dim), v_grid, 'DisplayName', ['v_',dim_name{i_dim}, ' Code'] );
%     vdfscale(i_dim) = max(h{i_dim}.Values, [], 'all') / max(vdf(i_dim, :), [], 'all');
%     vdfplot = plot(v_range, vdf(i_dim, :) * vdfscale(i_dim), 'DisplayName', ['v_',dim_name{i_dim}, 'Theory']);
end
legend;
xlabel('Velocity [m s^{-1}]');
ylabel('Counts');
set(gca, 'FontSize', 14);
% title(ax, ['Vel Dist through ', fName{i_pl_face},' Face']);

%% Calculate weight of SW pls
%Through x bound:
i_dim = 2;
i_inc = v_rand(:,i_dim) > 0;
nPlsInc = sum(i_inc);

Weight = n_SW * v_rand(i_inc,i_dim) * 1e2/ nPlsInc; 
TotWeight = sum(Weight);
%%
nBins = 5e2;
v_range = linspace(-v_sw * 20, v_sw * 20, nBins+1);

vdfx = (m/(2 * pi * kBT))^0.5 * exp(- m .* (v_range - v_sw).^2 ./ ( 2 * kBT) );% m^-1 s (if nSW added)
vdfy = (m/(2 * pi * kBT))^0.5 * exp(- m .* (v_range).^2 ./ ( 2 * kBT) );% m^-1 s
vdfz = (m/(2 * pi * kBT))^0.5 * exp(- m .* (v_range).^2 ./ ( 2 * kBT) );% m^-1 s

%% Weight against velocity
ixBin = discretize(v_rand(i_inc, i_dim), v_range);
WeightBinned = zeros(nBins,1);
for iBin = 1:nBins
    WeightBinned(iBin) = sum(Weight(ixBin==iBin));
end

%%
dv = v_range(2:end) - v_range(1:end-1);
if i_dim == 1
vdfPlot = vdfx;
else
    vdfPlot = vdfy;
end
f = figure;
ax = axes(f);
histogram(ax, 'BinEdges',v_range, 'BinCounts', WeightBinned./dv');
hold on 
rat = max(WeightBinned./dv',[],'all') / max(v_range .* vdfPlot,[],'all');
plot(ax, v_range, v_range .* vdfPlot * rat);
xlabel(ax, 'Velocity [ms^{-1}]');
ylabel(ax, 'Weight through Bin');

%% Theoretical Values - Integral
vTh = sqrt(2 * kBT / m);
alpha = v_sw / vTh;

xFluxTheory = (vTh / (2 * sqrt(pi))) *...
    ( exp(- alpha^2) + sqrt(pi) * alpha * ...
    (1 + erf(alpha)));
yFluxTheory = sqrt(kBT / (2 * pi * m));

%% Peka
v_peak = zeros(3,1);
v_peak(1) = 0.5 * v_sw + sqrt((0.25 * v_sw^2) + (kBT/m));
v_peak(2:3) = sig_v;






