% $Author: Yasuko Matsubara 
% $Date: 2014-04-30

%% demo: fitting sequence "./sequence.dat"
fn='./six_type_data/d5.txt';
% duration of sequence
T=24*4;
% # of max iteration
ITER=20;
% daily periodicity (24hours)
pfreq=24;


dat=load(fn);
% dat=dat(1:T);
outfn='output';
wantPlot=1; % Yes!
disp('===================================');
disp('DEMO - fitting sample sequence');
disp('-----------------------------------');
disp(['- filename = ', fn]);
disp(['- duration = ', num2str(T)]);
disp(['- max iteration = ', num2str(ITER)]);
disp('===================================');
disp(' ');
[RSE, params]=M_spikeMfit(dat, pfreq, outfn, ITER, wantPlot);

