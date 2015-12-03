%% test STC with DATA
% run first estAcorrSTA_Data_3_2
%% Load data

cpxproc = process;
indx = 3; % good indx
% indx = 2; % not good indx
stim = cpxproc.stimulus.Yuncorr;
dt = cpxproc.stimulus.dt;
CP = cpxproc.spiketimes{indx};
%% Pre-Process

% downsample CP
if( length(cpxproc.CP{indx}) > length(cpxproc.stimulus.Yuncorr) )
    t = linspace(0, cpxproc.stimulus.dt * (length(stim)-1),length(stim));
    CP = histc(cpxproc.spiketimes{indx}, t);
end
%% estimate STA and STC and measure significance 
warning off;
[ sta1 ] = compute_sta( stim, cpxproc.spiketimes{indx}, 35, 0 , dt);
[ sta, w_sta ] = compute_white_sta( stim, CP, 35 );
method = 'suppress sta';
[e_vec, est_sta, var_prec, e_val, X, rawStimuli] = compute_stc( stim, CP, 35, sta, method);
%% plot module

figure(2);
% plot(normalize( est_sta,1 ) );
plot(normalize( sta1 ),'or')
hold on;
% plot(normalize( sta1 ),'or')
plot(normalize( sta ),'ok')
plot(normalize( w_sta ),'--k')
hold off;
legend('sta1', 'sta', 'w sta','Location','northwest');
% legend('mu', 'sta1', 'sta', 'w sta','Location','northwest');

figure(3);
k = 34;
plot(normalize(sta,1),'-or')
hold on;
plot(normalize( e_vec(:,k),1),'--k')
hold off;
legend('sta', ['e vector(',num2str(k),')'],'Location','northwest');

figure(1);
plot(var_prec','o','MarkerFaceColor',[0 0.447058826684952 0.74117648601532]);
xlabel('Eigenvalue number');
ylabel('Explanied precentage of Variance');
title('STC Analysis');
%% test Joint probability 

fig = figure(5);
JointFilterResponseProb( sta, e_vec(:,end), 39, X ,fig );
%% Significance test
tic;
% numSimulation = 1e3;
numSimulation = 500;
alpha = .05;
eigen_struct = struct('eigVec',e_vec,'eigVal',e_val,'opVec',sta,'eigInd',(1:numel(e_val)));
stc_significanceTest( rawStimuli , CP , alpha, numSimulation, eigen_struct );
toc
%% results
% it seems like the stimulus ensemble subspace covariance matrix has a
% unique pattern where the last eigen vector always explaines less variance
% compared to all the other eigen vectors. This essentially means