% test STC

%% generate a dummy process from
% incoparate 2 different linear filters.

% simulate process
Fs = 50;
amp = 0.8;
% nl = @(x) exp(0 + 1 * x);
nl = @(x) x.^2;
filtSupp = linspace( 0, 2*pi, 3/(1/Fs));
filterMatrix = [ cos(pi/8 + filtSupp) .* 2.4 .* exp(-(filtSupp- pi).^2)/4 ;...
                2* exp(-(filtSupp - pi).^2)/5 ]';
filterRelationfunc = @(y) y(:, 1) + 2 * y(:, 2);
procSettings = struct('Fs',Fs,...
    'size', 5*1e4,...
    'noiseMean',0,...
    'noiseSTD',2,...
    'nl',nl,...
    'filterMatrix',filterMatrix,...
    'filterRelationfunc',filterRelationfunc);
cpxproc = ComplexCellSimulation( procSettings );
cpxproc = generateProcess( cpxproc, amp );

%% estimate STA and STC and measure significance 
[ sta1 ] = compute_sta( cpxproc.stimulus, cpxproc.spiketimes, 3/(1/Fs), 0 , cpxproc.dt);

[ sta, w_sta ] = compute_white_sta(cpxproc.stimulus, cpxproc.CP, 3/(1/Fs) );

[e_vec, est_sta, var_prec] = compute_stc( cpxproc.stimulus, cpxproc.CP, 3/(1/Fs), sta, 'suppress sta');
%% plots

figure;
hold on;
plot(sta1)
plot(sta,'--r')
plot(w_sta,'ok')
legend('sta','sta2','white sta');
hold off;

figure;
hold on;
plot(sta)
plot(filterMatrix(:,1),'--r')
plot(filterMatrix(:,2),'ok')
hold off;

figure;
hold on;
plot(est_sta);
plot(filterMatrix(:,1),'or')
plot(filterMatrix(:,2),'ok')
plot(e_vec(:,1),'--k')
plot(e_vec(:,2),'--r')
hold off;

figure();
plot(var_prec','o','MarkerFaceColor',[0 0.447058826684952 0.74117648601532]);
xlabel('Eigenvalue number');
ylabel('Explanied precentage of Variance');
title('STC Analysis');