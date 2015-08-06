%% Local Sensetivity Analysis for rho
% Analysis of Variance

% General simulation parameters
N = 5e4;

% parameter space and range
space_sig = linspace(1e-1, 1 , 200);      
space_rho = linspace(0, 2*pi, 200);         space_f = linspace(0, 8 ,200); 

% define models and cost
X = linspace(-2,2,500);
[sig_ITER, f_ITER] = deal( zeros(N,1) );
model_sparse = @(s,f) (exp((-X.^2)./(2*s.^2)).*cos((2*pi*f).*X))-exp((-X.^2-(2*pi*f)^2.*s.^4)/(2*s.^2)).*cos(2*(2*pi*f)*0);
model_full = @(s,f,rho) bsxfun(   @minus, (exp((-X.^2)/(2*s.^2)).*cos((2*pi*f).*X))',...
    bsxfun(@times , exp((-X.^2-(2*pi*f)^2*s.^4)./(2*s.^2))', cos(2*(2*pi*f).*rho))   );

% here we calc for random f and sig the abs differences of the extended
% model and the suggested model. 
try 
    parfor n = 1:N
        randIdx = randi(200,2,1);
        f_ITER(n) = space_f( randIdx(1) ); 
        sig_ITER(n) = space_sig( randIdx(2) );
        Y_sparse = model_sparse(sig_ITER(n),f_ITER(n));
        Y_full =  model_full(sig_ITER(n),f_ITER(n), space_rho);
        res_AE = abs( bsxfun(@minus, Y_full, Y_sparse') ); % residuals absolute error
        
        varOfDiff(:,n) = var(res_AE); % variance of the error for each possible rho
        totalVar(n) = var( varOfDiff(:,n) ); % variance of all sampled f, sig 
    end
catch exception 
    disp(exception);
end
    %% Plots

    figure(3); plot( (1:N), totalVar );
    xlabel('iter'); ylabel('Var [ Var(res_i) ]');
    title('Total Variance for changes in \rho, \sigma and  f'); 
    
    k = 114;
    figure(2); plot(space_rho, varOfDiff(:,k));
    xlabel('\rho'); ylabel('Var of residuals');
    set(gca,'XTick',0:pi/2:2*pi)
    set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    title('Var (|Y_{full} - Y_{sparse}|) '); 
    
    k_rho1 = 0; k_rho2 = pi;
    figure(1);plot(X, model_sparse(sig_ITER(k),f_ITER(k)));hold on; plot(X, model_full(sig_ITER(k),f_ITER(k), k_rho1),'--r');
    plot(X ,model_full(sig_ITER(k),f_ITER(k), k_rho2),'--k');hold off;
    legend('model sparse',['model full rho = ',num2str(k_rho1)],['model full rho = ',num2str(k_rho2)])
    title('Model evalutaion for different \rho');
    xlabel('time'); ylabel('Model Value');
    %% density analysis of faulty parameters
    
    faultyIter = totalVar >  0.001;
    faultyParam = [sig_ITER(faultyIter) f_ITER(faultyIter)];
    hist3(faultyParam,[50 50])
    set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
    xlabel('\sigma'); ylabel('f');