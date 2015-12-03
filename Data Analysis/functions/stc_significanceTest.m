function [ h0 ] = stc_significanceTest( rawStimuli ,CP , alpha, numIter, eig_struct )
%STC_significanceTest preforms significance test to STC as described in
%Schwartz, Pillow et al. 2006
% target - test if the importance of the significant kernels are greater
% than wht would expected by chance.
%         input: alpha - confidence interval (CI)
%                numIter - number of simulation for single CI estimation
%                rawStimuli - the stimulus matrix as described and calculated in compute_stc
%                eig_struct - a structure of eigVecors, eigValues, out
%                             projected vector and left indexes on eigVal
% struct('eigVec',e_vec,'eigVal',e_val,'opVec',sta,'eigInd',(1:numel(e_val)));
% TODO: check why the random eigenVal are not the same scale as the real
% e_val

% create Shift matrix for the CP shift
s_mat = eye(length(CP)); % shift mat
s_ind = randsample(length(s_mat),numIter,false);
[ eigMax_val, eigMin_val ]  = deal( zeros(numIter,1) );
h0 = true; % stops when h0 is false - case where no axis are significant

% out project vector in opVec
rawStimuli = rawStimuli - bsxfun(@times, eig_struct.opVec',...
    rawStimuli * eig_struct.opVec)./(norm(eig_struct.opVec).^2); % shift X out projected

% iterative Monte Carlo simulation of CP shift and STE extraction
for k = 1:numIter
    s_vec = s_mat(:, s_ind(k));
    s_direc = randi(2,1)-1; % shift direction
        if s_direc; s_vec = flipud(s_vec);end;
    var1 = conv(CP,s_vec,'full'); % full conv implementation
    s_CP = var1(length(CP):end);
    s_X_squeeze = rawStimuli( logical(s_CP),:);
    s_CP_squeeze = CP( logical(s_CP) );
    s_X = bsxfun(@times, s_X_squeeze ,s_CP_squeeze);
    s_X_opsta = s_X - bsxfun(@times, mean(s_X,1),...
            s_X * mean(s_X,1)')./(norm(mean(s_X,1)).^2); % shift X out projected
    [~, ~, s_eig_val, ~, var_prec] = pca(s_X_opsta);
    if ~any(s_X)
        continue;
    end
    [ eigMax_val(k), eigMin_val(k) ] = deal( max(s_eig_val), min(s_eig_val) );
    
    if mod(k,100) == 0
        disp(['** loop - finished iteration k - ',num2str(k)])       
%         figure(1);
        drawnow
        plot(s_eig_val','o','MarkerFaceColor',[0 0.447058826684952 0.74117648601532]);
        hold on;
        plot(eig_struct.eigVal','o','MarkerFaceColor',[0 0.7 0.8117648601532]);
        hold off;
        xlabel('Eigenvalue number');
        ylabel('Eigen Value');
        title('STC Analysis');
    end
end

% clean from nans and find the confidence interval
nanInd = isnan(eigMax_val);
eigMin_val(nanInd) = [];
eigMax_val(nanInd) = [];
ci = [findCI( eigMax_val, alpha), findCI( eigMin_val, alpha)];

% check real eig_val
e_check = [ eig_struct.eigVal(1) eig_struct.eigVal(end)];
EVIL = eig_struct.eigInd; % Eigen Values Indexes Left (in each nest we prune some of the eig)

    %if( ci(1) > e_check(1) && ci(1) > e_check(end) && ci(end) > e_check(1) || ci(end) > e_check(end)  )
    if( ci(1) > e_check(1) && ci(end) <  e_check(end) )
        disp('H0 is accepted - non of our axis are significant');
        h0 = false;
        
    elseif max( abs( [e_check(1) e_check(1)] - [ ci(1) ci(end)]) ) >... 
            max( abs( [e_check(end) e_check(end)] - [ ci(1) ci(end)]) ) % the max eig_vec is more distant from the CI
        disp(' ** nest - H0 is rejected - prunning from top...');
        eig_struct = struct('eigVec',eig_struct.eigVec,'eigVal',eig_struct.eigVal,...
            'opVec',eig_struct.eigVec( :, EVIL(1) ),...
            'eigInd', EVIL(2:end) ); % update eig_struct
        clearvars -except h0 eig_struct rawStimuli CP alpha numIter
        [ h0 ] = stc_significanceTest( rawStimuli ,CP , alpha, numIter, eig_struct );    
        
    else                                                                % the min eig_vec is more distant from the CI
        disp(' ** nest - H0 is rejected - prunning from bottom...');
        
        eig_struct = struct('eigVec',eig_struct.eigVec,'eigVal',eig_struct.eigVal,...
            'opVec',eig_struct.eigVec( :, EVIL(end) )...
            ,'eigInd', EVIL(1:end-1) ); % update eig_struct
        clearvars -except h0 eig_struct rawStimuli CP alpha numIter
        [ h0 ] = stc_significanceTest( rawStimuli ,CP , alpha, numIter, eig_struct );
        
    end
end
%% Helper Func

function [ out ] = shiftfun(a,b)
% full conv implementation
    var2 = conv(a,b,'full');
    out = var2(length(a):end);
end

function [ ci ] = findCI( samples, alpha)
    SEM = std(samples)/sqrt(length(samples));               % Standard Error
    ts = tinv([1 - alpha/2 alpha/2]  ,length(samples)-1);      % T-Score
    ci = mean(samples) + ts*SEM ;  
end

