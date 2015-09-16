%% stim
dt_stim = 0.01;
N = 1e5;
noiseVar = 1;
noiseMean = 0;
stim = noiseMean + noiseVar * randn(N,1);
t_stim = linspace(0,N*dt_stim,N);

%% lin k
space_sig = linspace(0.1, 0.15, 200);       space_mu = linspace(-0.6, -0.2, 200); 
space_fi = linspace(pi, 2*pi, 200);         space_f = linspace(8, 15 ,200);
p_kernel = [space_sig(randIdx(1)) space_mu(randIdx(2)) space_fi(randIdx(3)) space_f(randIdx(4))];
dt_kernel = 1/30;
amplitude = 0.1;
model_kernel = @(t) amplitude * exp(-(t-p_kernel(2)).^2/p_kernel(1)^2).*sin( p_kernel(4) *(t-p_kernel(3)) );
t_kernel = linspace(-30*dt_kernel, 0*dt_kernel, 1/dt_stim); 

%% nl
theta1_spc = linspace(0.0001, 0.0004, 200);       
theta2_spc = linspace(0.1, 0.3, 200); 
theta3_spc = linspace(2,5, 200);
randIdx = randi(200,4,1);
% draw values and find kernel
p_nl = [theta1_spc(randIdx(1)) theta2_spc(randIdx(2)) theta3_spc(randIdx(3))];
model_nl = @(x) p_nl(1)* x.^2 + p_nl(2) * x + p_nl(3);  % define parameter space for nl
    
amp = 1e4;
%     CGP = conv( stim,amp * model_kernel(fliplr(t_kernel)),'full');
CGP = conv( stim,amp * fliplr( model_kernel(t_kernel) ),'full');
lambda = model_nl(CGP); 
spiketimes = simpp( lambda , dt_stim);
spiketimes(spiketimes > dt_stim * length(stim)) = [];

%% Shy's code
len = length(t_kernel);
time_v = (0:dt_stim:length(stim)*dt_stim-dt_stim);
r = histc(spiketimes,time_v);
S = [zeros(len-1,1); stim];
S_len = zeros(length(stim),length(t_kernel));
len = length(t_kernel);
for k = 1:length(r)
    S_len( k,: ) =  S( k:k+len-1 );
%     disp(['iteration ',num2str(k)]);
end

h=normc(pinv(S_len'*S_len)*(S_len'*(r))); %find kernel using linear regression
u=S_len*h;%linear prediction
u=u/std(u);%normalize the linear prediction
p=prctile(u',[0:0.01:1]*100);%calculate 100 precentiles of the linear prediction
p(end)=p(end)+eps;
[N,bin]=histc(u,p);
N=N(1:end-1);
m1=[];
x=[];
for i=1:length(N)
    temp=find(bin==i);
    x(i)=mean(u(temp));
    m1(i)=mean(r_cell(temp));
end
f=polyfit(x,m1,polyorder);%calculate polynomial fit
h=(pinv(S_len'*S_len)*(S_len'*u)); %filter
z=polyval(f,S_len*h); %predicted rate
%%
domain = (min(CGP):max(CGP));
figure();
scatter(x,m1,'r');
hold on;
plot(x,polyval(f,x))
hold off;

figure()
plot(z,'c');hold on;plot(lambda/std(lambda),'-r'); hold off; axis tight;