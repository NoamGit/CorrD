%make_LNtestplot

%[filename,pathname]=uigetfile('bin*.mat', 'Choose Binned Data File');
%load([pathname filename])
global S r
h_truth=([7.9101   -9.9789    1.5759   -8.4679   21.9350  -14.1285]');
len=3;
S_3=S_hold(:,[2:(len+1) 12:(11+len)]);
u=S_3*h_truth;
h_truth=h_truth/std(u);
u=u/std(u);
real_LN=1.1+4*tanh(S_3*h_truth*2)/4;
r_cell=poissrnd(real_LN);
error=[];
poisson_like=[];
%r_cell=(real_LN);

%try different kernel lengths
for len=1:8
   S_len=S_hold(:,[2:(len+1) 12:(11+len)]);
   h=normc(pinv(S_len'*S_len)*(S_len'*(r_cell))); %filter
   u=S_len*h;%linear prediction
   u=u/std(u);
   p=prctile(u',[0:0.01:1]*100);%calculate prencetiles, 100 bins in each
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
   oldh=h;
   f=polyfit(x,m1,9);
   u=u/std(u);
   h=(pinv(S_len'*S_len)*(S_len'*u)); %filter
   temp=find(bin>2&bin<(max(bin)-1));%exclude data from boundaries
   error(len)=norm(r_cell(temp)-polyval(f,u(temp)));
   poisson_like(len)=sum(log(poisspdf(r_cell(temp),max(polyval(f,u(temp)),eps))));
end

%use optimal kernel length
len=3;
S_len=S_hold(:,[2:(len+1) 12:(11+len)]);
h=normc(pinv(S_len'*S_len)*(S_len'*(r_cell))); %filter
u=S_len*h;%linear prediction
u=u/std(u);
p=prctile(u',[0:0.01:1]*100);%calculate prencetiles, 100 bins in each
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
oldh=h;
f=polyfit(x,m1,9);
u=u/std(u);
h=(pinv(S_len'*S_len)*(S_len'*u)); %filter
temp=find(bin>2&bin<(max(bin)-1));%exclude data from boundaries
error(len)=norm(r_cell(temp)-polyval(f,u(temp)));
poisson_like(len)=sum(log(poisspdf(r_cell(temp),max(polyval(f,u(temp)),eps))));


figure(2)
clf

subplot(2,2,1)
cla
hold on
bar(linspace(0,1.999,40),r_cell(1:40),'w')
plot(linspace(0,1.999,40),real_LN(1:40),'k')
legend('Firing rate','Spike counts');
xlabel('time(sec)')
axis([-0.1 2.1 0 max(r_cell(1:40))+0.5])
ttl=title('');
pos1=get(gca,'xlim');
pos2=get(ttl,'position');
tex=text(pos1(1),pos2(2),'a','VerticalAlignment', 'bottom','FontSize',14);
box on


subplot(2,2,2)
cla
x_min=min(S_3*h_truth);
x_max=max(S_3*h_truth);
x_LN=linspace(x_min,x_max,100);
y_LN=1.1+4*tanh(x_LN*2)/4;
plot(x_LN,y_LN,'r','MarkerSize',2);
hold on
plot(x,m1,'.')
plot(x,polyval(f,x));
legend('f','Estimated f','Polynomial fit')
xlabel('Linear output');
ylabel('Nonlinear output');
axis ('tight')
ttl=title('');
pos1=get(gca,'xlim');
pos2=get(ttl,'position');
tex=text(pos1(1),pos2(2),'b','VerticalAlignment', 'bottom','FontSize',14);

subplot(2,2,3)
cla
plot(error,'k*-')
xlabel('Kernel Length')
ylabel('Error')
ttl=title('');
pos1=get(gca,'xlim');
pos2=get(ttl,'position');
tex=text(pos1(1),pos2(2),'c','VerticalAlignment', 'bottom','FontSize',14);



z=polyval(f,S_len*h); %predicted r
x_z=prctile(u',[1:1:99]);
x_z(end)=x_z(end)+eps;
[N,bin]=histc(u,x_z);
x_z=x_z(1:end-1);
means=[];vars=[];d0=[];d1=[];d2=[];d3=[];
for i=1:length(x_z)
   temp=find(bin==i);
   x_z(i)=mean(z(temp));
   means(i)=mean(r_cell(temp));
   vars(i)=var(r_cell(temp));
   d0(i)=sum(r_cell(temp)==0)/length(temp);
   d1(i)=sum(r_cell(temp)==1)/length(temp);
   d2(i)=sum(r_cell(temp)==2)/length(temp);
   d3(i)=sum(r_cell(temp)==3)/length(temp);
end

[m,I]=sort(m1);

subplot(2,2,4)
cla
hold on
plot(x_z,d0,'.',x_z,d1,'.',x_z,d2,'.',x_z,d3,'.') 
plot(x_z,poisspdf(0,x_z),'k--',x_z,poisspdf(1,x_z),'k--',x_z,poisspdf(2,x_z),'k--',x_z,poisspdf(3,x_z),'k--')
legend('P(0)','P(1)','P(2)','P(3)','Poisson');
xlabel('Expected firing rate')
ylabel('Conditional distributions')
axis ('tight')
ttl=title('');
pos1=get(gca,'xlim');
pos2=get(ttl,'position');
tex=text(pos1(1),pos2(2),'d','VerticalAlignment', 'bottom','FontSize',14);
box on
