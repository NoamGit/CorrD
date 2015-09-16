function [f,x,m1,means,vars,x_z,d0,d1,d2,d3,poisson_like]=LN_bin(cell, len, polyorder)

global S r

r_cell=r(:,cell);
S_len=S(:,[2:(len+1) 12:(11+len)]);
h=normc(pinv(S_len'*S_len)*(S_len'*(r_cell))); %find kernel using linear regression
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
poisson_like=sum(log(max(real(poisspdf(r_cell,polyval(f,u))),0.02)));

%noise plot calculation
temp=find(bin>2&bin<(max(bin)-1));%exclude data from boundaries
x_z=prctile(u',[1:1:99]);%calculate precentiles of the nonlinear prediction
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
