%make_LNcomparison_all
%compares linear vs. nonlinear models for all cells in a bin* file.
%using penalized likelihood measures

% if ~exist('S')|~exist('r') % load data if no data currently loaded
cd c:/shy/data/Brown
[filename,pathname]=uigetfile('bin*.mat', 'Choose Binned Data File');
load([pathname filename])
cd c:/shy/population
% end

global S r
global d0 d1 d2 d3 x_z

S_hold=S;
error=[];
poisson_like=[];
NG_like=[];
len=4;

% letters=['a','b','c','d','e','f'];
figure(1)
clf
plot_num=0;
for cell=1:size(r,2)
    plot_num=plot_num+1;
    r_cell=r(:,cell);
    %    len=len_cells(cell);
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
    for poly_order=7:-1:0;
        f=polyfit(x,m1,poly_order);
        h=(pinv(S_len'*S_len)*(S_len'*u)); %filter
        temp=find(bin>2&bin<(max(bin)-1));%exclude data from boundaries
        error(len)=norm(r_cell(temp)-polyval(f,u(temp)));
        res=m1-polyval(f,x);
        z=polyval(f,u); %predicted r  
        %neural noise calculation
        x_z=prctile(z',[1:99]);
        x_z(end)=x_z(end)+eps;
        %[N,bin]=histc(z,x_z);
        if poly_order>0
            means=[];vars=[];d0=[];d1=[];d2=[];d3=[];
            for i=1:length(x_z)-1
                temp1=find(z>x_z(i)&z<x_z(i+1));
                x_z(i)=mean(z(temp1));
                means(i)=mean(r_cell(temp1));
                vars(i)=var(r_cell(temp1));
                d0(i)=sum(r_cell(temp1)==0)/length(temp1);
                d1(i)=sum(r_cell(temp1)==1)/length(temp1);
                d2(i)=sum(r_cell(temp1)==2)/length(temp1);
                d3(i)=sum(r_cell(temp1)==3)/length(temp1);
                d4(i)=sum(r_cell(temp1)==4)/length(temp1);
            end
            x_z=x_z(1:end-1);
            sigma(cell,poly_order+1)=fminbnd('fit_NGsigma',0.5,1.5);
        else
            sigma(cell,poly_order+1)=sigma(cell,poly_order+2);
        end
        
        poisson_like(cell,poly_order+1)=sum(log(max(real(poisspdf(r_cell(temp),polyval(f,u(temp)))),0.02)))-poly_order/2*log(length(temp));
        NG_like(cell,poly_order+1)=sum(log(max(NGpdf1(r_cell(temp)',polyval(f,u(temp))',sigma(cell,poly_order+1)),0.02)))-(poly_order+1)/2*log(length(temp));
        %gauss_like(cell,poly_order+1)=-sum((res).^2/(2*var(res)))-poly_order/2*log(length(res));
        z=polyval(f,S_len*h); %predicted r   
    end
    
    figure(1);
    subplot(6,6,plot_num) %plot nonlinearity 
    cla
    hold on
    plot(x,m1*20,'k.')
    plot(x,20*polyval(f,x),'k');
    if plot_num==1
        legend('Estimated f','Polynomial fit')
    end
    xlabel('Linear prediction');
    ylabel('Expected rate');
    ttl=title(['']);
    pos1=get(gca,'xlim');
    pos2=get(ttl,'position');
    %tex=text(pos1(1),pos2(2),letters(plot_num),'VerticalAlignment', 'bottom','FontSize',14);
    box on
end
[max_Poisson,optimal_order]=max(poisson_like(:,1:7)');
[max_NG,optimal_order2]=max(real(NG_like(:,1:7)'));
optimal_order=optimal_order-1
optimal_order2=optimal_order2-1