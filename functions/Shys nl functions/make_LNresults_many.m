%make_LNresults_many
if ~exist('S')|~exist('r') % load data if no data currently loaded
   cd c:/shy/data/Brown
   [filename,pathname]=uigetfile('bin*.mat', 'Choose Binned Data File');
   load([pathname filename])
   cd c:/shy/population
end

global S r
global d0 d1 d2 d3 x_z

S_hold=S;
error=[];
poisson_like=[];
len_cells=[4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4];
poly_order=[0 1 2 1 2 1 4 3 2 2 1 1 0 1 1 1 1 1 0 1 2 1 1];

quality=[0.08 0.39 0.43 0.15 0.26 0.65 0.65 0.52 0.56...
      0.27 0.07 0.23 0.48 0.20 0.16 0.19 0.31 0.27];
letters=['a','b','c','d','e','f'];
figure(1)
clf
plot_num=0;
for cell=[3 7 8 11 16  14]
   plot_num=plot_num+1;
   r_cell=r(:,cell);
   len=len_cells(cell);
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
   f=polyfit(x,m1,poly_order(cell));
   u=u/std(u);
   h=(pinv(S_len'*S_len)*(S_len'*u)); %filter
   temp=find(bin>2&bin<(max(bin)-1));%exclude data from boundaries
   error(len)=norm(r_cell(temp)-polyval(f,u(temp)));
   poisson_like(len)=sum(log(poisspdf(r_cell(temp),max(polyval(f,u(temp)),eps))));
   z=polyval(f,S_len*h); %predicted r   
   
   figure(1);
   subplot(2,3,plot_num) %plot nonlinearity 
   cla
   hold on
   plot(x,m1*20,'k.')
   plot(x,20*polyval(f,x),'k');
   if plot_num==1
	   legend('Estimated f','Polynomial fit')
   end
   xlabel('Linear prediction');
   ylabel('Expected firing rate(Hz)');
   ttl=title(['']);
   pos1=get(gca,'xlim');
   pos2=get(ttl,'position');
   tex=text(pos1(1),pos2(2),letters(plot_num),'VerticalAlignment', 'bottom','FontSize',14);
   box on
   
   %neural noise calculation
   x_z=prctile(z',[1:99]);
   x_z(end)=x_z(end)+eps;
   %[N,bin]=histc(z,x_z);
   
   means=[];vars=[];d0=[];d1=[];d2=[];d3=[];
   for i=1:length(x_z)-1
      temp=find(z>x_z(i)&z<x_z(i+1));
      x_z(i)=mean(z(temp));
      means(i)=mean(r_cell(temp));
      vars(i)=var(r_cell(temp));
      d0(i)=sum(r_cell(temp)==0)/length(temp);
      d1(i)=sum(r_cell(temp)==1)/length(temp);
      d2(i)=sum(r_cell(temp)==2)/length(temp);
      d3(i)=sum(r_cell(temp)==3)/length(temp);
            d4(i)=sum(r_cell(temp)==4)/length(temp);
   end
   x_z=x_z(1:end-1);
   
   sigma=fminbnd('fit_NGsigma',0.5,1.5)
   
   allx=[NGpdf(0,x_z',sigma) NGpdf(1,x_z',sigma) NGpdf(2,x_z',sigma) NGpdf(3,x_z',sigma) NGpdf(4,x_z',sigma)]';   
   figure(2)
   subplot(2,3,plot_num)
   cla
   
   ER=x_z*20;
   plot(ER,d0,'k.',ER,d1,'k+',ER,d2,'k*',ER,d3,'ko')
   hold on
   plot(ER,allx,'k')
   axis tight
   if plot_num==1
		legend('P(0)','P(1)','P(2)','P(3)','NG');
   end
   xlabel('Expected firing rate(Hz)','FontSize',14)
	ylabel('Conditional probability','FontSize',14)
   ttl=title(['']);
   pos1=get(gca,'xlim');
   pos2=get(ttl,'position');
   tex=text(pos1(1),pos2(2),letters(plot_num),'VerticalAlignment', 'bottom','FontSize',14);
   box on
end