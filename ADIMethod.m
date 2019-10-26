clear all

DA=1e-11;
DC=1e-11;
kr=6e4;
L=1e-6;
deltax=L/20; 
deltay=L/20;
deltaz=L/20;
deltat=5e-4;
tmax=1;
AvNumb=6.022e23;
NA0=600;
ConcA0=(NA0/(L^3/2))/AvNumb;
NC0=0;
ConcC0=0;
x=[1:L/deltax+2]';
y=[1:L/deltay+2]';
z=[1:L/deltaz+2]';
t=[0:deltat:tmax];

ConcA=zeros(length(t),length(x),length(y),length(z));
ConcC=zeros(length(t),length(x),length(y),length(z));

ConcAs=zeros(length(t),length(x)-2,length(y)-2,length(z)-2);
ConcAs2=zeros(length(t),length(x)-2,length(y)-2,length(z)-2);
ConcAs3=zeros(length(t),length(x)-2,length(y)-2,length(z)-2);

ConcCs=zeros(length(t),length(x)-2,length(y)-2,length(z)-2);
ConcCs2=zeros(length(t),length(x)-2,length(y)-2,length(z)-2);
ConcCs3=zeros(length(t),length(x)-2,length(y)-2,length(z)-2);


    Ax=zeros(length(t),length(x),length(y),length(z));
    Ay=zeros(length(t),length(x),length(y),length(z));
    Az=zeros(length(t),length(x),length(y),length(z));
    
    Cx=zeros(length(t),length(x),length(y),length(z));
    Cy=zeros(length(t),length(x),length(y),length(z));
    Cz=zeros(length(t),length(x),length(y),length(z));    

g1=zeros(length(t),length(x)-2,length(y)-2,length(z)-2);
g2=zeros(length(t),length(x)-2,length(y)-2,length(z)-2);
g3=zeros(length(t),length(x)-2,length(y)-2,length(z)-2);

g1c=zeros(length(t),length(x)-2,length(y)-2,length(z)-2);
g2c=zeros(length(t),length(x)-2,length(y)-2,length(z)-2);
g3c=zeros(length(t),length(x)-2,length(y)-2,length(z)-2);

SumA=zeros(length(t),1);
SumC=zeros(length(t),1);

SumA(1,1)=NA0;
SumC(1,1)=NC0;

nux=DA*deltat/deltax^2;
nuy=DA*deltat/deltay^2;
nuz=DA*deltat/deltaz^2;

nuxc=DC*deltat/deltax^2;
nuyc=DC*deltat/deltay^2;
nuzc=DC*deltat/deltaz^2;

for i=2:(length(x))/2
    for j=2:length(y)-1
        for k=2:length(z)-1
            ConcA(1,i,j,k)=ConcA0;
        end
    end
end


for i=2:length(x)-1
    for j=2:length(y)-1
        for k=2:length(z)-1
            ConcC(1,i,j,k)=ConcC0;
        end
    end
end


for k=1:length(t)-1
    
       for jj=3:length(y)-2
            for kk=3:length(z)-2 
                ConcA(k,1,jj,kk)=ConcA(k,3,jj,kk);
                ConcA(k,length(x),jj,kk)=ConcA(k,length(x)-2,jj,kk);  
                ConcC(k,1,jj,kk)=ConcC(k,3,jj,kk);
                ConcC(k,length(x),jj,kk)=ConcC(k,length(x)-2,jj,kk);                 
            end
        end
        
        for ii=3:length(x)-2
            for jj=3:length(y)-2
                ConcA(k,ii,jj,1)=ConcA(k,ii,jj,3);
                ConcA(k,ii,jj,length(z))=ConcA(k,ii,jj,length(z)-2);
                ConcC(k,ii,jj,1)=ConcC(k,ii,jj,3);
                ConcC(k,ii,jj,length(z))=ConcC(k,ii,jj,length(z)-2);                
            end
        end
        
        for ii=3:length(x)-2
            for kk=3:length(z)-2
                ConcA(k,ii,1,kk)=ConcA(k,ii,3,kk); 
                ConcA(k,ii,length(y),kk)=ConcA(k,ii,length(y)-2,kk);
                ConcC(k,ii,1,kk)=ConcC(k,ii,3,kk); 
                ConcC(k,ii,length(y),kk)=ConcC(k,ii,length(y)-2,kk);                
            end
        end
        
        for jj=3:length(y)-2
            for kk=3:length(z)-2               
                ConcA(k,1,2,kk)=ConcA(k,3,2,kk);
                ConcA(k,1,jj,2)=ConcA(k,3,jj,2);
                ConcA(k,1,length(y)-1,kk)=ConcA(k,3,length(y)-1,kk);
                ConcA(k,1,jj,length(z)-1)=ConcA(k,3,jj,length(z)-1);
 
                ConcA(k,length(x),2,kk)=ConcA(k,length(x)-2,2,kk);
                ConcA(k,length(x),jj,2)=ConcA(k,length(x)-2,jj,2);
                ConcA(k,length(x),length(y)-1,kk)=ConcA(k,length(x)-2,length(y)-1,kk);
                ConcA(k,length(x),jj,length(z)-1)=ConcA(k,length(x)-2,jj,length(z)-1);  

                ConcC(k,1,2,kk)=ConcC(k,3,2,kk);
                ConcC(k,1,jj,2)=ConcC(k,3,jj,2);
                ConcC(k,1,length(y)-1,kk)=ConcC(k,3,length(y)-1,kk);
                ConcC(k,1,jj,length(z)-1)=ConcC(k,3,jj,length(z)-1);
 
                ConcC(k,length(x),2,kk)=ConcC(k,length(x)-2,2,kk);
                ConcC(k,length(x),jj,2)=ConcC(k,length(x)-2,jj,2);
                ConcC(k,length(x),length(y)-1,kk)=ConcC(k,length(x)-2,length(y)-1,kk);
                ConcC(k,length(x),jj,length(z)-1)=ConcC(k,length(x)-2,jj,length(z)-1);   
            end
        end
     for ii=3:length(x)-2
            for kk=3:length(z)-2              
                ConcA(k,2,1,kk)=ConcA(k,2,3,kk);
                ConcA(k,ii,1,2)=ConcA(k,ii,3,2);
                ConcA(k,length(x)-1,1,kk)=ConcA(k,length(x)-1,3,kk);
                ConcA(k,ii,1,length(z)-1)=ConcA(k,ii,3,length(z)-1);
                
                ConcA(k,2,length(y),kk)=ConcA(k,2,length(y)-2,kk);
                ConcA(k,ii,length(y),2)=ConcA(k,ii,length(y)-2,2);
                ConcA(k,length(x)-1,length(y),kk)=ConcA(k,length(x)-1,length(y)-2,kk);
                ConcA(k,ii,length(y),length(z)-1)=ConcA(k,ii,length(y)-2,length(z)-1);  
                
                ConcC(k,2,1,kk)=ConcC(k,2,3,kk);
                ConcC(k,ii,1,2)=ConcC(k,ii,3,2);
                ConcC(k,length(x)-1,1,kk)=ConcC(k,length(x)-1,3,kk);
                ConcC(k,ii,1,length(z)-1)=ConcC(k,ii,3,length(z)-1);
                
                ConcC(k,2,length(y),kk)=ConcC(k,2,length(y)-2,kk);
                ConcC(k,ii,length(y),2)=ConcC(k,ii,length(y)-2,2);
                ConcC(k,length(x)-1,length(y),kk)=ConcC(k,length(x)-1,length(y)-2,kk);
                ConcC(k,ii,length(y),length(z)-1)=ConcC(k,ii,length(y)-2,length(z)-1);                 
            end
     end

     for ii=3:length(x)-2
        for jj=3:length(y)-2
                ConcA(k,2,jj,1)=ConcA(k,2,jj,3);
                ConcA(k,ii,2,1)=ConcA(k,ii,2,3);
                ConcA(k,length(x)-1,jj,1)=ConcA(k,length(x)-1,jj,3);
                ConcA(k,ii,length(y)-1,1)=ConcA(k,ii,length(y)-1,3);
                
                ConcA(k,2,jj,length(z))=ConcA(k,2,jj,length(z)-2);
                ConcA(k,ii,2,length(z))=ConcA(k,ii,2,length(z)-2);
                ConcA(k,length(x)-1,jj,length(z))=ConcA(k,length(x)-1,jj,length(z)-2);
                ConcA(k,ii,length(y)-1,length(z))=ConcA(k,ii,length(y)-1,length(z)-2);
                
                ConcC(k,2,jj,1)=ConcC(k,2,jj,3);
                ConcC(k,ii,2,1)=ConcC(k,ii,2,3);
                ConcC(k,length(x)-1,jj,1)=ConcC(k,length(x)-1,jj,3);
                ConcC(k,ii,length(y)-1,1)=ConcC(k,ii,length(y)-1,3);
                
                ConcC(k,2,jj,length(z))=ConcC(k,2,jj,length(z)-2);
                ConcC(k,ii,2,length(z))=ConcC(k,ii,2,length(z)-2);
                ConcC(k,length(x)-1,jj,length(z))=ConcC(k,length(x)-1,jj,length(z)-2);
                ConcC(k,ii,length(y)-1,length(z))=ConcC(k,ii,length(y)-1,length(z)-2);                
        end
     end
     
                ConcA(k,1,2,2)=ConcA(k,3,2,2);
                ConcA(k,2,1,2)=ConcA(k,2,3,2);
                ConcA(k,2,2,1)=ConcA(k,2,2,3);
                
                ConcA(k,1,length(y)-1,2)=ConcA(k,3,length(y)-1,2);
                ConcA(k,2,length(y),2)=ConcA(k,2,length(y)-2,2);
                ConcA(k,2,length(y)-1,1)=ConcA(k,2,length(y)-1,3);
                
                ConcA(k,1,2,length(z)-1)=ConcA(k,3,2,length(z)-1);
                ConcA(k,2,1,length(z)-1)=ConcA(k,2,3,length(z)-1);
                ConcA(k,2,2,length(z))=ConcA(k,2,2,length(z)-2);
                
                ConcA(k,1,length(y)-1,length(z)-1)=ConcA(k,3,length(y)-1,length(z)-1);
                ConcA(k,2,length(y),length(z)-1)=ConcA(k,2,length(y)-2,length(z)-1);
                ConcA(k,2,length(y)-1,length(z))=ConcA(k,2,length(y)-1,length(z)-2);
                
                ConcA(k,length(x),2,2)=ConcA(k,length(x)-2,2,2);
                ConcA(k,length(x)-1,1,2)=ConcA(k,length(x)-1,3,2);
                ConcA(k,length(x)-1,2,1)=ConcA(k,length(x)-1,2,3);    
                
                ConcA(k,length(x),length(y)-1,2)=ConcA(k,length(x)-2,length(y)-1,2);
                ConcA(k,length(x)-1,length(y),2)=ConcA(k,length(x)-1,length(y)-2,2);
                ConcA(k,length(x)-1,length(y)-1,1)=ConcA(k,length(x)-1,length(y)-1,3);
                
                ConcA(k,length(x),2,length(z)-1)=ConcA(k,length(x)-2,2,length(z)-1);
                ConcA(k,length(x)-1,1,length(z)-1)=ConcA(k,length(x)-1,3,length(z)-1);
                ConcA(k,length(x)-1,2,length(z))= ConcA(k,length(x)-1,2,length(z)-2);    
                
                ConcA(k,length(x),length(y)-1,length(z)-1)=ConcA(k,length(x)-2,length(y)-1,length(z)-1);
                ConcA(k,length(x)-1,length(y),length(z)-1)=ConcA(k,length(x)-1,length(y)-2,length(z)-1);
                ConcA(k,length(x)-1,length(y)-1,length(z))=ConcA(k,length(x)-1,length(y)-1,length(z)-2);             
    
                ConcC(k,1,2,2)=ConcC(k,3,2,2);
                ConcC(k,2,1,2)=ConcC(k,2,3,2);
                ConcC(k,2,2,1)=ConcC(k,2,2,3);
                
                ConcC(k,1,length(y)-1,2)=ConcC(k,3,length(y)-1,2);
                ConcC(k,2,length(y),2)=ConcC(k,2,length(y)-2,2);
                ConcC(k,2,length(y)-1,1)=ConcC(k,2,length(y)-1,3);
                
                ConcC(k,1,2,length(z)-1)=ConcC(k,3,2,length(z)-1);
                ConcC(k,2,1,length(z)-1)=ConcC(k,2,3,length(z)-1);
                ConcC(k,2,2,length(z))=ConcC(k,2,2,length(z)-2);
                
                ConcC(k,1,length(y)-1,length(z)-1)=ConcC(k,3,length(y)-1,length(z)-1);
                ConcC(k,2,length(y),length(z)-1)=ConcC(k,2,length(y)-2,length(z)-1);
                ConcC(k,2,length(y)-1,length(z))=ConcC(k,2,length(y)-1,length(z)-2);
                
                ConcC(k,length(x),2,2)=ConcC(k,length(x)-2,2,2);
                ConcC(k,length(x)-1,1,2)=ConcC(k,length(x)-1,3,2);
                ConcC(k,length(x)-1,2,1)=ConcC(k,length(x)-1,2,3);    
                
                ConcC(k,length(x),length(y)-1,2)=ConcC(k,length(x)-2,length(y)-1,2);
                ConcC(k,length(x)-1,length(y),2)=ConcC(k,length(x)-1,length(y)-2,2);
                ConcC(k,length(x)-1,length(y)-1,1)=ConcC(k,length(x)-1,length(y)-1,3);
                
                ConcC(k,length(x),2,length(z)-1)=ConcC(k,length(x)-2,2,length(z)-1);
                ConcC(k,length(x)-1,1,length(z)-1)=ConcC(k,length(x)-1,3,length(z)-1);
                ConcC(k,length(x)-1,2,length(z))= ConcC(k,length(x)-1,2,length(z)-2);    
                
                ConcC(k,length(x),length(y)-1,length(z)-1)=ConcC(k,length(x)-2,length(y)-1,length(z)-1);
                ConcC(k,length(x)-1,length(y),length(z)-1)=ConcC(k,length(x)-1,length(y)-2,length(z)-1);
                ConcC(k,length(x)-1,length(y)-1,length(z))=ConcC(k,length(x)-1,length(y)-1,length(z)-2);             
                               
    
    for ii=2:length(x)-1
        for jj=2:length(y)-1
            for kk=2:length(z)-1
                Ax(k,ii,jj,kk)=nux*(ConcA(k,ii-1,jj,kk)-2*ConcA(k,ii,jj,kk)+ConcA(k,ii+1,jj,kk));
                Ay(k,ii,jj,kk)=nuy*(ConcA(k,ii,jj-1,kk)-2*ConcA(k,ii,jj,kk)+ConcA(k,ii,jj+1,kk));
                Az(k,ii,jj,kk)=nuz*(ConcA(k,ii,jj,kk-1)-2*ConcA(k,ii,jj,kk)+ConcA(k,ii,jj,kk+1));
                g1(k,ii-1,jj-1,kk-1)=ConcA(k,ii,jj,kk)+Ax(k,ii,jj,kk)/2+Ay(k,ii,jj,kk)+Az(k,ii,jj,kk)-kr*deltat*(ConcA(k,ii,jj,kk))^2;

                Cx(k,ii,jj,kk)=nuxc*(ConcC(k,ii-1,jj,kk)-2*ConcC(k,ii,jj,kk)+ConcC(k,ii+1,jj,kk));
                Cy(k,ii,jj,kk)=nuyc*(ConcC(k,ii,jj-1,kk)-2*ConcC(k,ii,jj,kk)+ConcC(k,ii,jj+1,kk));
                Cz(k,ii,jj,kk)=nuzc*(ConcC(k,ii,jj,kk-1)-2*ConcC(k,ii,jj,kk)+ConcC(k,ii,jj,kk+1));
                g1c(k,ii-1,jj-1,kk-1)=ConcC(k,ii,jj,kk)+Cx(k,ii,jj,kk)/2+Cy(k,ii,jj,kk)+Cz(k,ii,jj,kk)+kr*deltat*(ConcA(k,ii,jj,kk))^2;
                
            end
        end
    end   
    
    a1=(1+nux)*ones(1,length(x)-2);
    b1=[-nux/2*ones(1,length(x)-3),-nux];
    c1=[-nux,-nux/2*ones(1,length(x)-3)];
    
    a1c=(1+nuxc)*ones(1,length(x)-2);
    b1c=[-nuxc/2*ones(1,length(x)-3),-nuxc];
    c1c=[-nuxc,-nuxc/2*ones(1,length(x)-3)];    
    
    for jj=1:length(y)-2
        for kk=1:length(z)-2
            ConcAs(k,:,jj,kk)=tridiag(a1,b1,c1,g1(k,:,jj,kk));
            ConcCs(k,:,jj,kk)=tridiag(a1c,b1c,c1c,g1c(k,:,jj,kk));
        end
    end
   

    %%%
    
    for ii=2:length(x)-1
        for jj=2:length(y)-1
            for kk=2:length(z)-1
                g2(k,ii-1,jj-1,kk-1)=ConcAs(k,ii-1,jj-1,kk-1)-Ay(k,ii,jj,kk)/2;
                g2c(k,ii-1,jj-1,kk-1)=ConcCs(k,ii-1,jj-1,kk-1)-Cy(k,ii,jj,kk)/2;
            end
        end
    end 
    
    a2=(1+nuy)*ones(1,length(y)-2);
    b2=[-nuy/2*ones(1,length(y)-3),-nuy];
    c2=[-nuy,-nuy/2*ones(1,length(y)-3)];
    
    a2c=(1+nuyc)*ones(1,length(y)-2);
    b2c=[-nuyc/2*ones(1,length(y)-3),-nuyc];
    c2c=[-nuyc,-nuyc/2*ones(1,length(y)-3)];
    
    for ii=1:length(x)-2
        for kk=1:length(z)-2
            ConcAs2(k,ii,:,kk)=tridiag(a2,b2,c2,g2(k,ii,:,kk));
            ConcCs2(k,ii,:,kk)=tridiag(a2c,b2c,c2c,g2c(k,ii,:,kk));
        end
    end    
    
   %%%
   
    for ii=2:length(x)-1
        for jj=2:length(y)-1
            for kk=2:length(z)-1
                g3(k,ii-1,jj-1,kk-1)=ConcAs2(k,ii-1,jj-1,kk-1)-Az(k,ii,jj,kk)/2+kr*(ConcA(k,ii,jj,kk))^2*deltat;
            end
        end
    end    
    
    for ii=1:length(x)-2
        for jj=1:length(x)-2
            for kk=1:length(x)-2
                a3(1,ii,jj,kk)=1+nuz+kr*ConcA(k,ii+1,jj+1,kk+1)*deltat;
            end
        end
    end
    b3=[-nuz/2*ones(1,length(z)-3),-nuz];
    c3=[-nuz,-nuz/2*ones(1,length(z)-3)];
   
    

    for ii=1:length(x)-2
        for jj=1:length(y)-2
            aa3=squeeze(a3(1,ii,jj,:));
            ConcAs3(k,ii,jj,:)=tridiag(aa3',b3,c3,g3(k,ii,jj,:));
        end
    end

    for ii=1:length(x)-2
        for jj=1:length(y)-2
            for kk=1:length(z)-2
                ConcA(k+1,ii+1,jj+1,kk+1)=ConcAs3(k,ii,jj,kk);
            end
        end
    end
            
      
 SumA(k+1,1)=sum(sum(sum(ConcA(k+1,2:(length(x)/2),2:length(y)-1,2:length(z)-1)*deltax*deltay*deltaz*AvNumb)));   
     

   %%%
   
    for ii=2:length(x)-1
        for jj=2:length(y)-1
            for kk=2:length(z)-1
                g3(k,ii-1,jj-1,kk-1)=ConcCs2(k,ii-1,jj-1,kk-1)-Cz(k,ii,jj,kk)/2+kr*((ConcA(k+1,ii,jj,kk))^2-(ConcA(k,ii,jj,kk))^2)*deltat/2;
            end
        end
    end    
    
    a3=(1+nuz)*ones(1,length(z)-2);
    b3=[-nuz/2*ones(1,length(z)-3),-nuz];
    c3=[-nuz,-nuz/2*ones(1,length(z)-3)];
   
    

    for ii=1:length(x)-2
        for jj=1:length(y)-2
            ConcCs3(k,ii,jj,:)=tridiag(a3,b3,c3,g3(k,ii,jj,:));
        end
    end

    for ii=1:length(x)-2
        for jj=1:length(y)-2
            for kk=1:length(z)-2
                ConcC(k+1,ii+1,jj+1,kk+1)=ConcCs3(k,ii,jj,kk);
            end
        end
    end
            
      
 SumC(k+1,1)=sum(sum(sum(ConcC(k+1,2:(length(x)/2),2:length(y)-1,2:length(z)-1)*deltax*deltay*deltaz*AvNumb)));   
    
end

figure(1)
plot(t,SumC,'k','Linewidth',1.5) 
axis([0 tmax 0 300])
hold on 
legend({'PDE Solution'},'Location','southeast','NumColumns',1)

   
xlabel(['Time (s)'])
ylabel(['# Molecules of C'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
    print('ADIMethod', '-dpdf', '-r300');
