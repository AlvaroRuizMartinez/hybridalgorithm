clear all

% Diffusion-Reaction Algorithm - Brownian motion / Gillespie Algorithm

for I=1:1000
I
DA=1e-11;
DB=1e-11;
DC=1e-11;
kr=6e4;
tmax=1;
L=1e-6;
AvNumb=6.022e23;
tt=0;
dd=[];
ddt=[];
ddtc=[];
Das=[];

NA0=600;
NB0=600;
dd(1,1)=600;
ddt(1,1)=600;
ddtc(1,1)=0;
Das(1,1)=0;
NC0=600;
DaL=kr*NA0/(6*DA*L*AvNumb);

xa=zeros(NA0,1);
ya=zeros(NA0,1);
za=zeros(NA0,1);
xb=zeros(NB0,1);
yb=zeros(NB0,1);
zb=zeros(NB0,1);
xc=zeros(NC0,1);
yc=zeros(NC0,1);
zc=zeros(NC0,1);

% Initial position of particles

for na=1:NA0
    xa(na,1)=L/2*rand;
    ya(na,1)=L*rand;
    za(na,1)=L*rand;
end

for nb=1:NB0
    xb(nb,1)=L/2*rand;
    yb(nb,1)=L*rand;
    zb(nb,1)=L*rand;
end

% Diffusion-Reaction
t2=0;
t2p=0;
t=0;
l=0;

while t<tmax
  
    l=l+1;  
    kmp=floor(nthroot(dd(l,1),3));

    deltax=L/kmp;
    Dass(l,1)=kr/(6*DA*deltax*AvNumb);
    
    %%%
   
    if Dass(l,1)>1
        kmp=2*dd(l,1);
        if kmp~=0
            kmp=kmp;
        else
            kmp=2;
        end
        km0=kmp;
        km1=1;
    else
        km0=kmp;
        km1=kmp;
    end
    
    deltat(l,1)=0;
    
    %%%
    
    while deltat(l,1)==0
        
    deltax=L/km0; 
    deltay=L/km1;
    deltaz=L/km1;
    
    x=[0:deltax:L]';
    y=[0:deltay:L]';
    z=[0:deltaz:L]';
    
    Numbcell=zeros(length(x)-1,length(y)-1,length(z)-1);
    count=1;
    Ncelltot=(length(x)-1)*(length(y)-1)*(length(z)-1);

    for i=1:length(x)-1
        for j=1:length(y)-1
            for k=1:length(z)-1
                Numbcell(i,j,k)=count;
                count=count+1;
            end
        end
    end

    nacell=zeros(Ncelltot,1);
    nbcell=zeros(Ncelltot,1);
    ConcAt=zeros(Ncelltot,1);
    ConcBt=zeros(Ncelltot,1);
    nprod=zeros(Ncelltot,1);
    for na=1:NA0
        if xa(na,l)~=0
            icellna(na,l)=floor((xa(na,l))/deltax)+1;
            jcellna(na,l)=floor((ya(na,l))/deltay)+1;
            kcellna(na,l)=floor((za(na,l))/deltaz)+1;
            Cellal(na,1)=Numbcell(icellna(na,l),jcellna(na,l),kcellna(na,l));
            nacell(Cellal(na,1),1)=nacell(Cellal(na,1),1)+1;
             else
            continue
        end
    end
    for nb=1:NB0
        if xb(nb,l)~=0 
            icellnb(nb,l)=floor((xb(nb,l))/deltax)+1;
            jcellnb(nb,l)=floor((yb(nb,l))/deltay)+1;
            kcellnb(nb,l)=floor((zb(nb,l))/deltaz)+1;
            Cellbl(nb,1)=Numbcell(icellnb(nb,l),jcellnb(nb,l),kcellnb(nb,l));
            nbcell(Cellbl(nb,1),1)=nbcell(Cellbl(nb,1),1)+1;
            else
            continue
        end
    end 
    
   for i=1:Ncelltot
       if nacell(i,1)~=0 && nbcell(i,1)~=0
          nprod(i,1)=nacell(i,1)*nbcell(i,1);
       else
          continue
       end
   end
   
   sn=sum(nprod(nprod==1));
   sn2=sum(nprod(nprod>1));  
   
    if sn<0.5*(sn+sn2) && Dass(l,1)<=1
        km1=km1+2;
        km0=km0+2;
    elseif sn<0.1*(sn+sn2) && Dass(l,1)>1
        km1=km1;
        km0=km0+2;   
    else           
        deltat(l,1)=0.63*deltax*deltay*deltaz*AvNumb/kr;
    end

    
    end

    rxa=zeros(na,1);
    rya=zeros(na,1);
    rza=zeros(na,1);
    rxb=zeros(nb,1);
    ryb=zeros(nb,1);
    rzb=zeros(nb,1);

    % Diffusion step
    
   nacellt=nacell;
   nbcellt=nbcell;   
 
        for na=1:NA0
            if xa(na,l)==0 
                xa(na,l+1)=0;
                ya(na,l+1)=0;
                za(na,l+1)=0;
            else
                rxa(na,1)=randn;
                rya(na,1)=randn;
                rza(na,1)=randn;
                xa(na,l+1)=xa(na,l)+sqrt(2*DA*deltat(l,1))*rxa(na,1);
                ya(na,l+1)=ya(na,l)+sqrt(2*DA*deltat(l,1))*rya(na,1);
                za(na,l+1)=za(na,l)+sqrt(2*DA*deltat(l,1))*rza(na,1);
                if xa(na,l+1)>0 && xa(na,l+1)<=L
                    xa(na,l+1)=xa(na,l+1);
                else
                    xa(na,l+1)=abs(L*(1-4*abs(1/2-(xa(na,l+1)/(4*L)+1/4-floor(xa(na,l+1)/(4*L)+1/4)))));
                end
                if ya(na,l+1)>0 && ya(na,l+1)<=L
                    ya(na,l+1)=ya(na,l+1);
                else
                    ya(na,l+1)=abs(L*(1-4*abs(1/2-(ya(na,l+1)/(4*L)+1/4-floor(ya(na,l+1)/(4*L)+1/4)))));
                end
                if za(na,l+1)>0 && za(na,l+1)<=L
                    za(na,l+1)=za(na,l+1);
                else
                    za(na,l+1)=abs(L*(1-4*abs(1/2-(za(na,l+1)/(4*L)+1/4-floor(za(na,l+1)/(4*L)+1/4)))));
                end
            end
        end
    
        for nb=1:NB0
            if xb(nb,l)==0
                xb(nb,l+1)=0;
                yb(nb,l+1)=0;
                zb(nb,l+1)=0;
            else
                rxb(nb,1)=randn;
                ryb(nb,1)=randn;
                rzb(nb,1)=randn;
                xb(nb,l+1)=xb(nb,l)+sqrt(2*DB*deltat(l,1))*rxb(nb,1);
                yb(nb,l+1)=yb(nb,l)+sqrt(2*DB*deltat(l,1))*ryb(nb,1);
                zb(nb,l+1)=zb(nb,l)+sqrt(2*DB*deltat(l,1))*rzb(nb,1);
                if xb(nb,l+1)>0 && xb(nb,l+1)<=L
                    xb(nb,l+1)=xb(nb,l+1);
                else
                    xb(nb,l+1)=abs(L*(1-4*abs(1/2-(xb(nb,l+1)/(4*L)+1/4-floor(xb(nb,l+1)/(4*L)+1/4)))));
                end
                if yb(nb,l+1)>0 && yb(nb,l+1)<=L
                    yb(nb,l+1)=yb(nb,l+1);
                else
                    yb(nb,l+1)=abs(L*(1-4*abs(1/2-(yb(nb,l+1)/(4*L)+1/4-floor(yb(nb,l+1)/(4*L)+1/4)))));
                end
                if zb(nb,l+1)>0 && zb(nb,l+1)<=L
                    zb(nb,l+1)=zb(nb,l+1);
                else
                    zb(nb,l+1)=abs(L*(1-4*abs(1/2-(zb(nb,l+1)/(4*L)+1/4-floor(zb(nb,l+1)/(4*L)+1/4)))));
                end
            end
        end
        
        for nc=1:NC0
            if xc(nc,l)==0
                xc(nc,l+1)=0;
                yc(nc,l+1)=0;
                zc(nc,l+1)=0;
            else
                rxc(nc,1)=randn;
                ryc(nc,1)=randn;
                rzc(nc,1)=randn;
                xc(nc,l+1)=xc(nc,l)+sqrt(2*DC*deltat(l,1))*rxc(nc,1);
                yc(nc,l+1)=yc(nc,l)+sqrt(2*DC*deltat(l,1))*ryc(nc,1);
                zc(nc,l+1)=zc(nc,l)+sqrt(2*DC*deltat(l,1))*rzc(nc,1);
                if xc(nc,l+1)>0 && xc(nc,l+1)<=L
                    xc(nc,l+1)=xc(nc,l+1);
                else
                    xc(nc,l+1)=abs(L*(1-4*abs(1/2-(xc(nc,l+1)/(4*L)+1/4-floor(xc(nc,l+1)/(4*L)+1/4)))));
                end
                if yc(nc,l+1)>0 && yc(nc,l+1)<=L
                    yc(nc,l+1)=yc(nc,l+1);
                else
                    yc(nc,l+1)=abs(L*(1-4*abs(1/2-(yc(nc,l+1)/(4*L)+1/4-floor(yc(nc,l+1)/(4*L)+1/4)))));
                end
                if zc(nc,l+1)>0 && zc(nc,l+1)<=L
                    zc(nc,l+1)=zc(nc,l+1);
                else
                    zc(nc,l+1)=abs(L*(1-4*abs(1/2-(zc(nc,l+1)/(4*L)+1/4-floor(zc(nc,l+1)/(4*L)+1/4)))));
                end
            end
        end
      
        
   nacelldt=zeros(Ncelltot,1);
   nbcelldt=zeros(Ncelltot,1);
   ConcAdt=zeros(Ncelltot,1);
   ConcBdt=zeros(Ncelltot,1);
   
    for na=1:NA0
        if xa(na,l+1)~=0
            icellna(na,l+1)=floor((xa(na,l+1))/deltax)+1;
            jcellna(na,l+1)=floor((ya(na,l+1))/deltay)+1;
            kcellna(na,l+1)=floor((za(na,l+1))/deltaz)+1;
            Cellalp(na,1)=Numbcell(icellna(na,l+1),jcellna(na,l+1),kcellna(na,l+1));
            nacelldt(Cellalp(na,1),1)=nacelldt(Cellalp(na,1),1)+1;
              else
            continue
        end
    end
    for nb=1:NB0
        if xb(nb,l+1)~=0
            icellnb(nb,l+1)=floor((xb(nb,l+1))/deltax)+1;
            jcellnb(nb,l+1)=floor((yb(nb,l+1))/deltay)+1;
            kcellnb(nb,l+1)=floor((zb(nb,l+1))/deltaz)+1;
            Cellblp(nb,1)=Numbcell(icellnb(nb,l+1),jcellnb(nb,l+1),kcellnb(nb,l+1));  
            nbcelldt(Cellblp(nb,1),1)=nbcelldt(Cellblp(nb,1),1)+1;
           else
            continue
        end
    end 
    
    ConcAt(:,1)=nacellt(:,1)/(AvNumb*deltax*deltay*deltaz);
    ConcBt(:,1)=nbcellt(:,1)/(AvNumb*deltax*deltay*deltaz);
    ConcAdt(:,1)=nacelldt(:,1)/(AvNumb*deltax*deltay*deltaz);
    ConcBdt(:,1)=nbcelldt(:,1)/(AvNumb*deltax*deltay*deltaz);  

   % Reaction step
   
   xap=zeros(NA0,1);
   yap=zeros(NA0,1);
   zap=zeros(NA0,1);
   xbp=zeros(NB0,1);
   ybp=zeros(NB0,1);
   zbp=zeros(NB0,1);  
   Ar=zeros(Ncelltot,1);
   Br=zeros(Ncelltot,1);
  

      rr=isempty(nonzeros(nprod));
   if rr==true
       trmin(l,1)=10*deltax^2/(2*DA*3);
       Dasp(l,1)=10;
       %ratio(l,1)=1;
       
   else
       trmin(l,1)=deltax*deltay*deltaz*AvNumb/(kr*max(nonzeros(nprod)));
       Dasp(l,1)=(sqrt(2*3*DA*trmin(l,1))/deltax)^2;
   end

   na=0;
   nb=0;
   ttim=t*ones(Ncelltot,1);
   told=ttim;
   deltat(l,1);
   Das(l,1)=1/Dasp(l,1);

   while na+nb<NA0+NB0
       RN=rand;
   if RN<(NA0-na)/(NA0-na+NB0-nb) && na<NA0
     na=na+1;
         if xa(na,l)~=0
             ran=rand;
             if Das(l,1)>0.1
                deltatpa=ran/(AvNumb*deltax*deltay*deltaz*(0.37*ConcAt(Cellal(na,1),1)+0.63*ConcAdt(Cellalp(na,1),1))*kr*(0.37*ConcBt(Cellal(na,1),1)+0.63*ConcBdt(Cellalp(na,1),1)));                   
             else
                deltatpa=log(1/ran)/(AvNumb*deltax*deltay*deltaz*(ConcAt(Cellal(na,1),1)+ConcAdt(Cellalp(na,1),1))/2*kr*(ConcBt(Cellal(na,1),1)+ConcBdt(Cellalp(na,1),1))/2);                   
             end
             if deltatpa<deltat(l,1) && deltatpa>0
                xap(na,1)=xa(na,l)+sqrt(2*DA*deltatpa)*rxa(na,1);
                yap(na,1)=ya(na,l)+sqrt(2*DA*deltatpa)*rya(na,1);
                zap(na,1)=za(na,l)+sqrt(2*DA*deltatpa)*rza(na,1);
                if xap(na,1)>=0 && xap(na,1)<=L
                   xap(na,1)=xap(na,1);
                else
                   xap(na,1)=abs(L*(1-4*abs(1/2-(xap(na,1)/(4*L)+1/4-floor(xap(na,1)/(4*L)+1/4)))));
                end
                if yap(na,1)>=0 && yap(na,1)<=L
                   yap(na,1)=yap(na,1);
                else
                   yap(na,1)=abs(L*(1-4*abs(1/2-(yap(na,1)/(4*L)+1/4-floor(yap(na,1)/(4*L)+1/4)))));
                end
                if zap(na,1)>=0 && zap(na,1)<=L
                   zap(na,1)=zap(na,1);
                else
                   zap(na,1)=abs(L*(1-4*abs(1/2-(zap(na,1)/(4*L)+1/4-floor(zap(na,1)/(4*L)+1/4)))));
                end
             elseif deltatpa>=deltat(l,1)
                xap(na,1)=xa(na,l+1);
                yap(na,1)=ya(na,l+1);
                zap(na,1)=za(na,l+1);
                deltatpa=deltat(l,1);  
             else
                 continue
             end
                icellnap(na,1)=floor((xap(na,1))/deltax)+1;
                jcellnap(na,1)=floor((yap(na,1))/deltay)+1;
                kcellnap(na,1)=floor((zap(na,1))/deltaz)+1;
                Cella=Numbcell(icellnap(na,1),jcellnap(na,1),kcellnap(na,1));
                if ttim(Cella,1)<deltat(l,1)+told(Cella,1) %&& deltatpa==deltat(l,1)
                Ar(Cella,1)=Ar(Cella,1)+1;
                nareaccellp(Cella,Ar(Cella,1),1)=na;
                timea(Cella,Ar(Cella,1),1)=deltatpa;
                if Ar(Cella,1)>0 && Br(Cella,1)>0 
                    nnb=nbreaccellp(Cella,Br(Cella,1),1);
                    ts1=min(timea(Cella,Ar(Cella,1),1),timeb(Cella,Br(Cella,1),1));
                    ttim(Cella,1)=ts1+ttim(Cella,1);
                    xc(na,l+1)=(xap(na,1)+xbp(nnb,1))/2;
                    yc(na,l+1)=(yap(na,1)+ybp(nnb,1))/2;
                    zc(na,l+1)=(zap(na,1)+zbp(nnb,1))/2;  
                    xa(na,l+1)=0;
                    ya(na,l+1)=0;
                    za(na,l+1)=0;
                    xb(nnb,l+1)=0;
                    yb(nnb,l+1)=0;
                    zb(nnb,l+1)=0;
                    Ar(Cella,1)=Ar(Cella,1)-1;
                    Br(Cella,1)=Br(Cella,1)-1;                    
                    ConcAdt(Cellalp(na,1),1)=(nacelldt(Cellalp(na,1),1)-1)/(AvNumb*deltax*deltay*deltaz);
                    ConcBdt(Cellblp(nnb,1),1)=(nbcelldt(Cellblp(nnb,1),1)-1)/(AvNumb*deltax*deltay*deltaz); 
                    nacelldt(Cellalp(na,1),1)=nacelldt(Cellalp(na,1),1)-1;
                    nbcelldt(Cellblp(nnb,1),1)=nbcelldt(Cellblp(nnb,1),1)-1;                
                    
                    else
                    continue
                end
                else
                    continue
                end

         else
           continue
         end
   elseif RN>(NA0-na)/(NA0-na+NB0-nb) && nb<NB0
       nb=nb+1;

         if xb(nb,l)~=0
             ran=rand;
            if Das(l,1)>0.1
                deltatpb=ran/(AvNumb*deltax*deltay*deltaz*(0.37*ConcAt(Cellbl(nb,1),1)+0.63*ConcAdt(Cellblp(nb,1),1))*kr*(0.37*ConcBt(Cellbl(nb,1),1)+0.63*ConcBdt(Cellblp(nb,1),1))); 
            else
                deltatpb=log(1/ran)/(AvNumb*deltax*deltay*deltaz*(ConcAt(Cellbl(nb,1),1)+ConcAdt(Cellblp(nb,1),1))/2*kr*(ConcBt(Cellbl(nb,1),1)+ConcBdt(Cellblp(nb,1),1))/2);                   
            end
             if deltatpb<deltat(l,1) && deltatpb>0
                xbp(nb,1)=xb(nb,l)+sqrt(2*DB*deltatpb)*rxb(nb,1);
                ybp(nb,1)=yb(nb,l)+sqrt(2*DB*deltatpb)*ryb(nb,1);
                zbp(nb,1)=zb(nb,l)+sqrt(2*DB*deltatpb)*rzb(nb,1);
                if xbp(nb,1)>=0 && xbp(nb,1)<=L
                   xbp(nb,1)=xbp(nb,1);
                else
                   xbp(nb,1)=abs(L*(1-4*abs(1/2-(xbp(nb,1)/(4*L)+1/4-floor(xbp(nb,1)/(4*L)+1/4)))));
                end
                if ybp(nb,1)>=0 && ybp(nb,1)<=L
                   ybp(nb,1)=ybp(nb,1);
                else
                   ybp(nb,1)=abs(L*(1-4*abs(1/2-(ybp(nb,1)/(4*L)+1/4-floor(ybp(nb,1)/(4*L)+1/4)))));
                end
                if zbp(nb,1)>=0 && zbp(nb,1)<=L
                   zbp(nb,1)=zbp(nb,1);
                else
                   zbp(nb,1)=abs(L*(1-4*abs(1/2-(zbp(nb,1)/(4*L)+1/4-floor(zbp(nb,1)/(4*L)+1/4)))));
                end 
             elseif deltatpb>=deltat(l,1)
                xbp(nb,1)=xb(nb,l+1);
                ybp(nb,1)=yb(nb,l+1);
                zbp(nb,1)=zb(nb,l+1);
                deltatpb=deltat(l,1); 
             else
                 continue
             end
                icellnbp(nb,1)=floor((xbp(nb,1))/deltax)+1;
                jcellnbp(nb,1)=floor((ybp(nb,1))/deltay)+1;
                kcellnbp(nb,1)=floor((zbp(nb,1))/deltaz)+1;
                Cellb=Numbcell(icellnbp(nb,1),jcellnbp(nb,1),kcellnbp(nb,1));
                if ttim(Cellb,1)<deltat(l,1)+told(Cellb,1) %&& deltatpb==deltat(l,1)      
                Br(Cellb,1)=Br(Cellb,1)+1; 
                nbreaccellp(Cellb,Br(Cellb,1),1)=nb;
                timeb(Cellb,Br(Cellb,1),1)=deltatpb;
                if Ar(Cellb,1)>0 && Br(Cellb,1)>0
                    nna=nareaccellp(Cellb,Ar(Cellb,1),1);            
                    ts2=min(timea(Cellb,Ar(Cellb,1),1),timeb(Cellb,Br(Cellb,1),1));
                    ttim(Cellb,1)=ts2+ttim(Cellb,1); 
                    xc(nna,l+1)=(xbp(nb,1)+xap(nna,1))/2;
                    yc(nna,l+1)=(ybp(nb,1)+yap(nna,1))/2;
                    zc(nna,l+1)=(zbp(nb,1)+zap(nna,1))/2; 
                    xb(nb,l+1)=0;
                    yb(nb,l+1)=0;
                    zb(nb,l+1)=0;
                    xa(nna,l+1)=0;
                    ya(nna,l+1)=0;
                    za(nna,l+1)=0;
                    Ar(Cellb,1)=Ar(Cellb,1)-1;
                    Br(Cellb,1)=Br(Cellb,1)-1;                   
                    ConcAdt(Cellalp(nna,1),1)=(nacelldt(Cellalp(nna,1),1)-1)/(AvNumb*deltax*deltay*deltaz);
                    ConcBdt(Cellblp(nb,1),1)=(nbcelldt(Cellblp(nb,1),1)-1)/(AvNumb*deltax*deltay*deltaz); 
                    nacelldt(Cellalp(nna,1),1)=nacelldt(Cellalp(nna,1),1)-1;
                    nbcelldt(Cellblp(nb,1),1)=nbcelldt(Cellblp(nb,1),1)-1;
                    else
                    continue
                end
               else
                    continue
                end
         else
             continue
         end
   end
   end

   t=t+deltat(l,1);
   t2=t2+Das(l,1);
   t2p=t2p+Dass(l,1);   
   tt=[tt;t];
   dd(l+1,1)=nnz(find(xa(find(xa(:,l+1)>0),l+1)<L));
   ddt(l+1,1)=nnz(find(xa(find(xa(:,l+1)>0),l+1)<L/2));
   ddtc(l+1,1)=nnz(find(xc(find(xc(:,l+1)>0),l+1)<L/2));

end
   Das(l+1,1)=0;  

tdiff=min([(L/10)^2/(2*DA*3),(L/10)^2/(2*DB*3),(L/10)^2/(2*DC*3)]);

jj=1;
for ind=1:floor(0.01/(jj*tdiff))+1
    vv=find(jj*(ind-1)*tdiff<=tt);
    vvp=find(tt<jj*ind*tdiff); 
    int=intersect(vv,vvp);
    if isempty(int)
        vv1=0;
        vv2=0;
        vv3=0;
        be=1;
    else
        vv1=ddt(int,1);
        vv2=ddtc(int,1);
        vv3=Das(int,1);
        be=length(int);
    end
    ww1(ind,I)=sum(vv1)/be;
    ww2(ind,I)=sum(vv2)/be;
    ww3(ind,I)=sum(vv3)/be;
    ttf(ind,1)=jj*(ind-1)*tdiff+(jj*ind*tdiff-jj*(ind-1)*tdiff)/2;
end

ss=floor(0.01/(jj*tdiff))+1;
jj=10;
t11=ttf(ind,1);
t22=tmax;

for ind=ss+1:ss+floor(t22/(jj*tdiff))+1
    vv=find(jj*(ind-1-ss)*tdiff+t11<=tt);
    vvp=find(tt<jj*(ind-ss)*tdiff+t11);   
    int=intersect(vv,vvp);
    if isempty(int)
        vv1=0;
        vv2=0;
        vv3=0;
        be=1;
    else
        vv1=ddt(int,1);
        vv2=ddtc(int,1);
        vv3=Das(int,1);
        be=length(int);
    end
    ww1(ind,I)=sum(vv1)/be;
    ww2(ind,I)=sum(vv2)/be;
    ww3(ind,I)=sum(vv3)/be;
    ttf(ind,1)=jj*(ind-1-ss)*tdiff+t11+(jj*(ind-ss)*tdiff-jj*(ind-1-ss)*tdiff)/2;
end

for bb=1:length(ww1(:,I))
    if sum(ww1(bb,:))==0
        ww1p(bb,1)=0;
    else
        bb2=length(nonzeros(ww1(bb,:)));
        ww1p(bb,1)=sum(nonzeros(ww1(bb,:)))/bb2;
    end
    if sum(ww2(bb,:))==0
        ww2p(bb,1)=0;
    else
        bb3=length(nonzeros(ww2(bb,:)));
        ww2p(bb,1)=sum(nonzeros(ww2(bb,:)))/bb3;
    end
    if sum(ww3(bb,:))==0
        ww3p(bb,1)=0;
    else
        bb4=length(nonzeros(ww3(bb,:)));
        ww3p(bb,1)=sum(nonzeros(ww3(bb,:)))/bb4;
    end    
end

end

timeAA=[0;ttf(ww1p~=0)];
timeCC=[0;ttf(ww2p~=0)];
AA=[NA0;ww1p(ww1p~=0)];
CC=[0;ww2p(ww2p~=0)];
timeDAS=[ttf(ww3p~=0.1)];
DAS=[ww3p(ww3p~=0.1)];
timeDAS=[timeDAS(DAS~=0)];
DAS=[DAS(DAS~=0)];


figure(1)
plot(timeCC,CC,'r-.','Linewidth',1.5)
axis([0 tmax 0 300])
hold on
legend({'Our Method'},'Location','southeast','NumColumns',1)

xlabel(['Time (s)'])
ylabel(['# Molecules of C'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
    print('Our method', '-dpdf', '-r300');


figure(2)
plot(timeDAS,DAS,'b')
hold on
xlabel(['Time (s)'])
ylabel(['Da_{\Delta s_k}'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
    print('DA', '-dpdf', '-r300');

