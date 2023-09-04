function Result=AKPSOHHS(probtype,dmin,uxdoe,lxdoe,initsample,beta1,beta2,UOriginp,Unewinp,GOriginp,Gnewinp,PRJfile,HDFfile)
format long
regfunc=@regpoly0; corrfunc=@corrgauss;


%% Initial Parameter

[mu,sigma,nRV,dist]=problem(probtype);

nlevel=ceil(nRV*2);
modelimit=0.3;
nfailc=2;nnn=0;

theta = ones(1,size(dist,1));
lob = ones(1,size(dist,1))*1e-2; upb = ones(1,size(dist,1))*1e2;

[xdoe,G_xdoe,nKs1]=LHS(nRV,1,initsample,uxdoe,lxdoe,probtype,mu,sigma,dist,UOriginp,Unewinp,GOriginp,Gnewinp,PRJfile,HDFfile);FE=nKs1;

dmodel=dacefit(xdoe,G_xdoe,regfunc,corrfunc,theta,lob,upb);
CP=0;expansion=0;levelindex=1;plotfig=0;stopcrit=1;
ite=nRV*5;ecosize=nRV*15;
Ulimit=1.4;REIFlimit=0.05;REIFhist=[];Uhist=[];
recCP=[];recFELS=[];

indexloop=1;
while true
    if isempty(G_xdoe(G_xdoe<0))~=1 && expansion==0
        [beta1,beta2,recCP,CP,FELS]=linesearch2(xdoe,G_xdoe,recCP,expansion,levelindex,nlevel,dmodel);
        recFELS(indexloop)=FELS;
    elseif isempty(G_xdoe(G_xdoe<0))~=1 && expansion==1
        if mod(indexloop,5)==0
            [beta1,beta2,recCP,CP,FELS]=linesearch2(xdoe,G_xdoe,recCP,expansion,levelindex,nlevel,dmodel);
            recFELS(indexloop)=FELS;
        end
    end

    %Check to enter leveling stage
    error=1.5;
    nDATA=4;
    if expansion==0 && length(G_xdoe(G_xdoe<0))>nRV*nfailc && length(recCP)>nDATA && 100*(abs(mean(recCP(end-nDATA:end))-recCP(end)))/recCP(end)<error && 100*(abs(max(recCP(end-nDATA:end))-recCP(end)))/recCP(end)<error && 100*(abs(min(recCP(end-nDATA:end))-recCP(end)))/recCP(end)<error
        expansion=1;modelimit=0.85;finalCP=recCP(end);
    end
    %Check if require go back to CP stage
    if expansion==1 && 100*abs(finalCP-CP)/finalCP>10
        expansion=0;modelimit=0.3;beta1=CP+3;beta2=CP;stopcrit=0;
    end
    %Check to enter next level
    nDATA=max(10,nRV^2-15);
    if stopcrit>1 || (length(Uhist)>nDATA && min(Uhist(end-nDATA:end))>5 && max(REIFhist(end-nDATA:end))<0 && length(G_xdoe(G_xdoe<0))>nRV*nnn)
        levelindex=levelindex+1;stopcrit=0;
    end

    if CP<1
        Ulimit=0.3;REIFlimit=0.45;modelimit=0.85;expansion=1;finalCP=0;
    end

    if rand()<modelimit || length(G_xdoe(G_xdoe<0))<1
        mode=3;
    else
        mode=1;
    end

    dummy1=unifrnd(-1,1,ecosize,nRV);
    eco=dummy1.*(unifrnd(beta2,beta1,ecosize,1)./sqrt(sum(dummy1.^2,2)));
    [bestFitness,bestOrganism]=PSO1(ite,ecosize,eco,beta1,beta2,nRV,xdoe,G_xdoe,mode,dmin,CP,dmodel);

    [y,s2]=predictor(repmat(bestOrganism,2,1),dmodel);
    y=y(1);
    sig=sqrt(s2(1));
    w=2;
    if mode==3
        Uhist=vertcat(Uhist,abs(y)/sig);
        REIFhist=vertcat(REIFhist,(y*(1-2*normcdf(y/sig))+sig*(w-sqrt(2/pi)*exp(-0.5*(y/sig)^2))));
        if (Uhist(end)>Ulimit && expansion==1) || (REIFhist(end)<=REIFlimit && expansion==1)
            stopcrit=stopcrit+1;
        end
    end

    xdoe=vertcat(xdoe,bestOrganism);
    [sample,~,~,~]=summonsample(1,mu,sigma,nRV,dist,bestOrganism);G_xdoe=vertcat(G_xdoe,G(sample,probtype,UOriginp,Unewinp,GOriginp,Gnewinp,PRJfile,HDFfile));FE=FE+1;
    dmodel=dacefit(xdoe,G_xdoe,regfunc,corrfunc,theta,lob,upb);
    if levelindex==nlevel+1
        break
    end
    indexloop=indexloop+1;

    disp([levelindex length(G_xdoe) CP])
end

ypred=[];nMCS=1e5;ntMCS=100;
for i=1:ntMCS
    r1sample=normrnd(0,1,nMCS,nRV);
    r2sample=r1sample;
    Rr2sample=rssq(r2sample,2);
    r2sample(Rr2sample<0.8*CP,:)=[];
    [dummy,~]=predictor(r2sample,dmodel);
    ypred=vertcat(ypred,dummy);
    clear dummy
end
Pf=length(ypred(ypred<0))/(nMCS*ntMCS);

fprintf( 'Failure Probability= %d\n',Pf)
fprintf( 'Number of Function Evaluation= %d\n',FE)
fprintf( 'MPP Value= %d\n',CP)

Result=[mean(Pf) FE CP];
end
