function g=G(sample,probtype,UOriginp,Unewinp,GOriginp,Gnewinp,PRJfile,HDFfile)
if probtype==1

    Rainfall=sample(1);
    ManningsC=sample(2);
    
    alpha1=0;
    alpha2=1;
    alpha3=sample(3);
    alpha4=sample(4);
    
    while true
        [x,y,~,~]=pearcdf(alpha3,alpha4);
        y=y.*alpha2+alpha1;
        if sum(isnan(y))==0 && sum(find(y==inf))==0 && sum(y)~=0 && sum(find(y==0))<5
            break
        else
            alpha4=alpha4+0.05;
        end
    end

    if min(x)<-3
        minx=-3;
    else
        minx=min(x);
    end

    if max(x)>3
        maxx=3;
    else
        maxx=max(x);
    end

    idx=find(x>-3 & x<3);
    x=x(idx);
    y=y(idx);
    y=100.*y./sum(y);
    x=x.*24./(max(x)-min(x));
    x=x-min(x);

    R=zeros(24,1);
    for i=1:24
        dx=abs(x-i);
        idx=find(dx==min(dx));
        R(i)=y(idx(1));
    end

    R=R./sum(R).*100;

    Q=(Rainfall*0.01).*R;
    U=[0 25.849 60.314 77.031 56.395 20.759 12.123 0 0]';

    C=zeros(33,24);
    idx=1;
    for i=2:33-length(U)+1
        C(i:i+length(U)-1,idx)=U.*(Q(i-1).*0.1);
        idx=idx+1;
    end

    %Zhushei,Dashibeng,Baoli1,Houkei,Baoli2
    F=[0.02652 0.0472 0.01733 0.02912 0.02253];
    FinalQ=sum(C,2).*F;
    FinalQ=vertcat(FinalQ,zeros(16,length(F)));

    [size1,size2]=size(FinalQ);
    for i=1:size2
        dummy='';
        for j=1:size1
            dummy=append(dummy,sprintf('%8.3f',FinalQ(j,i)));
        end
        strF1{i}=dummy;
    end

    idx=1;
    for i=1:size2
        for j=1:5
            st=(j-1)*80+1;en=j*80;
            if j==5
                en=length(strF1{i});
            end
            strF2{idx,1}=strF1{i}(st:en);
            idx=idx+1;
        end
    end

    %% Create new input file

    Originp=UOriginp;
    newinp=Unewinp;
    readFid=fopen(Originp,'r');
    readinp=char(fread(readFid, 'uchar')');
    fclose(readFid);
    for i=1:5
        readinp=strrep(readinp,sprintf('row%dzhushei',i),strF2{i});
    end
    for i=1:5
        readinp=strrep(readinp,sprintf('row%ddashibeng',i),strF2{i+5});
    end
    for i=1:5
        readinp=strrep(readinp,sprintf('row%dbaoli1',i),strF2{i+10});
    end
    for i=1:5
        readinp=strrep(readinp,sprintf('row%dhoukei',i),strF2{i+15});
    end
    for i=1:5
        readinp=strrep(readinp,sprintf('row%dbaoli2',i),strF2{i+20});
    end
    writeFid=fopen(newinp,'w');
    fprintf(writeFid,'%s' ,readinp);
    fclose(writeFid);

    Originp=GOriginp;
    newinp=Gnewinp;
    readFid=fopen(Originp,'r');
    readinp=char(fread(readFid, 'uchar')');
    fclose(readFid);
    readinp=strrep(readinp,'ManningsC',num2str(ManningsC));
    writeFid=fopen(newinp,'w');
    fprintf(writeFid,'%s' ,readinp);
    fclose(writeFid);

    %% Analysis
    h=actxserver('RAS631.HECRASCONTROLLER');
    h.Project_Open(PRJfile)
    h.ShowRas
    h.Compute_CurrentPlan(0,0)
    while h.Compute_Complete()==0
        pause(60)
    end
    h.Project_Save
    h.QuitRas
    !taskkill /im ras.exe

    %% Extract Result
    cellnum=19903;
    data1=h5read(HDFfile,'/Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/Baoli Basin/Water Surface');
    data1=double(data1);
    ref1=data1(cellnum,1);
    ref2=max(data1(cellnum,2:end));
    g=3-(ref2-ref1);
end
end