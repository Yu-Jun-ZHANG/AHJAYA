function [optValue,bestP,conv]=AHJAYA(pType,Max_NFEs,NP)
    addpath('../');
    [LB ,UB,Dim ] = Parameter(pType);
    index=pType;

    MaMi=(repmat((UB-LB),NP,1));
    MiLB=repmat((LB),NP,1);
    X=MaMi.*rand(NP,Dim)+MiLB;
    NPg_current=NP;
    NPmin=3;
    Ra=0.3;
        
    initial_value=0.7;
    u1=1.8;
    u2=2.0;
    b=0.85;
    k=1;
    xlist(k) = initial_value;    

    for i=1:NP
        fitnessX(i)=ModelFunction(index,X(i,:));
    end
    NFEs=NP;
    chaos = Habird(NFEs);
    t=1;
    while NFEs<Max_NFEs
    [fitnessBestX,~]=min(fitnessX);
    [fitnessWorstX,~]=max(fitnessX);
    [~,sortIndexX]=sort(fitnessX);
    Xbest=X(sortIndexX(1),:);
    Xworst=X(sortIndexX(end),:);
        if fitnessBestX == 0
            W1=1;
        else
            W1=(mean(fitnessX)/fitnessBestX);
        end
        if fitnessWorstX == 0
            W2=1;
        else
            W2=(mean(fitnessX)/fitnessWorstX);
        end
        for i=1:NPg_current               
            rand1=rand;
            rand2=rand;
            for j=1:Dim 
                V(i,j) = X(i,j)+W1*rand1*(Xbest(j)-abs(X(i,j)))-W2*rand2*(Xworst(j)-abs(X(i,j)));
            end

            for j=1:Dim
                if V(i,j)>UB(j) || V(i,j)<LB(j)
                    V(i,j)=LB(j)+rand*(UB(j)-LB(j));
                end
            end

            fitnessV(i)=ModelFunction(index,V(i,:));
            NFEs=NFEs+1;
           
            if fitnessV(i)<fitnessX(i)
                X(i,:)=V(i,:);
                fitnessX(i)=fitnessV(i);
            end
        end
        
        NPg_current=round(((NPmin-NP)/Max_NFEs)*NFEs+NP);
        numP=size(X,1);
        if NPg_current<numP
            [~,indexSortP2]=sort(fitnessX);
            X(indexSortP2(NPg_current+1:numP),:)=[]; 
            fitnessX(indexSortP2(NPg_current+1:numP))=[]; 
        end

        if rand<Ra        
            maxBoun=max(X);
            minBoun=min(X);
            for i=1:NPg_current
                Xnew(i,:)=((maxBoun+minBoun)-X(i,:));
                GOP(i,:)=Xnew(i,:)*xlist(k);
                k=k+1;
            if xlist(k-1)<=0
                xlist(k)=(b*(1-u1.*(xlist(k-1).^2)));
            else
                xlist(k)= (1-u2.*xlist(k-1));
            end
                for j=1:Dim
                    if GOP(i,j)>maxBoun(j) || GOP(i,j)<minBoun(j)
                        GOP(i,j)=minBoun(j)+rand*(maxBoun(j)-minBoun(j));
                    end
                end
                fitnessGOP(i)=ModelFunction(index,GOP(i,:));
                NFEs=NFEs+1;
            end
            fitnessNew=[fitnessX,fitnessGOP];
            NewP=[X;GOP];
            [sortFitnessNew,sortIndex]=sort(fitnessNew);
            X=NewP(sortIndex(1:NPg_current),:);
            fitnessX=sortFitnessNew(1:NPg_current);
        end
        
        [fitnessBestX,recordIndex]=min(fitnessX);
        conv(t)=fitnessBestX;
        t=t+1;
    end
    
        %% 最终结果输出
    endNFEs=NFEs;
    optValue=fitnessBestX;
    bestP=X(recordIndex,:);

end

function [xlist] =Habird(maxIter)
    initial_value=0.7;
    u1=1.8;
    u2=2.0;
    b=0.85;
    x = initial_value;
    xlist = zeros(1,maxIter);
    for i=1:maxIter
        if i==1
            xlist(i)=x;
        else
            if xlist(i-1)<=0
                xlist(i)=b*(1-u1.*(xlist(i-1).^2));
            else
                xlist(i) = 1-u2.*xlist(i-1);
            end
        xlist(i) = abs(xlist(i));
        end
    end
end
