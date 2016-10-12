%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fixing each target in turn, then doing ten fold cross validation
%on the set of compounds, to predict which ones
%are/are not targeting the target in question:
function [myPredictions]=BLM_NII(y,theseIndices2,kCompound,kTarget,alpha)
% set parameters
gammab=1;
delta=1;
lengthKCompound = size(kCompound,1);
myPredictions=zeros(size(y));
% get the predictions of interaction for drugs with respect to each target
for i=1:size(y,1);
    % all drugs' labels with respect to this target
    currentY = y(i,:)';   
    for f=1:size(theseIndices2,1)
        %%%!%%%
        train=setdiff(1:lengthKCompound,theseIndices2(f,1):theseIndices2(f,2));
        test=theseIndices2(f,1):theseIndices2(f,2);
        Yleave=y;
        Yleave(i,test)=0;
        % calculate gamma
        gamma=(sum(sum(Yleave)))/lengthKCompound/gammab;
        % calculate kenel of training drugs based on the above representation (profile) of
        % the left lengthKconpound-1 drugs
        temp=pdist(Yleave').^2;
        temp=squareform(temp);
        Kgip=exp(-temp./gamma);
        
        dC=1./kCompound;
        gamma=0.2;
        kC=exp(-dC./(gamma*median(dC(:))));
       % K=(1-alpha)*Kgip+alpha*kCompound;
        % merge the two similarity matrices
        K=(1-alpha)*Kgip+alpha*kC;
       
        K1 = K(train,train); %similarities between training samples
        K2 = K(test,train);  % similarity between test and training samples
        %training:
        na = length(train);
        if sum(currentY(train) == 1) > 0  
        model =  (K1+delta*eye(na)) \ currentY(train);
        else
            inferredY=(kTarget(i,:))*Yleave;
            inferredY=inferredY(train);
            inferredY=(inferredY-min(inferredY))/(max(inferredY)-min(inferredY));
            model =  (K1+delta*eye(na)) \ inferredY';
        end
        %testing:
        % myPredictions(i,test)=(K2*model)';
         temp=(sum(Yleave(:,test))==0);
        for g=1:length(test)
        if temp(g)==0
            myPredictions(i,test(g))=(K2(g,:)*model)';
        else
            myPredictions(i,test(g))=(kCompound(test(g),train)*model)';
        end
        end

    end
end
end


