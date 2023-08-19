function  del = delays(nbits,dim)
% DESCRIPTION
%  Calculates delay of bits in a que along dimension dim.
%  The que is supposed to be empty when nbits is zero 
% INPUT
%  nbits -- A matrix indicating number of bits in the que
%  dim   -- delay dimension
% OUTPUT
%  del   -- A matrix of the same size as nbits with the number of
%           samples indicating how long the bits have been waiting
%           before the que was emptied.
% TRY
%  delays([0 20 10 0])
%  delays([20 10 0 0; 30 20 10 0],2)

% by Magnus Almgren 080220

%nbits = [0 0 1 1 0 1 1 0];
if isempty(nbits) % this is due to diff over an empty dimension
    del = nbits;
    return
end
setifnotexistorempty('dim',firstnonsing(nbits));
nbitsm = swdims(nbits,[dim 1]); % change to user perspective: move slot to first dimension: Tot_slot*1*1*Tot_user
% there are limits to what one can do with mplus et al
tot_user = numel(nbitsm)/size(nbitsm,1);
for user = 1:tot_user
    biq = nbitsm(:,user)>0; % bits in que 
    cnb  = cumsum(biq);
    % find places where the que has become empty
    dbiq = [diff(biq)<0; true]; %% places where there is 1 are emptied
    %%%% for debug for which has delay >0
     %if(ismember(1,dbiq))
     %plot(nbitsm(:,user));
     %title(['User: ',num2str(user)])
     %end
    %%%%%%%%%%%%%%%%%%%
    val = diff([false; cnb(dbiq)]); %% add 0 in the beginning cnb(dbiq) returns the values in cnb whose index places has 1 in dbiq
    delm(dbiq,user) = val; % use only two dimensions for now , delay for user m
end
 %delm(delm==0)=1; %% IS IT NECESSARY ???
 %plot_delay(delm,35);
 del = swdims(reshape(delm,size(nbitsm)),[dim 1]); % move dim back again
1; % for debug purposes
