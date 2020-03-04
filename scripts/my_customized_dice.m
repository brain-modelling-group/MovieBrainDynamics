function out = my_customized_dice(v1, v2)

%
%
% the use of this function to your own risk.
% we move v1 around a bit, and reclaculate a 'dice' index repeatedly
% then take the max of that.
%
%



td = @(v1, v2) sum(v1.*v2)/min(sum(v1)+sum(v2));

% move v1 around and calculate a couple of dices
dislocs=[3 2 1 0 -1 -2 -3];


v1=reshape(v1,1,numel(v1));
v2=reshape(v2,1,numel(v2));
tmp=[zeros(1,(numel(dislocs)-1)/2) v1 zeros(1,(numel(dislocs)-1)/2)];

dices=[];
check_to_do=[];
for i=1:numel(dislocs)
    to_do = tmp((numel(dislocs)-1)/2+1+dislocs(i) : (numel(dislocs)-1)/2+numel(v1)+dislocs(i));
    check_to_do=[check_to_do; to_do];
    dices(end+1) = td(to_do, v2);
end
    
out=max(dices);

