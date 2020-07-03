function [cnt, out] = grab_collection_of_ones(v)

% fetch indices and put em into a cell array
% input a logical vectors
% outputs indices indicating the start and count of contiguous_regions

if numel(v) == 0
    out={};
    cnt=[];
    return
end



if v(1)==0 
    out={};
    cnt=[];

elseif v(1)==1
    out={[1]};
    cnt=[1];
end
    


for i=2:numel(v)
    if v(i-1) == 1 && v(i) == 1
        out{end} = [out{end} i];
        cnt(end) = cnt(end) + 1;
    elseif v(i-1) == 1 && v(i) == 0
        
    elseif v(i-1) == 0 && v(i) == 1
        out{end+1} = [];
        cnt(end+1) = 0;
        out{end} = [out{end} i];
        cnt(end) = cnt(end) + 1;
    elseif v(i-1) == 0 && v(i) == 0
        
    end
end