function [cnt, out] = grab_collection_inbetweens(v)
% fetch the count and indices of zeros in a vector of ones

% check if there are at LEAST 2 visits:
[cnt, out] = grab_collection_of_ones(v);

switch numel(cnt)
    case 0
        cnt=0;
        out=0;
    case 1
        cnt=0;
        out=0;
    otherwise
        
        % calculate the *LIMITS* where we need to look:
        begin_look = find(v==1,1,'first');
        end_look = find(v==1,1,'last');
        
        [cnt, out] = grab_collection_of_ones(1-v(begin_look:end_look));

        
end
        
        