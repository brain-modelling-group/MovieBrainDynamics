function encmodel = standard_encoding(X,Y,T,options,binsize)
% Compute maps representing the "encoding" model for each
% time point (or window, see below) in the trial (i.e. a GML). 
% The reported statistic is the explained variance from
% regressing Y on the data at this sensor/voxel.
%
% INPUT
% X: Brain data, (time by regions) or (time by trials by regions)
% Y: Stimulus, (time by q); q is no. of stimulus features OR
%              (no.trials by q), meaning that each trial has a single
%              stimulus value
% T: Length of series
% options: structure with the preprocessing options - see documentation in 
%                       https://github.com/OHBA-analysis/HMM-MAR/wiki
% binsize: the trials will be segmented in windows of 'binsize' number of
%               time points, and one encoding model will be estimated per
%               window. By default, this is 1, such that one encoding model
%               is estimated per time point
%
% OUTPUT
% encmodel:  (windows by regions) maps of activation (explained variances)
% 
% Author: Diego Vidaurre, OHBA, University of Oxford (2018)

options.Nfeatures = 0; 
options.pca = 0;
options.embeddedlags = 0;
options.K = 1; 
[X,Y,T] = preproc4hmm(X,Y,T,options); 

N = length(T); ttrial = T(1); 
p = size(X,2); q = size(Y,2);
X = reshape(X,[ttrial N p]);
Y = reshape(Y,[ttrial N q]);

if any(T~=ttrial), error('All trials must have the same length'); end
if nargin<5, binsize = 1; end

nwin = round(ttrial / binsize);
binsize = floor(ttrial / nwin); 

encmodel = zeros(nwin,p);

for t = 1:nwin
    r = (1:binsize) + (t-1)*binsize;
    for j = 1:p
        Xj = X(r,:,j);
        Xj = zscore(Xj(:));
        Yj = reshape(Y(r,:,:),binsize*N,q);
        beta = (Yj' * Yj) \ (Yj' * Xj);
        res = (Xj - Yj * beta).^2;
        res0 = (Xj.^2);
        encmodel(t,j) = 1 - sum(res(:))/sum(res0(:));
    end
end
end