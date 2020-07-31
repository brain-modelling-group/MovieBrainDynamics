% we make the new analysis


s=load('sd.mat');


addpath('hmm_new');




% so clean up the following code a little bit -- then allow selection of N
% states - then allow determination of task, too.

cmat=[[230, 25, 75];...
    [60, 180, 75];...
    [255, 225, 25];...
    [0, 130, 200];...
    [245, 130, 48];...
    [70, 240, 240];...
    [240, 50, 230];...
    [250, 190, 190];...
    [0, 128, 128];...
    [230, 190, 255];...
    [170, 110, 40];...
    [255, 250, 200];...
    [128, 0, 0];...
    [170, 255, 195];...
    [0, 0, 128];...
    [128, 128, 128];...
    [255, 255, 255];...
    [0, 0, 0]];


%
%
% all this to extract subs_to_use -- which subject we will use in this
%
%
%

% analysis (i.e., 'all').
% deal with the subjects...
slist=struct();
slist.mov1a = [1:21];
slist.mov2a = [1:16 18:21];
slist.resta = [1:21];
slist.mov1b = [2:5 7:14 16:20];
slist.restb = [2:5 7:14 16:20];

% incorporate bad subjects -- that are OVERALL bad:
bads = [4 13 15];
fns=fieldnames(slist);
for i=1:numel(fns)
    for j=1:numel(bads)
        slist.(fns{i})(slist.(fns{i})==bads(j))=[];
    end
end

% then incorporate the occasioal bads:
slist.mov2a(slist.mov2a==11)=[];
slist.mov2a=sort([slist.mov2a 15]);
slist.resta(slist.resta==5)=[];
slist.resta(slist.restb==6)=[];


analyses = {...
    {'mov1a'},'mov1a'; ... % only mov1
    {'mov1b'},'mov1b'; ... % only mov2 day2
    {'mov1a','mov1b'},'mov1a-1b'; ... % mov1 day and and day 2
    {'resta'},'resta'; ...  % only resta
    {'restb'},'restb'; ...  % only rest day 2
    {'resta','restb'},'resta-b'; ...  % rest a and rest day 2
    {'mov1a','mov1b','resta','restb'},'all'; ...  % all movie and rest day 1 aand day 2 but not repeats
    };




load valid_inferences_all.mat
% pick any inference you wish to check here
this_inference=valid_inferences(1);


scan_names = analyses{7, 1};

analysis = analyses{7, 2};
% scan_names = {'mov1a','mov1b'};
% scan_names= {'mov1a'};

subs_to_use = [1:21]; % magic to make sure that IF you specified more than one functional task,
% subjects that do NOT have done them ALL are excluded.Otherwise
% hmm gets confused, I think.

nsubs=21;
for i=1:numel(scan_names)
    subs_to_use = intersect(subs_to_use, slist.(scan_names{i}));
    nsubs=numel(subs_to_use);
end


% grab the pupil Diameter
eph = load('../data/physiological/PD/s.mat');
eph_part = eph.s;



% EYE
% we go thourgh the states
all_pnts={};
t_pnts=[];
p_pnts=[];
npt_pnts=[];


my_s=s.sd(1).aroma.all.(sprintf('sm%d',this_inference));

for iState=1:10
    
    pnts=[];
    sd_pnts=[];
    n_pnts=[];
    % we go thorugh the subjects...
    for iSub=1:numel(subs_to_use)
        sub=subs_to_use(iSub);
        
        part=my_s.(sprintf('s%d',sub));
        
        % check mov1a
        hrdata=[];
        hmmdata=[];
%         if isfield(part.mov1a,'eyedata')
%             hrdata=[hrdata; part.mov1a.eyedata];
%             hmmdata=[hmmdata; part.mov1a.hmmmat(:, iState)];
%         end
        
        % the variable name is hrdata; but it is PD data.
        if isfield(eph_part.mov1b,sprintf('s%.3d',iSub))
            hrdata=[hrdata; eph_part.mov1b.(sprintf('s%.3d',iSub))(6:end)];
            hmmdata=[hmmdata; part.mov1b.hmmmat(:, iState)];
        end
        
        if sum(hmmdata) > 0 && numel(hrdata)>0
            hrdata=hrdata-mean(hrdata);
            pnts(end+1) = mean(hrdata(hmmdata==1));
            sd_pnts(end+1) = std(hrdata(hmmdata==1));
            % n_pnts(end+1) = numel(pnts);
        end
        
    end
    
    
    npt_pnts(end+1)=numel(pnts);
    all_pnts{end+1} = pnts;
    [H,P,CI,STATS] = ttest(pnts);
    % [P2,CI2,STATS2] = kruskalwallis(pnts);
    t_pnts(end+1)=H;
    p_pnts(end+1)=P;
    % npt_pnts(end+1) = P2;
end
%
%
y = num2cell(1:numel(all_pnts));
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], all_pnts, y, 'UniformOutput', 0); % adding labels to the cells
X = vertcat(x{:});
%
figure;
boxplot(X(:,1), X(:,2));
%
%
% anno_mat = load('/home/johan/mnt/hpcworking/projects/hmm-movie/rawdata/1_Data/Annotation/annobigmat.mat');
% am=anno_mat.bigmat(6:end,:);
% bigmat_labels={'language','changepoint','plottwist','faces   +','n','-','scenes   +','n','-'};
%
% am2=am(:, [4 6 7 9 1 2]);
% al2=bigmat_labels([4 6 7 9 1 2]);

subs_to_use=[2  3    7     8     9    10    11    12    14    16    17 18    19    20];

state_consistency = {};
for iState=1:10
    
    c=[];
    for iSub=1:numel(subs_to_use)
        sub=subs_to_use(iSub);
        part=my_s.(sprintf('s%d',sub));
        c=[c part.mov1a.hmmmat(:, iState)];
    end
    inds_to_check1 = sum(c')>=13;
    
    
    c=[];
    for iSub=1:numel(subs_to_use)
        sub=subs_to_use(iSub);
        part=my_s.(sprintf('s%d',sub));
        c=[c part.mov1b.hmmmat(:, iState)];
    end
    inds_to_check2 = sum(c')>=13;
    
    
    if sum(inds_to_check1)>0 || sum(inds_to_check2)>0
        state_consistency{iState} = logical(inds_to_check1) | logical(inds_to_check2);
    end
    
end


% check if things are different where there is a highly salient episode in
% the movie



meanvals=[];
for i=1:10
    meanvals(end+1) = mean(all_pnts{i});
end

s=sprintf('Deviation from Mean Value during brain state visits, values\n');
s=[s sprintf('%.4f\t',meanvals)];
s=[s sprintf('\n')];

s=[s sprintf('Deviation from Mean Value during brain state visits, p-values\n')];
s=[s sprintf('%.4f\t',p_pnts)];
s=[s sprintf('\n')];
s=[s sprintf('\n')];

fprintf(s);

