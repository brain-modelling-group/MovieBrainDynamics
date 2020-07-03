
% containing:
% 'mean_em',
% 'prob',
% 'FO',
% 'dign',
% 'avg_life'

ANALYSIS = 'all';


d=dir(['../results_10/aroma/' ANALYSIS '/Summary_measures_rep_*.mat']);
summary_vars={};
for i=1:numel(d) 
    summary_vars{end+1} = load([d(i).folder filesep d(i).name]);
end
    

% containing
% 'hmm',
% 'T',
% 'J',
% 'K',
% 'n_sub',
% 'fe'
d=dir(['../results_10/aroma/' ANALYSIS '/HMMrun_rep_*.mat']);
mark=[];
for i=1:numel(d)
    if numel(regexp(d(i).name,'data','match'))>0
        mark(end+1) = i;
    end
end
d(mark)=[];

        
hmm_vars={};
for i=1:numel(d) 
    hmm_vars{end+1} = load([d(i).folder filesep d(i).name]);
end

% grab the numbers in the files, since reading with a '*' will NOT sort the
% numbers.
real_nums = []; for i=1:numel(d); v=regexp({d(i).name},'[0-9]*','match');v=v{1}{1};real_nums(end+1) = str2double(v); end

disp('state + count; for each of the 10 states. This is followed by the i of the hmm run, and a number saying how many of the 10 states are actually expressed');
valid_inferences = [];
for i=1:numel(hmm_vars)
    this_path = hmm_vars{i}.vpath;
    s='\n';
    for j=1:10
        
        if any(j==this_path)
            s = [s '\t' sprintf('%d %d',j,sum(j==this_path))];
        else
            s = [s '\t'];
        end
    end
    
    fprintf([s sprintf('\t%d\t%d\n', real_nums(i), numel(unique(this_path)))]);
    if numel(unique(this_path)) == 10
        valid_inferences(end+1) = real_nums(i);
    end
    
end

save(['valid_inferences_' ANALYSIS '.mat'],'valid_inferences');