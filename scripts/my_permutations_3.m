function [fh, h, p_05, p_95, p_real, pseudo_z_score, vd_mean, vd_std, p_05_corr] = my_permutations_3(mat1,mat2, subs_for_this,xlims,ylims,centroids,centerit,print_n,NPERMS)


% mat2 is going to be permutified

% mat1 will remains constant.


% mat1 should be state information for each subj - across time. 0 == no
% state X (previously selected), 1 == state X, accordign to viterbi path..

% mat2 should be ephyus information, 0 == normal, 1 == it's 20/10 (?)
% percent high - or low. Determined by quantile function earlier on.
% mat2 cna ALSO be annotation information (!) - to be done a bit later
% on...

% rows are always subjects -- cols are time -- or permutations...

% NPERMS = 5000.
perms=NPERMS;


% these functions calculate dice and randomly shuffle shit:
% mydice = @(v1, v2) 2*sum(v1.*v2)/(sum(v1)+sum(v2));

myshuffle = @(v) v(randperm(numel(v)));


unshuffled_dices=[];
for j=1:size(mat1,1)
    unshuffled_dices(end+1) = my_customized_dice(mat1(j,:),mat2(j,:));
end



% DO ALL PERMS:
dicemat=[];

for i=1:perms
    
    
    this_dices=[];
    
    for j=1:size(mat1,1)
        
        this_dices(j) = my_customized_dice(mat1(j,:),myshuffle(mat2(j,:)));
    end
    
    dicemat(end+1,:)=this_dices;
end

 %keyboard;
% OPTIONAL PER SUBJECT, a CHECK:
PERSUB=0;
if PERSUB
    % PER SUBJECT:
    for j=1:size(dicemat,2)
        EDGES = unique(dicemat(:,j));
        EDGES = [EDGES;max(EDGES)+mean(diff(EDGES))];
        EDGES = [min(EDGES)-mean(diff(EDGES));EDGES];
        figure;
        histogram(dicemat(:,j),EDGES)
        
        % so what is the 95% interval?
        p_05=quantile(dicemat(:,j),0.05);
        p_95=quantile(dicemat(:,j),0.95);
        p_99=quantile(dicemat(:,j),0.99);
        p_01=quantile(dicemat(:,j),0.01);
        p_real=unshuffled_dices(j);
        
        vcolors={[0.5 0.5 0.5],[0.5 0.5 0.5],'r'};
        
        if p_real > p_01 && p_real < p_99
            vcolors{end} = [0.3 0.3 0.3];
        end
        
        v = [p_01, p_99, p_real];
        lwidths=[1 1 2];
        for i=1:numel(v)
            line([1 1]*v(i),get(gca,'ylim'),'color',vcolors{i},'linewidth',lwidths(i));
        end
        title(subs_for_this(j));
        
    end
end


% RUN THE AVERAGE OVER SUBJECTS:
vd=mean(dicemat,2);

% keyboard;


% so what is the 95% interval?
p_95=quantile(vd,0.95);
p_99=quantile(vd,0.99);
p_01=quantile(vd,0.01);
p_05=quantile(vd,0.05);
p_05_corr = quantile(vd,0.000833);

p_real=mean(unshuffled_dices);


% so we need to figure out how far p_real is from the mean, in units of
% std...

vd_std = std(vd);
vd_mean = mean(vd);

pseudo_z_score = (p_real - vd_mean) / vd_std;

pseudo_z_thresh = 3.1;

thresh_for_p_real = vd_mean + vd_std * pseudo_z_thresh;




fh=figure;
%EDGES = 0:0.01:0.4;
%h=histogram(vd,EDGES);

if centerit
    centerofmass = median(vd);
    centroids = centroids - mean(centroids) + centerofmass;
    xlims = xlims-mean(xlims) + centerofmass;
end

h=histogram(vd,centroids);

% keyboard;
vcolors={[0.5 0.5 0.5],[0.5 0.5 0.5],'r'};
% if p_real > p_01 && p_real < p_99
if p_real < p_99
    vcolors{end} = [0.3 0.3 0.3];
end
v = [p_01, p_99, p_real];
lwidths=[1 1 2];


% here, we test whether p_real is bigger than xlims(2) or smaller than
% xlims(1)
if p_real > xlims(2)
    xdiff = p_real-xlims(2);
    xlims = xlims + xdiff + diff(xlims)*0.05; % add a little bit to that, too.
    warning('changed the xlims - but keeping the interval');
end

% here, we test whether p_real is bigger than xlims(2) or smaller than
% xlims(1)
if p_real < xlims(1)
    xdiff = xlims(1) - p_real;
    xlims = xlims - xdiff - diff(xlims)*0.05; % add a little bit to that, too.
    warning('changed the xlims - but keeping the interval');
end





set(gca,'xlim',xlims,'ylim',ylims)



for i=1:numel(v)
    line([1 1]*v(i),get(gca,'ylim'),'color',vcolors{i},'linewidth',lwidths(i));
end
if print_n
    text(max(get(gca,'xlim')),max(get(gca,'ylim')),sprintf('N=%d',numel(subs_for_this)),'horizontalalignment','right','verticalalignment','top');
end







    
    
    
    



