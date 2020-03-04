addpath /home/johan/repos/edf-converter/


% to be run in the proper directory...

d=dir('*.mat');

for i=1:numel(d)
    
    edf1=load(d(i).name);
    edf1=edf1.edf1;

    
    plot(edf1);
    
    set(gcf,'paperunits','centimeters');
    set(gcf,'papersize',[25 25]);
    set(gcf,'paperposition',[0 0 25 25]);
    
    
    saveas(gcf,[d(i).name(1:end-4) '.jpg']);
    
    close(gcf);
    

end

