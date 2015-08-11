id = 3;  sum(CnvrgntTypes{id}.ACT_Record{1,1,1,2,3}.VL.WT.E.spktrain)

figure;
id = 1; 
C = 'bkr';
cnt = 0; 
for f = 1:2
    for a = 1:3
        cnt =cnt+1;
        for t = 1:3            
            subplot(2,3,cnt)
            plot(sum(CnvrgntTypes{t}.ACT_Record{1,1,1,f,a}.VL.WT.E.spktrain), 'color', C(t)); hold on;
            xlim([1000 1500])
            title([ 'Freq:' num2str(f) ',  Amp:' num2str(a)  ', Type:' num2str(t)])
            k = waitforbuttonpress;
        end
    end
end

