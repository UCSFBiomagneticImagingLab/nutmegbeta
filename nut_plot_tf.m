function nut_plot_tf(pop,roi,stat)


V=stroke_anatroivalue(pop,roi);
VV=shiftdim(mean(V,1));

% if nargin>2,
%     VV(stat.diff.p_corr>.05)=0;
% end

tfh = zeros(size(VV));

bands=pop.bands;
bands(bands(:,1)<45,1)=bands(bands(:,1)<45,1)-.5;
bands(bands(:,2)<45,2)=bands(bands(:,2)<45,2)+.5;
bands(bands<1)=1;
    
for timebin=1:size(pop.timewindow,1)
    for freqbin=1:size(bands,1)
        tfh(timebin,freqbin) = patch([pop.timewindow(timebin,1) pop.timewindow(timebin,1) pop.timewindow(timebin,2) pop.timewindow(timebin,2)],[bands(freqbin,[1 2]) bands(freqbin,[2 1])],VV(timebin,freqbin),'LineStyle','none');
    end
end

if nargin>2
    [t,f]=ind2sub(size(VV),find(stat.diff.p_corr<.05));
    patch([pop.timewindow(min(t),1) pop.timewindow(min(t),1) pop.timewindow(max(t),2) pop.timewindow(max(t),2)],[bands(min(f),[1 2]) bands(max(f),[2 1])],'w','facecolor','none','edgecolor','r','linewidth',3);
end

axis(gca,[min(pop.timewindow(:)) max(pop.timewindow(:)) 0 max(pop.bands(:))]);
tmp=bands(:);
tmp(abs(tmp-50)<6)=50;
tmp=unique(tmp);
set(gca,'YScale','log','YTick',tmp,'YTickLabel',num2str(tmp),'box','on')

MX=max(abs(VV(:)));
caxis(gca,[-MX MX]);

xlabel('Time (ms)'); 
ylabel('Frequency (Hz)');

cbar = colorbar;
ylabel(cbar,'Power Change (dB)');

