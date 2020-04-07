function beam=fcm_beam2varbeam(beam)

if size(beam.s{1},3)>1
    for s=1:length(beam.s)
        beam.s{s}=squeeze(var(beam.s{s},0,2));
    end
    beam.timewindow=beam.bands;
    beam.timepts=mean(beam.bands,2);
    beam.bands=[beam.bands(1,1) beam.bands(end,2)];
else
    for s=1:length(beam.s)
        beam.s{s}=var(beam.s{s},0,2);
    end
    beam.timewindow=[beam.timewindow(1,1) beam.timewindow(end,2)];
    beam.timepts=mean(beam.timewindow,2);
end
if isfield(beam,'connectionfile')
    beam=rmfield(beam,'connectionfile');
end
