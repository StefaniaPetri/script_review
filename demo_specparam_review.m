%% demo script for specparam algorhythm 
% adapted from fieltrip site:
% https://www.fieldtriptoolbox.org/example/spectral/fooof/ 

% see also:
% https://neuroimage.usc.edu/brainstorm/Tutorials/Fooof
% https://fooof-tools.github.io/fooof/


% add path to fieldtrip

clear all
close all

col_list = {'#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2','#7f7f7f', '#bcbd22', '#17becf'};
for nn =1:length(col_list)
    current_col = col_list{nn};
plot(1:10,(1:10)+(10*nn),"MarkerFaceColor",  current_col)
hold on
legend
end

% simulate data
F = 1; % weight of (F)ractal components of the simulated data
O = 1; % weight of (O)scillatory components of the simulated data
t = (1:60000)/1000; % time axis
for rpt = 1:100
    % use a simple method to make pink noise that does not rely on the digital signal processing toolbox
    fn = cumsum(randn(1,length(t)));
    fn = fn./max(abs(fn));

    sgn10 = ft_preproc_bandpassfilter(randn(1,length(t)),1000,[8 12],[],'firws');
    sgn10 = 0.15.*sgn10./max(abs(sgn10));

    sgn30 = ft_preproc_bandpassfilter(randn(1,length(t)),1000,[25 35],[],'firws');
    sgn30 = 0.05.*sgn30./max(abs(sgn30));

    % add a 10 Hz oscillation
    data.trial{1,rpt} = F * fn + O * sgn10 + O * sgn30;
    data.time{1,rpt}  = t;
    data.label{1}     = 'chan';
    data.trialinfo(rpt,1) = rpt;
end

% chunk into 2-second segments
cfg               = [];
cfg.length        = 2;
cfg.overlap       = 0.5;
data              = ft_redefinetrial(cfg, data);

% compute the fractal and original spectra
cfg               = [];
cfg.foilim        = [1 45];
cfg.pad           = 2;
cfg.tapsmofrq     = 2;
cfg.method        = 'mtmfft';
cfg.output        = 'fooof_aperiodic';
fractal = ft_freqanalysis(cfg, data);
cfg.output        = 'pow';
original = ft_freqanalysis(cfg, data);

% subtract the fractal component from the power spectrum
cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2-x1';
oscillatory = ft_math(cfg, fractal, original);

% original implementation by Donoghue et al. 2020 operates through the semilog-power
% (linear frequency, log10-power) space and transformed back into linear-linear space.
% thus defining an alternative expression for the oscillatory component as the quotient of
% the power spectrum and the fractal component
cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2./x1';  % equivalent to 10^(log10(x2)-log10(x1))
oscillatory_alt = ft_math(cfg, fractal, original);

% display the spectra on a log-log scale
fig1=figure("Color","white");
fig1.Position(3:4) = [fig1.Position(3)*1.5,fig1.Position(4)*1.5];
 hold on;
lwd = 3;
plot(log(original.freq), log(original.powspctrm),"Color",  '#000000','LineWidth',lwd);
plot(log(fractal.freq), log(fractal.powspctrm),"--","Color",  '#0000ff','LineWidth',lwd);
plot(log(fractal.freq), log(oscillatory.powspctrm),"Color", '#008000','LineWidth',lwd);
xlabel('log(Frequency)'); ylabel('log(Power)'); grid on;
legend({'Original Data','Aperiodic','Oscillatory = Original Data - Aperiodic'},'location','southwest','Box','off');
if F~=0 && O==0
    title('pure fractal signal');
elseif F==0 && O~=0
    title('pure oscillatory signal');
elseif F~=0 && O~=0
    title('mixed signal');
end
fontsize(16,"points")


exportgraphics(fig1,'demo_fooof_review_ft1.png','Resolution',300)

hold off 
fig2=figure("Color","white");
hold on;
plot(log(original.freq), log(original.powspctrm),"Color",  '#000000','LineWidth',lwd);
plot(log(fractal.freq), log(fractal.powspctrm),"--", "Color",  '#0000ff','LineWidth',lwd);
plot(log(oscillatory_alt.freq), log(oscillatory_alt.powspctrm),"Color",  '#008000','LineWidth',lwd);
xlabel('log(Frequency)'); ylabel('log(Power)'); grid on;
legend({'Original Data','Aperiodic','Oscillatory = Original Data / Aperiodic'},'location','southwest','Box','off');
title('Oscillatory = Spectrum / Aperiodic');
fontsize(16,"points")
% saveas(fig, 'myfigure.png', 'png');
fig2.Position(3:4) = [fig2.Position(3)*1.5,fig2.Position(4)*1.5];

exportgraphics(fig2,'demo_fooof_review_ft2.png','Resolution',300)
