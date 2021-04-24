%Parameter sweep data analysis
load('TD_Parameter_Sweep_4_20.mat')
load('TD_Parameter_Sweep_4_20_FrqMgs.mat')

figure(1)
clf
hold on
for i = 1:length(Mags)
    plot(Freqs,MeanXs(:,i))
    MagLegend{i} = sprintf('%.1f',Mags(i));
end
MagLegend{1} = sprintf('%.2f',Mags(1));
xlabel('Freqency (Hz)')
ylabel('Mean Final X Position')
title('Mean Final X Position vs Wind Change Frequency')
legend(MagLegend)
hold off

figure(2)
clf
hold on
for i = 1:length(Mags)
    plot(Freqs,MeanYs(:,i))
end
xlabel('Freqency (Hz)')
ylabel('Mean Final Y Position')
title('Mean Final Y Position vs Wind Change Frequency')
legend(MagLegend)
hold off

figure(3)
clf
hold on
for i = 1:length(Mags)
    plot(Freqs,StdXs(:,i))
end
xlabel('Freqency (Hz)')
ylabel('Stdev Final X Position')
title('Stdev Final X Position vs Wind Change Frequency')
legend(MagLegend,'Location','northwest')
hold off

figure(4)
clf
hold on
for i = 1:length(Mags)
    plot(Freqs,StdYs(:,i))
end
xlabel('Freqency (Hz)')
ylabel('Stdev Final Y Position')
title('Stdev Final Y Position vs Wind Change Frequency')
legend(MagLegend,'Location','northwest')
hold off

figure(5)
clf
hold on
for i = 1:length(Freqs)
    plot(Mags,MeanXs(i,:))
    FreqLegend{i} = sprintf('%d',Freqs(i));
end
xlabel('Magnitude (m/s)')
ylabel('Mean Final X Position')
title('Mean Final X Position vs Wind Change Magnitude')
legend(FreqLegend)
hold off

figure(6)
clf
hold on
for i = 1:length(Freqs)
    plot(Mags,MeanYs(i,:))
end
xlabel('Magnitude (m/s)')
ylabel('Mean Final Y Position')
title('Mean Final Y Position vs Wind Change Magnitude')
legend(FreqLegend)
hold off

figure(7)
clf
hold on
for i = 1:length(Freqs)
    plot(Mags,StdXs(i,:))
end
xlabel('Magnitude (m/s)')
ylabel('Stdev Final X Position')
title('Stdev Final X Position vs Wind Change Magnitude')
legend(FreqLegend)
hold off

figure(8)
clf
hold on
for i = 1:length(Freqs)
    plot(Mags,MeanXs(i,:))
end
xlabel('Magnitude (m/s)')
ylabel('Stdev Final Y Position')
title('Stdev Final Y Position vs Wind Change Magnitude')
legend(FreqLegend)
hold off