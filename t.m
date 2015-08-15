ebsd = createTestSample('Isl2A.txt', 100, 100, 'sin', 'Amplitude', 5, 'Period', 10, 'Noise', 0.5, 'Coeff', [1 1; 1 1; 1 1]);
figure(); plot(ebsd, 'property');
xlabel('���','FontSize',14);
ylabel('���','FontSize',14);



grains0 = calcGrains(ebsd, 'threshold', 5*degree);
grains = calcGrains(grains0(grainSize(grains0) > 5), 'threshold', 2*degree);
figure();
plot(get(grains,'ebsd'),'property',KAM0);
xlabel('���','FontSize',14);
ylabel('���','FontSize',14);