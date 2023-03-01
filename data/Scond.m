hH = [
    1/(5-1)
    1/(10-1)
    1/(15-1)
    1/(20-1)
];
conds = [
    3.6e3
    7.9e3
    1.2e4
    1.7e4
];

figure
plot(...
    hH,conds,'k--o',...
    'LineWidth',2,...
    'MarkerSize',15,...
    'MarkerEdgeColor','k')
title('Conditionnement de Sp avec le rapport h/H');
xlabel('h/H');
ylabel('Cond(Sp)');
