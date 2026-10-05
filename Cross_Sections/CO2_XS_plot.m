xs_CO2_Carbon=csvread('/Users/pps18/Documents/Cross_Sections/sig_CICII_fromCO2.txt',1,0)';
xs_CO2_4pg=csvread('/Users/pps18/Documents/Cross_Sections/CO2_4PG_02_XS.txt',1,0);

xs_CO2_4PG=interp1(xs_CO2_4pg(:,1),xs_CO2_4pg(:,2),xs_CO2_Carbon(1,:))*1e-18;
xs_CO2_1657=xs_CO2_Carbon(3,:)-xs_CO2_4PG;
xs_CO2_1657(or(isnan(xs_CO2_1657),xs_CO2_1657<0))=0;
xs_CO2_1657(xs_CO2_Carbon(1,:)<23)=0;

colorset=get(gca,'ColorOrder');

figure;
plot(xs_CO2_Carbon(1,:),xs_CO2_Carbon(3,:),'DisplayName','CO2 Total','LineWidth', 2,'color','k');
hold on
plot(xs_CO2_Carbon(1,:),xs_CO2_1657,'DisplayName','CO2 CI1657','LineWidth', 2,'color',colorset(2,:),'LineStyle','--');
hold on
plot(xs_CO2_Carbon(1,:),xs_CO2_4PG,'DisplayName','CO2 4PG','LineWidth', 2,'color',colorset(2,:),'LineStyle',':');
hold on
% set(gca,{'xscale','yscale'},{'log','log'});
scatter([30,100],[30.9*1.65/23.5,18.6*1.11/17.8]*1e-19,'DisplayName',['CO2 4PG Ajello+19'],...
    'MarkerEdgeColor',colorset(2,:),'MarkerFaceColor',colorset(2,:),'Marker','*',...
    'SizeData',150,'LineWidth',1.2);%,'HandleVisibility','off'); 
hold on
scatter([30,100],[30.9*1.65/23.5+6.75,18.6*1.11/17.8+13.9]*1e-19,'DisplayName',['CO2 Total Ajello+19'],...
    'MarkerEdgeColor','k','MarkerFaceColor','k','Marker','*',...
    'SizeData',150,'LineWidth',1.2);%,'HandleVisibility','off'); 
% plot(xs_CO2_Carbon(1,:),xs_CO2_4PG*1.88/1.11,'DisplayName','CO2 4PG','LineWidth', 2,'color',colorset(2,:),'LineStyle',':')
set(gca, {'XScale', 'YScale'},{'log','log'});
set(gca, {'XLim','Fontsize'},{[8,5e2],18});

xlabel('Energy$/eV$','Interpreter','latex');
ylabel('Cross Section $/cm^{-2}$','Interpreter','latex');
grid on
legend('location','southeast');