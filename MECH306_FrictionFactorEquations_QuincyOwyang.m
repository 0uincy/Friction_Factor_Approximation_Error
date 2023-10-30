% Friction Factor Plots of Colebrook Equation vs. Swaimee-Jain & Colebrook

% Input Variables
eDc = [0.00001 0.00005 0.0001 0.0002 0.0004 0.0006 0.0008 0.001...
    0.002 0.004 0.006 0.008 0.01 0.015 0.02 0.03 0.04 0.05];
ReMin = 4000; ReMax = 10^8; n = 50;
Re = logspace(log10(ReMin),log10(ReMax),n);

% Swamee-Jain
SJ = zeros(length(Re),length(eDc));
for i = 1:length(eDc)
    for j = 1:length(Re)
        SJ(j,i) = 1.325.*(log(abs(eDc(i)./3.7+5.74./Re(j).^0.9))).^-2
    end
end

% Haaland
H = zeros(length(Re),length(eDc));
for i = 1:length(eDc)
    for j = 1:length(Re)
        H(j,i) = (-1.8.*log10(abs((eDc(i)./3.7).^1.11+6.9/Re(j))))^-2
    end
end

% Colebrook
C = zeros(length(Re),length(eDc));
syms f
for i = 1:length(eDc)
    for j = 1:length(Re)
        C(j,i) = vpasolve(1/sqrt(f)==-2*log10(eDc(i)/3.7+2.51/Re(j)./f),f)
    end
end

% Friction Factor Plots
lw = 1.7;
figure(1)
loglog(Re,SJ(:,1),'Color','#D95319','LineWidth',lw)
axis([0 Re(n) 0 0.1]); yticks(linspace(0.01,0.1,10))
hold on
loglog(Re,H(:,1),'Color','#4DBEEE','LineWidth',lw)
loglog(Re,C(:,1),'Color','#77AC30','LineWidth',lw)
if eDc(1) < 0.0001
    text(Re(n),C(n,1),['  ' num2str(eDc(1),'%.5f')])
else
    text(Re(n),C(n,1),['  ' num2str(eDc(1))])
end
if (length(eDc) > 1)
    for i = 2:length(eDc)
        loglog(Re,SJ(:,i),'Color','#D95319','LineWidth',lw)
        loglog(Re,H(:,i),'Color','#4DBEEE','LineWidth',lw)
        loglog(Re,C(:,i),'Color','#77AC30','LineWidth',lw)
        if eDc(i) < 0.0001
            text(Re(n),C(n,i),['  ' num2str(eDc(i),'%.5f')])
        else
            text(Re(n),C(n,i),['  ' num2str(eDc(i))])
        end
    end
end
hold off
title('Friction Factor as a Function of Reynolds Number')
xlabel('Re'); ylabel('Friction Factor'); grid on
colororder({'k','k'}); yyaxis right; set(gca,'YTick',[])
yl = ylabel('$\frac{\epsilon}{D}$','Interpreter','latex','FontSize',18);
set(get(gca,'YLabel'),'Rotation',0)
yl.Position(1) = yl.Position(1) * 2;
legend('Swaimee-Jain','Haaland','Colebrook')

% Percent Difference from Colebrook
SJ = SJ'; H = H'; C = C';
pctSJ = zeros(length(eDc),length(Re));
for i = 1:length(eDc)
    for j = 1:length(Re)
        pctSJ(i,j) = abs(SJ(i,j)-C(i,j))./C(i,j).*100;
    end
end

pctH = zeros(length(eDc),length(Re));
for i = 1:length(eDc)
    for j = 1:length(Re)
        pctH(i,j) = abs(H(i,j)-C(i,j))./C(i,j).*100;
    end
end

% Plot Percent Error
figure(2)
pctSJ_surf = surf(Re,eDc,pctSJ,'FaceColor','#D95319');
hold on
pctH_surf = surf(Re,eDc,pctH,'FaceColor','#4DBEEE');
hold off
legend('Swaimee-Jain','Haaland')
camorbit(225,-9); set(gca,'xScale','log')
axis([Re(1) Re(n) eDc(1) eDc(end) 0 40])
title('Percent Difference Between Colebrook & Swaimee-Jain & Haaland')
ylabel('$\frac{\epsilon}{D}$','Interpreter','latex','FontSize',18)
xlabel('Re'); zlabel('Percent Error')