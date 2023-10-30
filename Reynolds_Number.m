% Input Variables
eDc = [10^-6 10^-5 10^-4 10^-3];
ReMin = 4000; ReMax = 10^8;

% Define Array
min = log10(ReMin); max = log10(ReMax);
n = 50; Re = logspace(min,max,n);

% Colebrook
C = zeros(length(Re),length(eDc));
syms f
for i = 1:length(eDc)
   for j = 1:length(Re)
       C(j,i) = vpasolve(1/sqrt(f)==-2*log10(eDc(i)/3.7+2.51/Re(j)./f),f)
   end
end


semilogx(Re,C(:,1))
hold on

if (length(eDc) > 1)
    for k = 2:length(eDc)
        semilogx(Re,C(:,k));
    end
end

title('Friction Factor as a Function of Reynold''s Number')
xlabel('Re'); ylabel('f');
grid on;
legend('10^-^6','10^-^5','0.0001','0.001')
