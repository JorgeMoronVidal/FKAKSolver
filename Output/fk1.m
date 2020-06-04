%
% execute C++ code
%

% !fk1

%
% use MATLAB routine to plot output
%

addpath('..');

nvert = 3;

del1 = cell(3,1);
del2 = cell(3,1);
var1 = cell(3,1);
var2 = cell(3,1);
kur1 = cell(3,1);
chk1 = cell(3,1);
cost = cell(3,1);
Eps  = cell(3,1);
mlmc_cost = cell(3,1);
std_cost = cell(3,1);
ls   = cell(3,1);
Nls  = cell(3,1);
L    = cell(3,1);

for o = 1:3
  filename = ['fk1_8'];
  [d1 d2 v1 v2 k1 c1 c2 eps mlmc std l Nl ] = mlmc_read(filename);
  del1{o} = d1;
  del2{o} = d2;
  var1{o} = v1;
  var2{o} = v2;
  kur1{o} = k1;
  chk1{o} = c1;
  cost{o} = c2 ./ (30*4.^(0:(length(c2)-1)));
  Eps{o} = eps;
  mlmc_cost{o} = mlmc;
  std_cost{o} = std;
  ls{o}  = l;
  Nls{o} = Nl;
  L{o} = length(d1) - 1;
end


%
% plot figures
%

close all

nvert = 2;

figs(1) = figure; 
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.75*nvert]; set(gcf,'pos',pos);

set(0,'DefaultAxesColorOrder',[0 0 0]);
set(0,'DefaultAxesLineStyleOrder','-*|--*|-.*|:*')

subplot(nvert,2,1)
plot(0:L{1},log2(var2{1}), ...
     1:L{2},log2(var1{2}(2:end)), ...
     1:L{1},log2(var1{1}(2:end)), ...
     1:L{3},log2(var1{3}(2:end)) )
xlabel('level $\ell$','Interpreter','latex'); ylabel('log_2 variance');
current = axis; axis([ 0 L{1} current(3:4)+[-6 0] ]);
legend({'$P_\ell$ (orig)','$\Delta P_\ell$ (orig)','$\Delta P_\ell$ (new1)', ...
       '$\Delta P_\ell$ (new2)'},'Location','SouthWest','Interpreter','latex')

subplot(nvert,2,2)
plot(0:L{1},log2(abs(del2{1})), ...
     1:L{2},log2(abs(del1{2}(2:end))), ...
     1:L{1},log2(abs(del1{1}(2:end))), ...
     1:L{3},log2(abs(del1{3}(2:end))) )
xlabel('level $\ell$','Interpreter','latex'); ylabel('log_2 |mean|');
current = axis; axis([ 0 L{1} current(3:4)+[-4 0] ]);
legend({'$P_\ell$ (orig)','$\Delta P_\ell$ (orig)','$\Delta P_\ell$ (new1)', ...
       '$\Delta P_\ell$ (new2)'},'Location','SouthWest','Interpreter','latex')

set(0,'DefaultAxesLineStyleOrder','--*|-.*|:*')

subplot(nvert,2,3)
plot(0:L{2},cost{2},0:L{1},cost{1},0:L{3},cost{3})
xlabel('level $\ell$','Interpreter','latex'); ylabel('normalised average cost');
current = axis; axis([ 0 L{1} current(3:4)+[0 .05] ]);
legend({'$\Delta P_\ell$ (orig)','$\Delta P_\ell$ (new1)','$\Delta P_\ell$ (new2)'}, ...
       'Location','NorthEast','Interpreter','latex')

subplot(nvert,2,4)
%semilogy(1:L{2},kur1{2}(2:end),1:L{1},kur1{1}(2:end),1:L{3},kur1{3}(2:end))
plot(1:L{2},kur1{2}(2:end),1:L{1},kur1{1}(2:end),1:L{3},kur1{3}(2:end))
current = axis; axis([ 0 L{1} current(3:4) ]);
xlabel('level $\ell$','Interpreter','latex'); ylabel('kurtosis');
legend({'$\Delta P_\ell$ (orig)','$\Delta P_\ell$ (new1)','$\Delta P_\ell$ (new2)'}, ...
       'Location','NorthWest','Interpreter','latex')


figs(2) = figure;
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.75]; set(gcf,'pos',pos);

set(0,'DefaultAxesLineStyleOrder','-.o|-.x|-.d|-.*|-.s|:o|:x|:d|:*|:s');

subplot(1,2,1)
semilogy(ls{1}, Nls{1}, ls{3}, Nls{3})
xlabel('level $\ell$','Interpreter','latex'); 
ylabel('$N_\ell$','Interpreter','latex'); 
%current = axis; axis([ 0 size(Nls{1},1)-1 current(3:4) ]);
axis([0 10 1 1e7])
for i=1:length(Eps{1})
  labels{i} = num2str(Eps{1}(i));
end
legend(labels,'Location','NorthEast')

set(0,'DefaultAxesLineStyleOrder','-*|--*|-.*|:*')

subplot(1,2,2)
loglog(Eps{2},Eps{2}.^2.*std_cost{2}(:)',...
       Eps{2},Eps{2}.^2.*mlmc_cost{2}(:)',...
       Eps{1},Eps{1}.^2.*mlmc_cost{1}(:)',...
       Eps{3},Eps{3}.^2.*mlmc_cost{3}(:)')
xlabel('accuracy $\varepsilon$','Interpreter','latex');
 ylabel('$\varepsilon^2$ Cost','Interpreter','latex');
%current = axis; axis([ Eps{1}(1) Eps{1}(end) 5*current(3:4) ]);
axis([5e-4 1e-2 10 1e7])
legend('Std MC','MLMC (orig)','MLMC (new1)','MLMC (new2)')


figure(1)
print('-deps2c','fk1_1.eps')
figure(2)
print('-deps2c','fk1_2.eps')
