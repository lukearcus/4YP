function [MaxStepsize] = calc_step_size(ev,deltas,coord)

A_0 = diag(coord.A_0);
Lf = 0;
alpha = inf;
for m = 1:size(deltas,2)
   A_m = diag(deltas(m).pr(:,1));
   grad = [];
   for i = 1:size(ev,2)
      gradRow = [];
      for j = 1:size(ev,2)
          if i == j
             gradRow = [gradRow, 2*A_0 + 2*A_m]; 
          else
             gradRow = [gradRow, A_0 + 2*A_m];
          end
      end
      grad = [grad;gradRow];
   end
   eigs = eig(grad);
   Lf = max(Lf,max(eigs));
   alpha = min(alpha,min(eigs));
end

% T = size(ev(1).schedule,1);
% ubs = zeros(T,size(ev,2));
% for i = 1:size(ev,2)
%    ubs(:,i) = ev(i).poly.ub; 
% end
% 
% sigmaub = sum(ubs,2);
% 
% g = zeros(1,size(deltas,2));
% for i = 1:size(deltas,2)
%    g(i) = sigmaub'*(deltas(i).pr(:,1).*sigmaub + deltas(i).pr(:,2));
% end
% 
% [~,m] = max(g);
% 
% UncertainGrad = (2*deltas(m).pr(:,1).*sigmaub+ deltas(m).pr(:,2))/size(ev,2);
% CommonGrad = UncertainGrad + coord.A_0.*sigmaub;
% 
% Fn = zeros(T*size(ev,2),1);
% for i = 1:size(ev,2)
%     Fn((((i-1)*T)+1):i*T) = CommonGrad + coord.A_0.*ubs(:,i);
% end
% S = svd(Fn);
% Lf = max(S);
% alpha = min(S);

MaxStepsize = (-Lf^2+sqrt(Lf^4+4*alpha^2))/(2*alpha);

end