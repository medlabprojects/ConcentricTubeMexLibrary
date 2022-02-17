function [ P ] = get_power( beta, Lt, L, k )

curv_start = beta + Lt;
curv_end = beta + L;

pts = [beta; beta + Lt; beta + L];
pts = sort(pts);
ivals = [pts(1:end-1), pts(2:end)];
P = zeros(size(ivals,1),1);
for i=1:length(ivals)
   ival = ivals(i,:);
   L = ival(2) - ival(1);
   k1 = 0; k2 = 0; k3 = 0;
   if (sum(ival)/2 >= curv_start(1) && sum(ival)/2 <= curv_end(1))
       k1 = k(1);
   end
   if (sum(ival)/2 >= curv_start(2) && sum(ival)/2 <= curv_end(2))
       k2 = k(2);
   end
   if (sum(ival)/2 >= curv_start(3) && sum(ival)/2 <= curv_end(3))
       k3 = k(3);
   end
 
   P(i) = L*L*(k1*k2 + k2*k3 + k1*k3);
end

P = sum(P);
end

