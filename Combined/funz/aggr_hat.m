% aggregator: computes (avg) demand at each time interval
function [demand] = aggr_hat(ev)

x = [ev.x_hat];

demand = sum(x,2);

end