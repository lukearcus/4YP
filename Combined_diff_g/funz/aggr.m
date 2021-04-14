% aggregator: computes (avg) demand at each time interval
function [demand] = aggr(ev)

x = [ev.schedule];

demand = sum(x,2);

end