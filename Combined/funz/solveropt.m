
opt_setup = optimoptions('quadprog');
opt_setup = optimoptions(opt_setup,'Display','none');

% QP solver options for maxnashgame.m
opt_setup_max = optimoptions('quadprog');
opt_setup_max = optimoptions(opt_setup_max,'Display','none');
