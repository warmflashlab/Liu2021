function [expr_mean,expr_err]=get_plot_dat(fn,expr_thresh_high)
load(fn)%pSmad2_inproducing
dat = cat(1,positions.nodist);
dat(dat>expr_thresh_high)=0;
expr_mean = mean(nonzeros(dat(isfinite(dat))));
expr_err = std(nonzeros(dat(isfinite(dat))))/power(size(nonzeros(dat(isfinite(dat))),1),0.5);
end
