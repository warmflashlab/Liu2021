function saveAllOpenFigures(saveInPath)
%% saving all open figures.
h = get(0,'children');

mkdir([saveInPath])
extn = {'.png', '.fig'}; %save as both pdf(for ppt) and fig(for later changes).

%for ii = 1
for ii=1:length(h)
    for jj = 1:2
        %for i = 27
        %stuff where i is the numbering of the figure *and* the handle to use,
        if jj == 1
          %saveas(h(ii), [saveInPath filesep 'figure_' num2str(ii) extn{jj}], 'epsc');
          saveas(h(ii), [saveInPath filesep 'figure_' num2str(ii) extn{jj}])
        else
        saveas(h(ii), [saveInPath filesep 'figure_' num2str(ii) extn{jj}]);
        end
        
    end
end
%%

