function makeplot(hObject,event,analysis,duration,x,j,hplot,hplot2)
trial = get(hObject,'Value');
set(hplot,'ydata',polyval(analysis.(duration{j}).SW.pop_fit_params(trial,:),x));
set(hplot2,'ydata',analysis.(duration{j}).SW.pop_est_m(trial,:));
title(num2str(trial))
drawnow;
end