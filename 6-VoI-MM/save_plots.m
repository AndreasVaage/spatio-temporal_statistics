function save_plots(name,dir,figure)
%save plots in PDF.

set(gcf, 'PaperSize', [15 15]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure,'-dpdf',strcat(dir,name));
system(['pdfcrop ',dir,name,'.pdf ',dir,name,'.pdf']);

end

