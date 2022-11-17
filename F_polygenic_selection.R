setwd('/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/selection_manuscript/')

data=read.delim('For_figures/30Aug22_polygenic_selection.txt')
data=subset(data,Biomarker!='Oestradiol')

library(ggplot2)
ggplot(data, aes(x=reorder(Biomarker,Mean.FCS..observed.), y=Mean.FCS..background.)) + 
  geom_line() +
  geom_point()+coord_flip()+theme_bw(13)+
  geom_errorbar(aes(ymin=Mean.FCS..background.-SD.FCS..background.*2, ymax=Mean.FCS..background.+SD.FCS..background.*2), width=.2,
                 position=position_dodge(0.05))+geom_point(aes(x=Biomarker, y=Mean.FCS..observed.,col='red'))