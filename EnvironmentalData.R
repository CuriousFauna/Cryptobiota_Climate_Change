
#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
#  Environmental data from mesocosm experiment - Figure 1
#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

library(gtable)
library(grid)
library(ggplot2)

table1<-read.csv("Daily_mean_temp_chem.csv", header=T)
colorOrder<-c("Heated","AcidifiedHeated","Control","Acidified")
table1$Treatment<-factor(table1$Treatment, levels = colorOrder)

table2<-read.csv("Diel_mean_temp_chem.csv", header=T)
table2$Treatment<-factor(table2$Treatment, levels = colorOrder)

p1<-ggplot(table1, aes(x=Day.past.1.1.16, y=Temp..in.situ., color=Treatment)) + 
  geom_errorbar(aes(ymin=Temp..in.situ.-Temp..in.situ..SD, ymax=Temp..in.situ.+Temp..in.situ..SD), color="black", width=10) +
  geom_line(size=0.5)+
  theme_classic()+theme(legend.position="none")+
  scale_color_manual(values = c(Acidified ="yellow3",
                                AcidifiedHeated ="#999999",
                                Heated ="#E69F00",
                                Control = "#56B4E9"))+
  theme(plot.tag=element_text(size=12))+theme(plot.tag=element_text(face="bold"))+
  geom_hline(aes(yintercept=27.98), linetype="dashed", size=0.5)+
  labs(x="", y="Temperature (°C)", tag="A")+
  scale_x_continuous(breaks=c(182,547,912), limits=c(0,1035), labels=c("", "", ""), expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), limits=c(21,31))+
  theme(plot.margin=margin(10,20,10,10))+
  theme(axis.line=element_line(size=0.001))+
  theme(axis.ticks=element_line(size=0.001))+
  theme(axis.line.x = element_line(color="black", size = .25),
        axis.line.y = element_line(color="black", size = .25),
        axis.title.y = element_text(size = 10),
        axis.text=element_text(size=8))


p2<-ggplot(table2, aes(x=Time.past.start, y=Temp.mean, color=Treatment)) + 
  geom_errorbar(aes(ymin=Temp.mean-Temp.SD, ymax=Temp.mean+Temp.SD), color="black", width=0.5) +
  geom_line(size=0.5)+
  theme_classic()+theme(legend.position="none")+
  scale_color_manual(values = c(Acidified ="yellow3",
                                AcidifiedHeated ="#999999",
                                Heated ="#E69F00",
                                Control = "#56B4E9"))+
  theme(plot.tag=element_text(size=12))+theme(plot.tag=element_text(face="bold"))+
  theme(axis.title.y=element_text(size=10))+
  theme(axis.text.y=element_blank())+
  theme(axis.text.x=element_text(size=8))+
  theme(axis.title.x=element_text(size=10))+
  labs(x="", y="", tag="B")+
  scale_x_continuous(breaks=c(-4,8,20), limits=c(-4,24), labels=c("","",""), expand=c(0,0))+
  scale_y_continuous(limits=c(21,31), expand=c(0,0))+
  theme(plot.margin=margin(10,20,10,10))+
  theme(axis.line=element_line(size=0.001))+
  theme(axis.ticks=element_line(size=0.001))+
  theme(axis.line.x = element_line(color="black", size = .25),
        axis.line.y = element_line(color="black", size = .25))

p3<-ggplot(table1, aes(x=Day.past.1.1.16, y=pH.out, color=Treatment)) + 
  geom_errorbar(aes(ymin=pH.out-pH.out.SD, ymax=pH.out+pH.out.SD), color="black", width=10) +
  geom_line(size=0.5)+
  theme_classic()+theme(legend.position="none")+
  scale_color_manual(values = c(Acidified ="yellow3",
                                AcidifiedHeated ="#999999",
                                Heated ="#E69F00",
                                Control = "#56B4E9"))+
  theme(plot.tag=element_text(size=12))+theme(plot.tag=element_text(face="bold"))+
  labs(x="Date", y=expression("pH"), tag="C")+
  scale_x_continuous(breaks=c(182,547,912), limits=c(182,912), labels=c("July\n2016", "July\n2017", "July\n2018"), expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), limits=c(7.5,8.3))+
  theme(plot.margin=margin(10,20,10,10))+
  theme(axis.line=element_line(size=0.001))+
  theme(axis.ticks=element_line(size=0.001))+
  theme(axis.line.x = element_line(color="black", size = .25),
        axis.line.y = element_line(color="black", size = .25),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text=element_text(size=8))


p4<-ggplot(table2, aes(x=Time.past.start, y=pH.mean, color=Treatment)) + 
  geom_errorbar(aes(ymin=pH.mean-pH.SD, ymax=pH.mean+pH.SD), color="black", width=0.5) +
  geom_line(size=0.5)+
  theme_classic()+theme(legend.position="none")+
  scale_color_manual(values = c(Acidified ="yellow3",
                                AcidifiedHeated ="#999999",
                                Heated ="#E69F00",
                                Control = "#56B4E9"))+
  theme(plot.tag=element_text(size=12))+theme(plot.tag=element_text(face="bold"))+
  theme(axis.title.y=element_text(size=10))+
  theme(axis.text.y=element_blank())+
  theme(axis.text.x=element_text(size=8))+
  theme(axis.title.x=element_text(size=10))+
  labs(x="Time", y="", tag="D")+
  scale_x_continuous(breaks=c(-4,8,20), limits=c(-4,24), labels=c("0800","2000","0800"), expand=c(0,0))+
  scale_y_continuous(limits=c(7.5,8.3), expand=c(0,0))+
  theme(plot.margin=margin(10,20,10,10))+
  theme(axis.line=element_line(size=0.001))+
  theme(axis.ticks=element_line(size=0.001))+
  theme(axis.line.x = element_line(color="black", size = .25),
        axis.line.y = element_line(color="black", size = .25),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text=element_text(size=8))

g1<-ggplotGrob(p1)
g3<-ggplotGrob(p3)
g2<-ggplotGrob(p2)
g4<-ggplotGrob(p4)
g<-rbind(g1, g3, size="first")
g$widths<-unit.pmax(g1$widths, g3$widths)
gg<-rbind(g2, g4, size="first")
gg$widths<-unit.pmax(g2$widths, g4$widths)
ggg<-cbind(g, gg, size="first")
grid.newpage()
grid.draw(ggg)
