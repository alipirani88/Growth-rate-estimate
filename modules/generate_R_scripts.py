__author__ = 'alipirani'
import os

def generate_perc_coverage_graph(csv_matrix_file, ptr_value):
    out_path = os.path.dirname(csv_matrix_file)
    header = os.path.basename(csv_matrix_file).replace('.csv', '')
    with open(out_path + "/perc_coverage_graph.R", 'w') as out:
        print_string = "library(reshape)\nlibrary(ggplot2)\n"
        print_string = print_string + "dat=read.csv(\"%s\")\n" % (os.path.basename(csv_matrix_file))
        print_string = print_string + "mdf=melt(dat,id.vars=\"bin\")\n"
        print_string = print_string + "pdf(\"perc_coverage_graph.pdf\", width = 15, height = 10)\n"
        print_string = print_string + "p1 <- ggplot(mdf,aes(x=bin,y=value,colour=variable,group=variable)) + geom_line(size=1.5) + ggtitle(\"%s\")  + labs(x=\"Genomic Location\",y=\"Coverage in Percentage (perc reads/total)\") + theme(legend.title = element_blank(), legend.position=\"none\") + annotate(\"text\", x=max(mdf$bin)/2, y=max(mdf$value), label= \"PTR = %s\", size=5) + theme(text = element_text(size=10), panel.background = element_rect(fill = 'white', colour = 'white'), plot.title = element_text(size=10, face=\"bold\", margin = margin(10, 0, 10, 0)), axis.text.x = element_text(colour = \"black\", face= \"bold.italic\", vjust = 3)) + scale_color_manual(values=c(\"#3288bd\"))\n" % (header, ptr_value)
        print_string = print_string + "p1\ndev.off()"
        out.write(print_string)



