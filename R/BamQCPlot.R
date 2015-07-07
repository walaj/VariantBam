library(optparse)
#RLIBDIR = '/cga/meyerson/home/marcin/Software/R/x86_64-unknown-linux-gnu-library/3.1/'
#GIT.HOME = '/cga/meyerson/home/marcin/DB/git'

option_list = list(
    make_option(c("-i", "--input"),  type = "character", default = "qcreport.txt",  help = "Input txt file from a snowman preprocess qcreport.txt"),
    make_option(c("-o", "--output"), type = "character", default = "qcreport.pdf",  help = "Output pdf to generate"),
    make_option(c("-r", "--readgroup"), type = "character", default = NULL,  help = "Read group")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$input))
    stop(print_help(parseobj))

print(getwd())
if (!file.exists(opt$input))
  stop(paste("Input file does not exist", opt$input))

print(opt)

require(ggplot2)
require(reshape2)
require(gridExtra)

## read the table
con  <- file(opt$input, open = "r")

rg <- df.mapq <- df.nm <- df.isize <- df.as <- df.xp <- df.len <- df.phred <- df.clip <- data.frame()

## function for parsing histogram strings
parse_hist <- function(ll, id, rg) {
  ss <- strsplit(ll[id],',')[[1]]
  df <- data.frame(matrix(as.numeric(unlist(strsplit(ss, '_'))), nrow=length(ss), byrow=T), stringsAsFactors=FALSE)
  colnames(df) <- c("Min", "Max", "Count")
  df$readgroup = rg
  return (df)
}

line <- readLines(con, n = 1, warn = FALSE)
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {

  ll <- strsplit(line, '\t')[[1]]
  
  thisrg = ll[1]
  rg <- rbind(rg, data.frame(readgroup = thisrg, ReadCount = as.numeric(ll[2]), Supplementary = as.numeric(ll[3]),
                             Unmapped = as.numeric(ll[4]), MateUnmapped = as.numeric(ll[5]),
                             QCFailed = as.numeric(ll[6]), Duplicate = as.numeric(ll[7])))

  df.mapq <-  rbind(df.mapq, parse_hist(ll, 8, thisrg))
  df.nm <-    rbind(df.nm, parse_hist(ll, 9, thisrg))
  df.isize <- rbind(df.isize, parse_hist(ll, 10, thisrg))
  df.clip <-  rbind(df.clip, parse_hist(ll, 11, thisrg))
  #df.phred <- rbind(df.phred, parse_hist(ll, 12, thisrg))
  df.length <- rbind(df.length, parse_hist(ll, 13, thisrg))
  
}

close(con)

if (is.null(opt$readgroup))
  opt$readgroup = unique(rg$readgroup)


cbbPalette <- rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), 50)

#g.mapq  <- ggplot(df.mapq, aes(x=as.numeric(Min), y=pmax(log(Count,10),0), group=readgroup, color=readgroup)) + geom_line() + theme_bw() + scale_x_continuous(limits=c(30,61))    + ylab('Reads') + xlab('Mapping Quality')
g.mapq  <- ggplot(df.mapq[df.mapq$readgroup %in% opt$readgroup, ], aes(x=Min, y=pmax(log(Count, 10),0), group=readgroup, color=readgroup, fill=readgroup)) + geom_bar(stat='identity', position='dodge') + theme_bw() + ylab('Reads') + xlab('Mapping Quality') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^'))) 
g.nm    <- ggplot(df.nm.m    <- melt(df.nm, id='readgroup'),    aes(x=as.numeric(variable), y=pmax(log(value,10),0), group=readgroup, color=readgroup)) + geom_line() + theme_bw() + scale_x_continuous(limits=c(0,50))     + ylab('Reads') + xlab('NM Tag') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^')))
g.isize <- ggplot(df.isize[df.isize$readgroup %in% opt$readgroup, ], aes(x=Min, y=pmax(log(Count,10),0), group=readgroup, color=readgroup, fill=readgroup)) + geom_line() + theme_bw() + ylab('Reads') + xlab('InsertSize') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^')))  + scale_color_manual("Read Group", values=cbbPalette) + theme(legend.position = 'none') + coord_cartesian(xlim=c(200,500))
g.xp    <- ggplot(df.xp.m    <- melt(df.xp, id='readgroup'),    aes(x=as.numeric(variable), y=pmax(log(value,10),0), group=readgroup, color=readgroup)) + geom_line() + theme_bw() + scale_x_continuous(limits=c(0,100))    + ylab('Reads') + xlab('XP Tag') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^'))) 
g.as    <- ggplot(df.as.m    <- melt(df.as, id='readgroup'),    aes(x=as.numeric(variable), y=pmax(log(value,10),0), group=readgroup, color=readgroup)) + geom_line() + theme_bw() + scale_x_continuous(limits=c(0,100))    + ylab('Reads') + xlab('AS Tag') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^'))) 
g.len   <- ggplot(df.len.m   <- melt(df.len, id='readgroup'),   aes(x=as.numeric(variable), y=pmax(log(value,10),0), group=readgroup, color=readgroup)) + geom_line() + theme_bw() + scale_x_continuous(limits=c(20,102))   + ylab('Reads') + xlab('Read Length') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^')))
g.clip <- ggplot(df.clip[df.clip$readgroup %in% opt$readgroup, ], aes(x=Min, y=pmax(log(Count,10),0), group=readgroup, color=readgroup, fill=readgroup)) + geom_line() + theme_bw() + ylab('Reads') + xlab('Clipped Bases') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^')))  + scale_color_manual("Read Group", values=cbbPalette)
g.phred <- ggplot(df.phred.m <- melt(df.phred, id='readgroup'), aes(x=as.numeric(variable), y=pmax(log(value,10),0), group=readgroup, color=readgroup)) + geom_line() + theme_bw() + scale_x_continuous(limits=c(0,43))     + ylab('Reads') + xlab('Mean read Phred quality') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^'))) 


#pdf(opt$output, width=15, height=8)
ppdf(print(grid.arrange(g.isize, g.clip, ncol=2)), width=4, height=6)
#dev.off()

