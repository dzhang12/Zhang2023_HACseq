
f.data <- read.table("all.f.coverage.txt")
f.data.filtered <- f.data[-read.table("all.f.coverage.quant.txt")[[2]],]
write.table(f.data.filtered, "f.coverage.filtered.txt", row.names = F, col.names = F, quote = F)

r.data <- read.table("all.r.coverage.txt")
r.data.filtered <- r.data[-read.table("all.r.coverage.quant.txt")[[2]],]
write.table(r.data.filtered, "r.coverage.filtered.txt", row.names = F, col.names = F, quote = F)

library(tidyverse)

f.data <- read.table("f.coverage.filtered.txt", header = T)
f.pct <- data.frame(WT1 = f.data$WT1.start/f.data$WT1,
                    WT2 = f.data$WT2.start/f.data$WT2,
                    WT3 = f.data$WT3.start/f.data$WT3,
                    WT4 = f.data$WT4.start/f.data$WT4,
                    KO3 = f.data$KO3.start/f.data$KO3,
                    KO4 = f.data$KO4.start/f.data$KO4,
                    KO5 = f.data$KO5.start/f.data$KO5)
f.pct <- t(f.pct)
f.pct <- as.data.frame(f.pct)

f.stat <- vector("numeric", ncol(f.pct))
for (i in 1:ncol(f.pct)){
  x <- f.pct[[i]]
  p <- try(t.test(x[1:4], x[5:7])$p.value, silent = T)
  if (is.numeric(p)){
    f.stat[i] <- p
  } else {
    f.stat[i] <- NA
  }
  if (i %% 10000 == 0){
    print(i)
  }
}
f.results <- data.frame(Chr = f.data$Chr,
                        Pos = f.data$End,
                        Pval = f.stat)
write.table(f.results, "f.results.txt", row.names = F, quote = F)




r.data <- read.table("r.coverage.filtered.txt", header = T)
r.pct <- data.frame(WT1 = r.data$WT1.start/r.data$WT1,
                    WT2 = r.data$WT2.start/r.data$WT2,
                    WT3 = r.data$WT3.start/r.data$WT3,
                    WT4 = r.data$WT4.start/r.data$WT4,
                    KO3 = r.data$KO3.start/r.data$KO3,
                    KO4 = r.data$KO4.start/r.data$KO4,
                    KO5 = r.data$KO5.start/r.data$KO5)
r.pct <- t(r.pct)
r.pct <- as.data.frame(r.pct)

r.stat <- vector("numeric", ncol(r.pct))
for (i in 1:ncol(r.pct)){
  x <- r.pct[[i]]
  p <- try(t.test(x[1:4], x[5:7])$p.value, silent = T)
  if (is.numeric(p)){
    r.stat[i] <- p
  } else {
    r.stat[i] <- NA
  }
  if (i %% 10000 == 0){
    print(i)
  }
}
r.results <- data.frame(Chr = r.data$Chr,
                        Pos = r.data$End,
                        Pval = r.stat)
write.table(r.results, "r.results.txt", row.names = F, quote = F)



f.data$strand <- "+"
r.data$strand <- "-"
f.data$pval <- f.stat
r.data$pval <- r.stat
all.data <- rbind(f.data, r.data)
all.data <- all.data[,-2]
colnames(all.data)[2] <- "Pos"
head(arrange(all.data, pval), 20)
write.table(all.data, "all.results.txt", row.names = F, quote = F)
