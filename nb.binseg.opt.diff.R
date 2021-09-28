data(neuroblastoma, package="neuroblastoma")
library(data.table)
nb.dt <- data.table(neuroblastoma[["profiles"]])
some <- nb.dt#[as.integer(profile.id)<100]
max.segments <- 5

loss.dt <- some[, {
  this.max <- if(.N < max.segments).N else max.segments
  binseg <- binsegRcpp::binseg_normal(logratio, this.max)
  opt <- jointseg::Fpsn(logratio, this.max)
  print(profile.id)
  data.table(
    coef(binseg)[, .(binseg_minSize=min(end-start+1)), by=segments],
    optimal_minSize=sapply(1:this.max, function(k){
      end <- opt$t.est[k,1:k]
      start <- c(1, end[-k]+1)
      min(end-start+1)
    }),
    binseg_loss=binseg$loss,
    optimal_loss=opt$J)
}, by=.(profile.id, chromosome)]

loss.dt[, diff := binseg_loss - optimal_loss]
loss.dt[order(-diff)][1:20]
## > loss.dt[order(-diff)][1:20]
##     profile.id chromosome segments    binseg   optimal      diff
##  1:          2          2        3 83.447805 10.497219 72.950586
##  2:         68          2        3 36.461920  5.627889 30.834031
##  3:         79          1        4 35.307299 20.601278 14.706021
##  4:         79          1        3 35.602470 21.841733 13.760737
##  5:         79          1        5 31.570036 18.279781 13.290255
##  6:         22         12        5 27.809012 16.965207 10.843805
##  7:         22         12        3 40.200220 30.471351  9.728870
##  8:          9          1        3 34.627397 24.999328  9.628069
##  9:         15          6        3  9.411855  1.483113  7.928742
## 10:         87          Y        3 21.763455 14.022949  7.740506
## 11:          9          1        4 17.210658  9.564846  7.645812
## 12:         15          6        4  7.393322  1.129333  6.263989
## 13:          9          Y        3 26.249142 20.088606  6.160537
## 14:         53          Y        3 29.470561 24.075606  5.394955
## 15:        101          Y        4 11.131917  6.999434  4.132483
## 16:         64          X        4  5.966889  2.169986  3.796904
## 17:         64          X        3  6.228569  2.520701  3.707868
## 18:         14          Y        4 14.355923 10.743628  3.612295
## 19:         98          X        3  4.349869  0.976401  3.373468
## 20:         14          Y        3 18.273989 15.477230  2.796759

one <- loss.dt[chromosome=="2" & profile.id=="2"]
one[, rev.cum.min := rev(cummin(rev(diff)))]
loss.dt[, .(
  min=min(diff),
  max=max(diff)
), by=segments]
diff.wide <- dcast(
  loss.dt, profile.id + chromosome ~ segments, value.var="diff")
library(ggplot2)

loss.dt[order(-diff)][1:20]
loss.dt[profile.id=="229" & chromosome=="X"]
ggplot()+
  geom_point(aes(
    position, logratio),
    data=nb.dt[profile.id=="229" & chromosome=="X"])
    
loss.dt[profile.id=="565" & chromosome=="X"]
ggplot()+
  geom_point(aes(
    position, logratio),
    data=nb.dt[profile.id=="565" & chromosome=="X"])

ggplot()+
  geom_point(aes(
    `3`, `4`),
    data=diff.wide)+
  scale_x_log10()+
  scale_y_log10()

ggplot()+
  geom_point(aes(
    `5`, `4`),
    data=diff.wide)+
  scale_x_log10()+
  scale_y_log10()

diff.wide[order(-`5`)][1:20]

ggplot()+
  geom_point(aes(
    position, logratio),
    data=nb.dt[chromosome=="1" & profile.id=="590"])



ggplot()+
  geom_point(aes(
    position, logratio),
    data=nb.dt[chromosome=="2" & profile.id=="2"])

max.dt <- loss.dt[, .SD[which.max(diff)], by=.(profile.id, chromosome)]
max.dt[order(-diff)][1:20]

one <- nb.dt[profile.id=="565" & chromosome=="X"]
logratio <- one[["logratio"]]
binseg <- binsegRcpp::binseg_normal(logratio, max.segments)
opt <- jointseg::Fpsn(logratio, max.segments)
cum.end <- cumsum(logratio)
cum.start <- c(0,cum.end)
end <- t(opt[["t.est"]])
start <- rbind(1, end[-max.segments,]+1)
n.data <- end-start+1
seg.sum <- cum.end[end]-cum.start[start]
segs.dt <- rbind(
  data.table(
    algorithm="optimal",
    segments=as.integer(col(end)),
    start=as.integer(start),
    end=as.integer(end),
    mean=as.numeric(seg.sum/n.data)
  )[!is.na(mean)],
  data.table(algorithm="binseg", coef(binseg)))
one[, data.i := .I]
change.dt <- segs.dt[1 < start]
ggplot()+
  geom_point(aes(
    data.i, logratio),
    shape=1,
    data=one)+
  geom_segment(aes(
    start-0.5, mean,
    color=algorithm,
    size=algorithm,
    xend=end+0.5, yend=mean),
    data=segs.dt)+
  geom_vline(aes(
    xintercept=start-0.5,
    size=algorithm,
    color=algorithm),
    data=change.dt)+
  scale_size_manual(values=c(optimal=2, binseg=1))+
  facet_grid(segments ~ ., labeller=label_both)
