data(neuroblastoma, package="neuroblastoma")
library(data.table)
nb.dt <- data.table(neuroblastoma[["profiles"]])
some <- nb.dt[as.integer(profile.id)<100]
max.segments <- 5

loss.dt <- some[, {
  this.max <- if(.N < max.segments).N else max.segments
  binseg <- binsegRcpp::binseg_normal(logratio, this.max)
  opt <- jointseg::Fpsn(logratio, this.max)
  data.table(
    segments=1:max.segments,
    binseg=binseg$loss,
    optimal=opt$J)
}, by=.(profile.id, chromosome)]
loss.dt[, diff := binseg - optimal]
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

