noise <- read.table(file='noise-data.csv', header=TRUE, sep=',')

scaled <- scale(noise)

write.table(x=scaled,file='noise-data-scaled.data')

library('rgl')

plot3d(scaled)