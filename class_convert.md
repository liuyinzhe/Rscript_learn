
#### matrix
```Rscript
test_matrix <- matrix(c(1,2,3,4,2,1,2,1,3,2,1,3,4,1,3,1), nrow=4)
```

#### matrix > data.frame
```
test_df = as.data.frame(test_matrix)
rownames(test_df) = paste("name_", 1:4, sep="")
colnames(test_df) = paste("name_", 1:4, sep="")
```

#### matrix > dist
```
as.dist(test_matrix)
```

#### data.frame > dist
```
as.dist(test_df)
```

#### dist > matrix
```
as.matrix(test_dist)
```
#### dist > data.frame
```
as.matrix(test_dist)
```
