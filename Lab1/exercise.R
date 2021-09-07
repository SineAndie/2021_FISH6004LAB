mat1 = matrix(rep(NA, 30), nrow = 3, ncol = 10)
mat2 = matrix(rep(NA, 30), nrow = 3, ncol = 10)

A = runif(30, min = 1, max = 5)
mat3 = matrix(A, nrow = 3, ncol = 10)
mat4 = matrix(runif(30, min = 1, max = 5), nrow = 3, ncol = 10)

for ( i in 1:3)
{
  for (j in 1:10){
    mat3[i,j] = runif (1, min = 1, max = 5)
  }
} 


mat3%*%mat4

mat3%*%t(mat4)

B = as.data.frame(mat3%*%t(mat4))
write.table(B, file = "result.csv")


Data <- read.csv('cod.csv', header = TRUE)
summary(Data)

C = subset(Data, YEAR >= "1990" & YEAR <= "1995")

E = cbind(rep(1,dim(C)[1]),C)
names(E)[1] = "id"

D = tapply(E$id, INDEX = E$YEAR, FUN = sum, na.rm = TRUE)
barplot(D)

F = tapply(C$WTCPUE, INDEX = E$YEAR, FUN = mean, na.rm = TRUE)



plot(F, type='l')

result <- lm(WTCPUE~BOT_TEMP, data = C)
result


