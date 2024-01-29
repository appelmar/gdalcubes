#f <- function() {
x <- read_chunk_as_array(with.dimnames = TRUE)
d = dim(x)
bnames = dimnames(x)[[1]]
nvar = d[1]
nobs = prod(d[2:4])
dim(x) = c(nvar, nobs)
df = data.frame(t(x))
colnames(df) <- bnames
idx = which(complete.cases(df))
df = df[idx,]
args$newdata = df
args$object = model
p = do.call("predict", args)
out = array(NA, dim=c(1,d[2:4]))
out[idx] = p
write_chunk_from_array(out)
#}