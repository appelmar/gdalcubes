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
if ("model_fit" %in% class(model)) {
  args$new_data = df # tidymodels fix
} else {
  args$newdata = df
}
args$object = model
p = do.call("predict", args)
out = array(NA, dim=c(length(output_names),nobs))
if (is.list(p)) {
  for (i in 1:length(output_names)) {
    z = p[[output_names[i]]]
    if (!is.numeric(z) && !is.factor(z)) {
      stop("Unexpected output from predict().")
    }
    out[i,idx] = as.numeric(z)
  }
} else if (is.matrix(p)) {
  for (i in 1:length(output_names)) {
    out[i,idx] = p[,output_names[i]]
  }
}  else if (is.numeric(p) || is.factor(p)) {
  out[1, idx] = p
} else {
  stop("Unexpected output from predict().")
}
dim(out) <- c(length(output_names),d[2:4])
write_chunk_from_array(out)
#}