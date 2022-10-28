zy_calc_tpm <- function(dt=NA, le=NA, ID="name",le_name="length"){
  rownames(le) = le[,ID]
  le = le[rownames(dt),]
  dt = dt/le[,le_name]
  dt = t(t(dt)/colSums(dt))
  dt
}
