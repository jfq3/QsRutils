# Prevent RCMD CHECK from reporting no visible binding for 
# variables created dynamically

utils::globalVariables("Max")
utils::globalVariables(c("seq_length", "taxon"))
