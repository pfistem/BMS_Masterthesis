# Drop the loser for MAMS | Binary Endpoint
# Bristol Myers Squibb
# Contact: manuel.pfister@bms.com

# -------------------------------------------------------------------------

finddesign_binary <- function(J, K, ns, delta1, delta0, requiredfwer, requiredpower, groupsizes, treatmentsigmas) {
  # Code specific to binary outcomes
  # ...
  return(list(n = n, c = c, totalSS = n * sum((ns + 1) * groupsizes / groupsizes[1])))
}

# -------------------------------------------------------------------------