options(warn = 1)

# Read config file
config <- yaml::yaml.load_file("config/config.yml")

library(data.table)
library(GenomicRanges)

#
#  Read the loop data
#

genome <- "D_mel"
other_genome <- "D_vir"

# read the "reverse" loop matching (from other_genome to genome)
load(paste0("data/loops/loops_", other_genome, "_to_", genome, ".Rdata")) # match
other_match <- match
# now mask the loaded variables with the "actual" loop matching (from genome to other_genome)
load(paste0("data/loops/loops_", genome, "_to_", other_genome, ".Rdata")) # match

#
#  Check for fidelity of the loop interactions
#

# extract the loops conserved in D. vir.
conserved <- match$candidates[!is.na(other_loop_id), ]
stopifnot(with(conserved, other_loop_id == A1_other_loop_id))
stopifnot(with(conserved, other_loop_id == A2_other_loop_id))

# extract the loops conserved in D. mel.
other_conserved <- other_match$candidates[!is.na(other_loop_id), ]
stopifnot(with(other_conserved, other_loop_id == A1_other_loop_id))
stopifnot(with(other_conserved, other_loop_id == A2_other_loop_id))

stopifnot(setequal(conserved$other_loop_id, other_conserved$loop_id))
stopifnot(setequal(other_conserved$loop_id, conserved$other_loop_id))

# here we look for pairs of conserved anchors that form a loop in D. mel.,
# the loop is not conserved in D. vir., but both anchors take part in some non-conserved D. vir. loops

nonlooping <- match$candidates[!loop_id %in% conserved$loop_id, ]
message("Possibly unfaithful D. mel. loops:")
print(nonlooping[!is.na(A1_other_loop_id) & !is.na(A2_other_loop_id)
  & !(A1_other_loop_id %in% other_conserved$loop_id) & !(A2_other_loop_id %in% other_conserved$loop_id)])
message()

# now we look for pairs of conserved anchors that form a loop in D. vir.,
# the loop is not conserved in D. mel., but both anchors take part in some non-conserved D. mel. loops

other_nonlooping <- other_match$candidates[!loop_id %in% other_conserved$loop_id, ]
message("Possibly unfaithful D. vir. loops:")
print(other_nonlooping[!is.na(A1_other_loop_id) & !is.na(A2_other_loop_id)
  & !(A1_other_loop_id %in% conserved$loop_id) & !(A2_other_loop_id %in% conserved$loop_id)])
message()
