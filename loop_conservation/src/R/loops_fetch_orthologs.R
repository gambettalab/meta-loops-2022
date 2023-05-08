options(warn = 1)

# Read config file
config <- yaml::yaml.load_file("config/config.yml")

library(data.table)
library(httr)


#
#  Read the loop data
#

genome <- "D_mel"
other_genome <- "D_vir"

load(paste0("data/loops/loops_", genome, "_to_", other_genome, ".Rdata")) # match


#
#  Function to fetch ortholog annotations from OrthoDB
#

find_orthologs <- function(gene_id, other_tax_id)
{
  message(gene_id)
  result <- NULL

  Sys.sleep(1) # sleep for 1 second
  resp <- GET(paste0("https://data.orthodb.org/current/search?query=", gene_id))
  stopifnot(!http_error(resp))
  stopifnot(http_type(resp) == "application/json")

  for (odb_id in unlist(content(resp)$data)) # for each OrthoDB cluster id
  {
    message("  - ", odb_id)
    resp <- GET(paste0("https://data.orthodb.org/current/orthologs?id=", odb_id))
    stopifnot(!http_error(resp))
    json <- jsonlite::parse_json(resp)

    # limit the orthologs to the ones from other_tax_id species
    sel <- sapply(json$data, function(l) l$organism$id == paste0(other_tax_id, "_0"))
    json$data <- json$data[sel]

    resp <- GET(paste0("https://data.orthodb.org/current/group?id=", odb_id))
    stopifnot(!http_error(resp))
    stopifnot(http_type(resp) == "application/json")
    level_name <- content(resp)$data$level_name

    stopifnot(length(json$data) <= 1)
    if (length(json$data) > 0)
    {
      for (other_gene in json$data[[1]]$genes)
      {
        message("      - ", other_gene$gene_id$id)
        param <- other_gene$gene_id$param

        resp <- GET(paste0("https://data.orthodb.org/current/ogdetails?id=", param))
        stopifnot(!http_error(resp))
        stopifnot(http_type(resp) == "application/json")

        for (xref in content(resp)$data$xrefs)
        {
          if (xref$type == "FlyBase")
          {
            result <- rbind(result, data.table(
              gene_id = gene_id,
              orthodb_cluster_id = odb_id,
              orthodb_level_name = level_name,
              other_orthodb_gene_id = param,
              other_gene_id = xref$id,
              other_gene_name = other_gene$gene_id$id
            ))
          }
        }
      }
    }
  }

  return(result)
}

if (genome == "D_mel") {
  tax_id <- 7227
} else {
  stop("unknown source species")
}

if (other_genome == "D_vir") {
  other_tax_id <- 7244
} else {
  stop("unknown target species")
}


#
#  Save the ortholog table for genes at D. mel. loop anchors
#

dt <- NULL
anchor_gene_ids <- unique(unlist(strsplit(match$lap$anchor_nearest_gene_id, ", ")))
for (gene_id in anchor_gene_ids)
  dt <- rbind(dt, find_orthologs(gene_id, other_tax_id))

write.table(dt, file = paste0("data/loops/genes_loops_", genome, "_to_", other_genome, ".tsv"),
  quote = F, sep = "\t", row.names = F, col.names = T)


#
#  Save the ortholog table for genes at D. vir. loop anchors
#

lap_other <- fread(paste0("data/loops/loops_", other_genome, "_annotated.tsv"), header = T, sep = "\t")

dt <- NULL
anchor_gene_ids <- unique(unlist(strsplit(lap_other$anchor_nearest_gene_id, ", ")))
for (gene_id in anchor_gene_ids)
  dt <- rbind(dt, find_orthologs(gene_id, tax_id))

write.table(dt, file = paste0("data/loops/genes_loops_", other_genome, "_to_", genome, ".tsv"),
  quote = F, sep = "\t", row.names = F, col.names = T)


#
#  Save the Metazoa-level orthologs for genes at D. vir. TSS-proximal and TSS-distal loop anchors
#

stopifnot(lap_other$anchor_TSS_proximal %in% c(0, 1))
proximal_anchor_gene_ids <- unique(unlist(strsplit(lap_other[anchor_TSS_proximal == 1, ]$anchor_nearest_gene_id, ", ")))
distal_anchor_gene_ids <- unique(unlist(strsplit(lap_other[anchor_TSS_proximal == 0, ]$anchor_nearest_gene_id, ", ")))

write.table(dt[orthodb_level_name == "Metazoa" & gene_id %in% proximal_anchor_gene_ids],
  file = paste0("data/loops/genes_proximal_loops_", other_genome, "_to_", genome, ".tsv"),
  quote = F, sep = "\t", row.names = F, col.names = T)

write.table(dt[orthodb_level_name == "Metazoa" & gene_id %in% distal_anchor_gene_ids],
  file = paste0("data/loops/genes_distal_loops_", other_genome, "_to_", genome, ".tsv"),
  quote = F, sep = "\t", row.names = F, col.names = T)
