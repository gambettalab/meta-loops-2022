library(CNEr)

assemblyDir <- "data/lastz/2bit"
axtDir <- "data/lastz/axt"
suppressWarnings(dir.create(assemblyDir))
suppressWarnings(dir.create(axtDir))

## convert FASTA to 2bit format
system("faToTwoBit data/lastz/fasta/D_mel.fa data/lastz/2bit/D_mel.2bit")
system("faToTwoBit data/lastz/fasta/D_vir.fa data/lastz/2bit/D_vir.2bit")

make_net_axt <- function(assemblyTarget, assemblyQuery)
{
  ## lastz aligner
  lavs <- lastz(assemblyTarget, assemblyQuery,
    outputDir=axtDir,
    distance="medium", mc.cores=64)

  ## lav files to psl files conversion
  psls <- lavToPsl(lavs, removeLav=FALSE, binary="lavToPsl")


  ## Join close alignments
  chains <- axtChain(psls, assemblyTarget=assemblyTarget,
    assemblyQuery=assemblyQuery, distance="medium",
    removePsl=FALSE, binary="axtChain")

  ## Sort and combine
  allChain <- chainMergeSort(chains, assemblyTarget, assemblyQuery,
    allChain=file.path(axtDir,
      paste0(sub("\\.2bit$", "", basename(assemblyTarget),
          ignore.case=TRUE), ".", 
        sub("\\.2bit$", "", basename(assemblyQuery), 
          ignore.case=TRUE), ".all.chain")),
      removeChains=FALSE, binary="chainMergeSort")


  ## Filtering out chains
  allPreChain <- chainPreNet(allChain, assemblyTarget, assemblyQuery,
    allPreChain=file.path(axtDir,
      paste0(sub("\\.2bit$", "", basename(assemblyTarget),ignore.case = TRUE), ".", 
        sub("\\.2bit$", "", basename(assemblyQuery), ignore.case = TRUE), ".all.pre.chain")),
    removeAllChain=FALSE, binary="chainPreNet")

  ## Keep the best chain and add synteny information
  netSyntenicFile <- chainNetSyntenic(allPreChain, assemblyTarget, assemblyQuery,
    netSyntenicFile=file.path(axtDir,
      paste0(sub("\\.2bit$", "", basename(assemblyTarget), ignore.case = TRUE), ".",
        sub("\\.2bit$", "", basename(assemblyQuery), ignore.case = TRUE), ".noClass.net")),
    binaryChainNet="chainNet", binaryNetSyntenic="netSyntenic")


  ## Create .net.axt file from the previous net and chain files
  netToAxt(netSyntenicFile, allPreChain, assemblyTarget, assemblyQuery,
    axtFile=file.path(axtDir,
      paste0(sub("\\.2bit$", "", basename(assemblyTarget), ignore.case = TRUE), ".",
        sub("\\.2bit$", "", basename(assemblyQuery), ignore.case = TRUE),
            ".net.axt")),
    removeFiles=FALSE,
    binaryNetToAxt="netToAxt", binaryAxtSort="axtSort")
}


make_net_axt("data/lastz/2bit/D_mel.2bit", "data/lastz/2bit/D_vir.2bit")
make_net_axt("data/lastz/2bit/D_vir.2bit", "data/lastz/2bit/D_mel.2bit")
