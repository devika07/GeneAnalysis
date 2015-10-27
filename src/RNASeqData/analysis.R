###################################################################################################################################
# PROJECT FOR GENE EXPRESSION ANALYSIS
# TOPIC 2 - from RNA Seq dataset
###################################################################################################################################


inputFile <- "/drosophiladata/ERR029113.fastq"
con  <- file(inputFile, open = "r")
fileConn <-file("4000head.fastq")
oneLine <- readLines(con, n = 4000, warn = FALSE)
writeLines(oneLine, fileConn)
close(fileConn)
close(con)
