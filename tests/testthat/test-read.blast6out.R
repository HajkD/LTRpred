context("Test: read.blast6out() ")

test_that("read.blast6out() correctly imports blast output file.",{
    # read example *.blast6out file
    expect_silent(test.blast6out <- read.blast6out(system.file("test.blast6out", package = "LTRpred")))
})

test_that("read.blast6out() correctly imports correct columns.",{
    # read example *.blast6out file
    test.blast6out <- read.blast6out(system.file("test.blast6out", package = "LTRpred"))
    expect_identical(colnames(test.blast6out), c("query", "subject", 
                                                 "perc_ident", "align_len", 
                                                 "n_mismatch", "n_gap_open", 
                                                 "start_q", "end_q", 
                                                 "start_s", "end_s", 
                                                 "evalue", "bit_score"))

})

test_that("read.blast6out() is a tibble.", {
    
    test.blast6out <- read.blast6out(system.file("test.blast6out", package = "LTRpred"))
    
    expect_true(tibble::is.tibble(test.blast6out))
})
          
