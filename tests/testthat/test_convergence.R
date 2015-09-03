message("=================")
message("Check convergence")

### Test totally euploid sample ###
file <- list.files(pattern='euploid_binned_data')
model <- findCNVs(file, ID='test', eps=0.1)
w <- model$weights
expect_that(w['disomy'], is_more_than(0.88))
expect_that(w['disomy'], is_less_than(0.90))

### Test aneuploid sample where trisomy is most frequent
file <- list.files(pattern='trisomy_most_frequent_binned_data')
model <- findCNVs(file, ID='test', eps=0.1)
w <- model$weights
expect_that(w['disomy'], is_more_than(0.30))
expect_that(w['disomy'], is_less_than(0.33))
expect_that(w['trisomy'], is_more_than(0.30))
expect_that(w['trisomy'], is_less_than(0.40))

