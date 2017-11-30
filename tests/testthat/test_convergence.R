message("=================")
message("Check convergence")

### Test totally euploid sample ###
file <- list.files(pattern='euploid_')
model <- findCNVs(file, ID='test', eps=0.1, most.frequent.state='2-somy', states=c("zero-inflation",paste0(1:10,'-somy')), num.trials=1, method = 'HMM')
w <- model$weights
expect_that(w['2-somy'], is_more_than(0.85))
expect_that(w['2-somy'], is_less_than(0.91))

# With 0-somy
file <- list.files(pattern='euploid_')
model <- findCNVs(file, ID='test', eps=0.1, most.frequent.state='2-somy', states=c("zero-inflation",paste0(0:10,'-somy')), num.trials=1, method = 'HMM')
w <- model$weights
expect_that(w['2-somy'], is_more_than(0.85))
expect_that(w['2-somy'], is_less_than(0.91))

### Test aneuploid sample where 3-somy is most frequent
file <- list.files(pattern='trisomy_')
model <- findCNVs(file, ID='test', eps=0.1, most.frequent.state='2-somy', states=c("zero-inflation",paste0(1:10,'-somy')), num.trials=1, method = 'HMM')
w <- model$weights
expect_that(w['2-somy'], is_more_than(0.30))
expect_that(w['2-somy'], is_less_than(0.33))
expect_that(w['3-somy'], is_more_than(0.30))
expect_that(w['3-somy'], is_less_than(0.40))

# With 0-somy
file <- list.files(pattern='trisomy_')
model <- findCNVs(file, ID='test', eps=0.1, most.frequent.state='2-somy', states=c("zero-inflation",paste0(0:10,'-somy')), num.trials=1, method = 'HMM')
w <- model$weights
expect_that(w['2-somy'], is_more_than(0.30))
expect_that(w['2-somy'], is_less_than(0.33))
expect_that(w['3-somy'], is_more_than(0.30))
expect_that(w['3-somy'], is_less_than(0.40))
