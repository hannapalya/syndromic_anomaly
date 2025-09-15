# Simple test script to check if R is working
print("Starting simple test...")

# Test basic R functionality
x <- 1:10
print(paste("Basic test passed. x =", length(x)))

# Test CSV reading
print("Testing CSV reading...")
bankholidays <- read.csv("Bankholidays.csv", sep=',')
print(paste("CSV read successfully. Rows:", nrow(bankholidays)))

# Test simple function
test_func <- function(n) {
  return(n * 2)
}
print(paste("Function test:", test_func(5)))

print("All basic tests passed!")
