# Define the parameters
f0 <- 3
x <- seq(0, 100, length.out = 1000)  # Time vector
y <- sin(2 * pi * f0 * x)  # Signal

# Perform FFT in R
r_fft <- fft(y)
N <- length(y)  # Number of points

# Calculate corresponding frequencies for the R FFT result
r_freq <- seq(0, N/2 - 1) * (1 / (x[2] - x[1])) / N

# Periodogram (Magnitude squared of the FFT)
r_periodogram <- Mod(r_fft[1:(N/2)])^2

# Install the reticulate package to call Python code from R
if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}

library(reticulate)

# Call Python's FFT function using reticulate
np <- import("numpy")
py_fft <- np$fft$fft(y)
py_freq <- np$fft$fftfreq(N, d = x[2] - x[1])

# Periodogram for Python (Magnitude squared of the FFT)
py_periodogram <- Mod(py_fft[1:(N/2)])^2

# Save the plot as a PNG
png("fft_periodogram_comparison.png", width = 800, height = 600)

# Plot the periodograms
plot(r_freq, r_periodogram, type = "l", col = "blue", main = "Periodogram Comparison (R vs Python)", 
     xlab = "Frequency (Hz)", ylab = "Power", lwd = 2)
lines(py_freq[1:(N/2)], py_periodogram, col = "red", lwd = 2)
legend("topright", legend = c("R FFT", "Python FFT"), col = c("blue", "red"), lty = 1, lwd = 2)

# Close the device and save the plot
dev.off()

# Inform the user where the plot is saved
cat("Periodogram plot saved as fft_periodogram_comparison.png\n")
