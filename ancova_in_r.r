df <- read.table(“W1_gatk_rates_for_ancova.txt”, header = TRUE)
attach(df)
nuclist <- unique(df$nucs)
quallist <- unique(df$qual)
for (i in nuclist) {
	nucdata<-df[df$nucs %in% i,]
	for (j in quallist) {
		nucqualdata<-nucdata[nucdata$qual %in% j,]
		polyfit <- function(i) x <- AIC(lm(nucqualdata$rate~poly(nucqualdata$pos,i, raw = TRUE)))
		a = as.integer(optimize(polyfit,interval = c(1,6))$minimum)
		fit <- lm(nucqualdata$rate~poly(nucqualdata$pos,a, raw = TRUE))
		b = coefficients(fit)
		x<-c(i,j,b)
		write(x, file = "W1_gatk_polynomial_regressions.out", append = TRUE, sep = "\t")
	}
}
