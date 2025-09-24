# Load necessary library
library(dplyr) # Optional, for easier data manipulation
library(ggplot2)

power_f <- read.csv('ss_f_0.02_800-2800_3centres.csv')
power_m <- read.csv('ss_m_0.01_1200-3200_3centres.csv')
all_data <- rbind(power_m, power_f)
##power curves by sex
ggplot(all_data, aes(x = ss, y = power, colour = factor(ni.margin, levels = c(0.02, 0.01)))) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "purple") +
  scale_x_continuous(breaks=seq(800, 3200, 400), lim=c(800, 3200)) +
  facet_wrap(~sex, labeller = as_labeller(c(f = "Female", m = "Male"))) +
  labs(x='Sample size', y='Power', color='NI margin') +
  theme_bw()
ggsave('powercurves.jpeg', height = 4, width = 8, dpi=300)
