library(tidyverse)

hurricane <- readr::read_csv("HourlyMaxWinds18512006.csv")
hurricane %>% select(Yr) %>% unique() %>% range()

# histogram for data
hurricane %>% select(Yr, Wmax) %>% subset(Yr == 2006) %>% 
  ggplot(aes(x = Wmax)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  labs(x = "Maximum Wind Speed (Wmax)", y = "Density") +
  theme_minimal() + theme_bw()
