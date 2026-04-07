
library(dplyr)

es %>%
  summarise(
    n_negative = sum(yi < -0.2),
    n_neutral  = sum(yi >= -0.2 & yi <= 0.2),
    n_positive = sum(yi > 0.2)
  )


es %>%
  filter(!is.na(yi)) %>%
  summarise(
    pct_negative = 100 * sum(yi < -0.2) / n(),
    pct_neutral  = 100 * sum(yi >= -0.2 & yi <= 0.2) / n(),
    pct_positive = 100 * sum(yi > 0.2) / n()
  )


length(unique(class.es$order_ncbi))
table(class.es$class_ncbi)

length(unique(es$Study_Country))
table(es$Study_Country)


es_summary <- data.frame(
  yi = es$yi,
  vi = es$vi,
  lower = es$yi - 1.96 * sqrt(es$vi),
  upper = es$yi + 1.96 * sqrt(es$vi)
)

# Add logical column indicating whether CI overlaps zero
es_summary$overlaps_zero <- es_summary$lower < 0 & es_summary$upper > 0

# Summary: counts
summary_table <- table(es_summary$overlaps_zero)
names(summary_table) <- c("Does NOT Overlap Zero", "Overlaps Zero")
print(summary_table)

# Optional: proportions
prop_table <- prop.table(summary_table)
print(round(prop_table, 3))
