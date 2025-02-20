setwd("C:/Users/andre/Box/New folder")

# Function to calculate FST(t) based on the provided equation
fst_drift <- function(t, FST0, Ne) {
  FST_t <- 1 - (1 - FST0) * exp(-t / (2 * Ne))
  return(FST_t)
}

# Parameters
Nc <- 13000     # Census population size

# Subgroup effective population sizes
Ne_central <- 0.15 * Nc  # Central subgroup: 15% of total population
Ne_west <- 0.50 * Nc     # West subgroup: 50% of total population
Ne_east <- 0.35 * Nc     # East subgroup: 35% of total population

# Pairwise harmonic means for effective population size (small and large Ne scenarios)
Ne_central_west_small <- 2 * (0.1 * Ne_central) * (0.1 * Ne_west) / ((0.1 * Ne_central) + (0.1 * Ne_west))
Ne_central_west_large <- 2 * (0.3 * Ne_central) * (0.3 * Ne_west) / ((0.3 * Ne_central) + (0.3 * Ne_west))

Ne_central_east_small <- 2 * (0.1 * Ne_central) * (0.1 * Ne_east) / ((0.1 * Ne_central) + (0.1 * Ne_east))
Ne_central_east_large <- 2 * (0.3 * Ne_central) * (0.3 * Ne_east) / ((0.3 * Ne_central) + (0.3 * Ne_east))

Ne_east_west_small <- 2 * (0.1 * Ne_east) * (0.1 * Ne_west) / ((0.1 * Ne_east) + (0.1 * Ne_west))
Ne_east_west_large <- 2 * (0.3 * Ne_east) * (0.3 * Ne_west) / ((0.3 * Ne_east) + (0.3 * Ne_west))

t <- c(0, 10, 30, 50)  # Time in generations: 0, 10, 30, and 50

# Observed pairwise genetic differentiation at time 0
FST0_central_west <- 0.0008
FST0_central_east <- 0.0004
FST0_east_west <- 0.0000

# Calculate FST(t) for small and large Ne scenarios
FST_central_west_small <- fst_drift(t, FST0_central_west, Ne_central_west_small)
FST_central_west_large <- fst_drift(t, FST0_central_west, Ne_central_west_large)
FST_central_west_avg <- (FST_central_west_small + FST_central_west_large) / 2

FST_central_east_small <- fst_drift(t, FST0_central_east, Ne_central_east_small)
FST_central_east_large <- fst_drift(t, FST0_central_east, Ne_central_east_large)
FST_central_east_avg <- (FST_central_east_small + FST_central_east_large) / 2

FST_east_west_small <- fst_drift(t, FST0_east_west, Ne_east_west_small)
FST_east_west_large <- fst_drift(t, FST0_east_west, Ne_east_west_large)
FST_east_west_avg <- (FST_east_west_small + FST_east_west_large) / 2

# Combine results into data frames for plotting
results_central_west <- data.frame(
  Generation = t,
  FST_Min = FST_central_west_small,
  FST_Max = FST_central_west_large,
  FST_Avg = FST_central_west_avg
)

results_central_east <- data.frame(
  Generation = t,
  FST_Min = FST_central_east_small,
  FST_Max = FST_central_east_large,
  FST_Avg = FST_central_east_avg
)

results_east_west <- data.frame(
  Generation = t,
  FST_Min = FST_east_west_small,
  FST_Max = FST_east_west_large,
  FST_Avg = FST_east_west_avg
)

# Plot FST(t) over time for each pairwise comparison
library(ggplot2)

# Central-West plot
ggplot(results_central_west, aes(x = Generation)) +
  geom_ribbon(aes(ymin = FST_Min, ymax = FST_Max), fill = "#3f8eaa", alpha = 0.3) +
  geom_line(aes(y = FST_Avg), color = "#3f8eaa", size = 1) +
  labs(title = "Change in Genetic Differentiation: Central-West",
       x = "Generation (t)",
       y = "FST(t)") +
  ylim(0, 0.1) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),  # Remove background grid
    panel.background = element_rect(fill = "white")  # Set background to white
  )

# Central-East plot
ggplot(results_central_east, aes(x = Generation)) +
  geom_ribbon(aes(ymin = FST_Min, ymax = FST_Max), fill = "#99cd91", alpha = 0.3) +
  geom_line(aes(y = FST_Avg), color = "#99cd91", size = 1) +
  labs(title = "Change in Genetic Differentiation: Central-East",
       x = "Generation (t)",
       y = "FST(t)") +
  ylim(0, 0.1) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),  # Remove background grid
    panel.background = element_rect(fill = "white")  # Set background to white
  )

# East-West plot
ggplot(results_east_west, aes(x = Generation)) +
  geom_ribbon(aes(ymin = FST_Min, ymax = FST_Max), fill = "#f06c45", alpha = 0.3) +
  geom_line(aes(y = FST_Avg), color = "#f06c45", size = 1) +
  labs(title = "Change in Genetic Differentiation: East-West",
       x = "Generation (t)",
       y = "FST(t)") +
  ylim(0, 0.1) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),  # Remove background grid
    panel.background = element_rect(fill = "white")  # Set background to white
  )
###################################################################################
##Togiak vs mch 
# Function to calculate FST(t) based on the provided equation
fst_drift <- function(t, FST0, Ne) {
  FST_t <- 1 - (1 - FST0) * exp(-t / (2 * Ne))
  return(FST_t)
}

# Parameters
Nc_mch <- 13000     # MCH census population size
Nc_nushagak <- 600   # Nushagak census population size

# Togiak total census size (MCH + Nushagak)
Nc_togiak <- Nc_mch + Nc_nushagak

# Effective population sizes (assume 20% for MCH and 30% for Togiak)
Ne_mch <- 0.2 * Nc_mch
Ne_togiak <- 0.3 * Nc_togiak

# Pairwise harmonic mean for MCH-Togiak effective population size
Ne_mch_togiak_small <- 2 * (0.1 * Ne_mch) * (0.1 * Ne_togiak) / ((0.1 * Ne_mch) + (0.1 * Ne_togiak))
Ne_mch_togiak_large <- 2 * (0.3 * Ne_mch) * (0.3 * Ne_togiak) / ((0.3 * Ne_mch) + (0.3 * Ne_togiak))

t <- c(0, 10, 30, 50)  # Time in generations: 0, 10, 30, and 50

# Observed pairwise genetic differentiation at time 0
FST0_mch_togiak <- 0.003

# Calculate FST(t) for small and large Ne scenarios
FST_mch_togiak_small <- fst_drift(t, FST0_mch_togiak, Ne_mch_togiak_small)
FST_mch_togiak_large <- fst_drift(t, FST0_mch_togiak, Ne_mch_togiak_large)
FST_mch_togiak_avg <- (FST_mch_togiak_small + FST_mch_togiak_large) / 2

# Combine results into a data frame for plotting
results_mch_togiak <- data.frame(
  Generation = t,
  FST_Min = FST_mch_togiak_small,
  FST_Max = FST_mch_togiak_large,
  FST_Avg = FST_mch_togiak_avg
)

# Plot FST(t) over time for MCH-Togiak comparison
ggplot(results_mch_togiak, aes(x = Generation)) +
  geom_ribbon(aes(ymin = FST_Min, ymax = FST_Max), fill = "#b99c76", alpha = 0.3) +
   geom_line(aes(y = FST_Avg), color = "#b99c76", size = 1) +
  # geom_line(aes(y = FST_Min), color = "#b99c76", size = 0.8) +
  # geom_line(aes(y = FST_Max), color = "#b99c76", size = 0.8) +
  labs(title = "Genetic Differentiation Between MCH and Togiak",
       x = "Generation (t)",
       y = "FST(t)") +
  ylim(0, 0.1) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),  # Remove background grid
    panel.background = element_rect(fill = "white")  # Set background to white
  )

# Combine all results into a single table
results_combined <- rbind(
  data.frame(Comparison = "Central-West", results_central_west),
  data.frame(Comparison = "Central-East", results_central_east),
  data.frame(Comparison = "East-West", results_east_west),
  data.frame(Comparison = "MCH-Togiak", results_mch_togiak
))

# Display the combined table
print(results_combined)

# Save to a CSV file
write.csv(results_combined, "C:/Users/andre/Desktop/FST_Combined_Results.csv", row.names = FALSE)
setwd("C:/Users/andre/Box/New folder")

View(results_combined)



#####################################

# Install necessary packages if not already installed
if (!require("kableExtra")) install.packages("kableExtra")

# Load the library
library(kableExtra)

# Create a pretty table with custom headings and black column lines
custom_headings <- c("Comparison", "Generation (t)", expression(F[ST]~"Min"), expression(F[ST]~"Max"), expression(F[ST]~"Average"))

results_combined %>%
  kbl(col.names = custom_headings, caption = "FST Results for Pairwise Comparisons") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE) %>%
  column_spec(1, border_left = TRUE, border_right = TRUE) %>%
  column_spec(2:5, border_right = TRUE) %>%
  pack_rows("Central-West", 1, 4, hline_after = TRUE) %>%
  pack_rows("Central-East", 5, 8, hline_after = TRUE) %>%
  pack_rows("East-West", 9, 12, hline_after = TRUE) %>%
  pack_rows("MCH-Togiak", 13, 16, hline_after = TRUE)

##########################################################

#incorporating migration 
setwd("C:/Users/andre/Box/New folder")

# Function to calculate FST(t) with migration (m)
fst_drift <- function(t, FST0, Ne, m) {
  FST_t <- 1 - (1 - FST0) * exp(-t / (2 * Ne * (1 - m)))
  return(FST_t)
}

# Parameters
Nc <- 13000     # Census population size

# Subgroup effective population sizes
Ne_central <- 0.15 * Nc  # Central subgroup: 15% of total population
Ne_west <- 0.50 * Nc     # West subgroup: 50% of total population
Ne_east <- 0.35 * Nc     # East subgroup: 35% of total population

# Migration rate (12% 0.12, 0.06, 0.01, 0.00)
m <- 0.06  

# Pairwise harmonic means for effective population size (small and large Ne scenarios)
Ne_central_west_small <- 2 * (0.1 * Ne_central) * (0.1 * Ne_west) / ((0.1 * Ne_central) + (0.1 * Ne_west))
Ne_central_west_large <- 2 * (0.3 * Ne_central) * (0.3 * Ne_west) / ((0.3 * Ne_central) + (0.3 * Ne_west))

Ne_central_east_small <- 2 * (0.1 * Ne_central) * (0.1 * Ne_east) / ((0.1 * Ne_central) + (0.1 * Ne_east))
Ne_central_east_large <- 2 * (0.3 * Ne_central) * (0.3 * Ne_east) / ((0.3 * Ne_central) + (0.3 * Ne_east))

Ne_east_west_small <- 2 * (0.1 * Ne_east) * (0.1 * Ne_west) / ((0.1 * Ne_east) + (0.1 * Ne_west))
Ne_east_west_large <- 2 * (0.3 * Ne_east) * (0.3 * Ne_west) / ((0.3 * Ne_east) + (0.3 * Ne_west))

t <- c(0, 10, 30, 50)  # Time in generations: 0, 10, 30, and 50

# Observed pairwise genetic differentiation at time 0
FST0_central_west <- 0.0008
FST0_central_east <- 0.0004
FST0_east_west <- 0.0000

# Calculate FST(t) with migration (m = 12%)
FST_central_west_small <- fst_drift(t, FST0_central_west, Ne_central_west_small, m)
FST_central_west_large <- fst_drift(t, FST0_central_west, Ne_central_west_large, m)
FST_central_west_avg <- (FST_central_west_small + FST_central_west_large) / 2

FST_central_east_small <- fst_drift(t, FST0_central_east, Ne_central_east_small, m)
FST_central_east_large <- fst_drift(t, FST0_central_east, Ne_central_east_large, m)
FST_central_east_avg <- (FST_central_east_small + FST_central_east_large) / 2

FST_east_west_small <- fst_drift(t, FST0_east_west, Ne_east_west_small, m)
FST_east_west_large <- fst_drift(t, FST0_east_west, Ne_east_west_large, m)
FST_east_west_avg <- (FST_east_west_small + FST_east_west_large) / 2

# Combine results into data frames for plotting
results_central_west <- data.frame(
  Generation = t,
  FST_Min = FST_central_west_small,
  FST_Max = FST_central_west_large,
  FST_Avg = FST_central_west_avg
)

results_central_east <- data.frame(
  Generation = t,
  FST_Min = FST_central_east_small,
  FST_Max = FST_central_east_large,
  FST_Avg = FST_central_east_avg
)

results_east_west <- data.frame(
  Generation = t,
  FST_Min = FST_east_west_small,
  FST_Max = FST_east_west_large,
  FST_Avg = FST_east_west_avg
)

# Load ggplot2 for visualization
library(ggplot2)

# Central-West plot
ggplot(results_central_west, aes(x = Generation)) +
  geom_ribbon(aes(ymin = FST_Min, ymax = FST_Max), fill = "#3f8eaa", alpha = 0.3) +
  geom_line(aes(y = FST_Avg), color = "#3f8eaa", size = 1) +
  labs(title = "Change in Genetic Differentiation: Central-West (12% Migration)",
       x = "Generation (t)",
       y = "FST(t)") +
  ylim(0, 0.1) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),  # Remove background grid
    panel.background = element_rect(fill = "white")  # Set background to white
  )

# Central-East plot
ggplot(results_central_east, aes(x = Generation)) +
  geom_ribbon(aes(ymin = FST_Min, ymax = FST_Max), fill = "#99cd91", alpha = 0.3) +
  geom_line(aes(y = FST_Avg), color = "#99cd91", size = 1) +
  labs(title = "Change in Genetic Differentiation: Central-East (12% Migration)",
       x = "Generation (t)",
       y = "FST(t)") +
  ylim(0, 0.1) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),  # Remove background grid
    panel.background = element_rect(fill = "white")  # Set background to white
  )

# East-West plot
ggplot(results_east_west, aes(x = Generation)) +
  geom_ribbon(aes(ymin = FST_Min, ymax = FST_Max), fill = "#b99c76", alpha = 0.3) +
  geom_line(aes(y = FST_Avg), color = "#b99c76", size = 1) +
  labs(title = "Change in Genetic Differentiation: East-West (12% Migration)",
       x = "Generation (t)",
       y = "FST(t)") +
  ylim(0, 0.1) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),  # Remove background grid
    panel.background = element_rect(fill = "white")  # Set background to white
  )
print(results_east_west)
