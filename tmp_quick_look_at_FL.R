library(tidyverse)
nbr_files = dir("output/04_merge/", pattern = "neighbor", full.names = TRUE)
names(nbr_files) = basename(nbr_files)
nbr_files = as.list(nbr_files)
nbr_df.l = lapply(nbr_files, readr::read_csv)

nbr_df = bind_rows(nbr_df.l, .id = "file")
nbr_df = nbr_df %>% mutate(fraction_kappa = n_kappa / n_cells)
nbr_df = nbr_df %>% mutate(name = sub(".+hoods_", "", file))

ggplot(nbr_df, aes(x = n_lambda, y = n_kappa)) +
  geom_point(size = .3) +
  facet_wrap(~name)

ggplot(nbr_df, aes(x = n_lambda, y = n_kappa, color = fraction_kappa)) +
  geom_point(size = .3) +
  facet_wrap(~name) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_viridis_c(option = "magma")


ggplot(nbr_df, aes(x = n_lambda, y = n_kappa, color = fraction_kappa > .7)) +
  geom_point(size = .3) +
  facet_wrap(~name) +
  scale_x_log10() +
  scale_y_log10()

best_df = nbr_df %>% filter(name == "r20um.csv" & fraction_kappa > .7)

cell_file = dir("output/04_merge/", pattern = "cells_r20", full.names = TRUE)
cell_df = readr::read_csv(cell_file)
cell_df %>% filter(k_ct > 0 & l_ct > 0) %>% arrange(-k_ct)
cell_df %>% filter(k_ct > 0 & l_ct > 0) %>% arrange(-l_ct)
cell_df$k_ct %>% sort %>% plot
cell_df$k_ct %>% log10() %>% hist
cell_df$l_ct %>% log10() %>% hist

nbr_df %>% filter(name == "r20um.csv") %>% select(nb, n_kappa, n_lambda, name)
best_df = merge(best_df, cell_df)
ggplot(best_df, aes(x = l_ct+10, y = k_ct+10)) +
  geom_jitter(size = .3) +
  facet_wrap(~nb) +
  scale_x_log10() +
  scale_y_log10()

ggplot(best_df, aes(x = l_ct+10, y = k_ct+10)) +
  geom_density2d_filled() +
  facet_wrap(~nb) +
  scale_x_log10() +
  scale_y_log10()

nbr_df.l$lightchain_neighborhoods_r20um.csv
