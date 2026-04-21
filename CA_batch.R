library(terra)
library(igraph)
library(nabor)
library(here)
library(sf)
library(data.table)
library(dplyr)
library(tidygraph)
library(sfnetworks)
library(readr)
library(psych)
library(future.apply)

source("Road cleaning.R")

plan(multisession, workers = parallelly::availableCores())

sf_to_tidygraph <- function(x, directed = FALSE, snap_tolerance = 100) {
  x <- st_cast(x, "LINESTRING")
  x <- st_make_valid(x)
  x <- st_set_precision(x, 1 / snap_tolerance)

  net <- sfnetworks::as_sfnetwork(x, directed = directed)

  list_col_names <- function(tbl) {
    nms <- names(tbl)
    nms[vapply(tbl, is.list, logical(1)) & nms != attr(tbl, "sf_column")]
  }

  node_drops <- list_col_names(sf::st_as_sf(net, "nodes"))
  edge_drops <- list_col_names(sf::st_as_sf(net, "edges"))

  if (length(node_drops) > 0) {
    net <- net %>% activate(nodes) %>% select(-all_of(node_drops))
  }
  if (length(edge_drops) > 0) {
    net <- net %>% activate(edges) %>% select(-all_of(edge_drops))
  }

  as_tbl_graph(net)
}

t_flow_fun <- function(years,
                       write = TRUE,
                       write_dir = tempdir()) {
  if (length(years) != 2) {
    stop("Years must be of length 2 :(")
  }

  ca_pa <- ca[species == "Ae aegypti" &
                year %in% years, .(p = as.integer(any(total > 0))),
              by = c("longitude", "latitude", "year")]

  pts <- terra::vect(ca_pa,
                     geom = c("longitude", "latitude"),
                     crs = "epsg:4326")
  pts <- terra::project(pts, CA_grid)

  pts_o <- pts[pts$year == years[1], ]
  pts_d <- pts[pts$year == years[2], ]

  r_o <- terra::rasterize(pts_o, CA_grid, field = "p", fun = max, background = NA)
  r_d <- terra::rasterize(pts_d, CA_grid, field = "p", fun = max, background = NA)

  r_o <- terra::mask(r_o, ca_boundary)
  r_d <- terra::mask(r_d, ca_boundary)

  cells_o <- which(values(r_o) == 1)
  cells_d <- which(values(r_d) %in% c(0, 1))

  edges <- expand.grid(from = cells_o, to = cells_d) %>%
    subset(from != to)

  nodes <- ca_graph %>%
    activate(nodes) %>%
    as_tibble() %>%
    st_as_sf()

  coords <- st_coordinates(nodes)
  edge_weights <- ca_graph %>% activate(edges) %>% pull(Shape_Leng)

  r_o_wrapped <- terra::wrap(r_o)
  r_d_wrapped <- terra::wrap(r_d)

  cell_metrics <- future_lapply(cells_d, function(target) {
    r_o_local <- terra::unwrap(r_o_wrapped)
    r_d_local <- terra::unwrap(r_d_wrapped)

    connected <- edges$from[edges$to == target]

    if (length(connected) == 0) {
      return(data.frame(mean_aadt = NA, median_aadt = NA,
                        harmonic_aadt = NA, dist_d = NA))
    }

    coords_o <- terra::xyFromCell(r_o_local, connected)
    coord_d  <- terra::xyFromCell(r_d_local, target)

    node_indices_o <- nabor::knn(data = coords, query = coords_o, k = 1)
    node_index_d   <- nabor::knn(data = coords,
                                 query = matrix(coord_d, ncol = 2), k = 1)

    origin_ids <- node_indices_o$nn.idx
    dest_id    <- node_index_d$nn.idx[1]

    aadt_vals <- sapply(origin_ids, function(o_idx) {
      path <- igraph::shortest_paths(
        graph   = ca_graph,
        from    = o_idx,
        to      = dest_id,
        output  = "both",
        weights = edge_weights
      )
      epath <- path$epath[[1]]
      if (length(epath) > 0) {
        edge_aadt <- igraph::edge_attr(ca_graph, "AADT", index = epath)
        scaled    <- edge_aadt / 1e6
        c(mean     = mean(scaled, na.rm = TRUE),
          median   = median(scaled, na.rm = TRUE),
          harmonic = psych::harmonic.mean(scaled, na.rm = TRUE, zero = FALSE))
      } else {
        c(mean = 0, median = 0, harmonic = 0)
      }
    })

    dist_d <- node_index_d$nn.dists[1, 1]

    if (is.matrix(aadt_vals)) {
      data.frame(
        mean_aadt     = sum(aadt_vals["mean", ],     na.rm = TRUE),
        median_aadt   = sum(aadt_vals["median", ],   na.rm = TRUE),
        harmonic_aadt = sum(aadt_vals["harmonic", ], na.rm = TRUE),
        dist_d        = dist_d
      )
    } else {
      data.frame(
        mean_aadt     = aadt_vals["mean"],
        median_aadt   = aadt_vals["median"],
        harmonic_aadt = aadt_vals["harmonic"],
        dist_d        = dist_d
      )
    }
  }, future.seed = TRUE)

  cell_metrics_df <- dplyr::bind_rows(cell_metrics)

  make_layer <- function(vals) {
    r <- CA_grid
    values(r) <- NA
    r[cells_d] <- vals
    terra::mask(r, ca_boundary)
  }

  out_raster <- c(
    make_layer(cell_metrics_df$mean_aadt),
    make_layer(cell_metrics_df$median_aadt),
    make_layer(cell_metrics_df$harmonic_aadt),
    make_layer(cell_metrics_df$dist_d)
  )
  names(out_raster) <- c("AADT_Mean", "AADT_Median", "AADT_Harmonic", "Dist_to_Road")

  if (write) {
    terra::writeRaster(out_raster, filename = write_dir, overwrite = TRUE)
  }

  return(out_raster)
}

##### SLURM SETUP #####
yr1_sequence <- 2013:2024
yr2_sequence <- 2014:2025

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) task_id <- 1

yr1 <- yr1_sequence[task_id]
yr2 <- yr2_sequence[task_id]

CA_grid <- terra::rast(here::here("California5km_raster.tif"))
ca <- readRDS(here::here("CA aegypti pseudo.rds"))

ca_boundary <- readRDS(here::here("ca_boundary.rds")) %>%
  terra::vect() %>%
  terra::project(CA_grid)

ca_major <- sf::read_sf(here::here("CA_avg_aadt_clean.shp"))
ca_graph <- sf_to_tidygraph(ca_major)

print(paste("Task ID:", task_id, "| Processing Years:", yr1, "and", yr2))

if (!dir.exists(here::here("results")))
  dir.create(here::here("results"))

output_path <- here::here("results", paste0("flow_output_", yr1, "_", yr2, ".tif"))

temp <- t_flow_fun(
  years     = c(yr1, yr2),
  write     = TRUE,
  write_dir = output_path
)

##### WILCOXON TEST #####
ca_pa <- ca[species == "Ae aegypti" &
              year == yr2, .(p = as.integer(any(total > 0))),
            by = c("longitude", "latitude", "year")]

pts <- terra::vect(ca_pa,
                   geom = c("longitude", "latitude"),
                   crs = "epsg:4326")
pts <- terra::project(pts, CA_grid)

pts_d <- pts[pts$year == yr2, ]

r_d <- terra::rasterize(pts_d, CA_grid, field = "p", fun = max, background = NA)

pa_aadt <- data.frame(
  aadt = terra::values(temp[["AADT_Mean"]]),
  max  = terra::values(r_d)
) %>%
  tidyr::drop_na()

wil <- wilcox.test(
  subset(pa_aadt, max == 1)$aadt,
  subset(pa_aadt, max == 0)$aadt
)
print(wil)

wil_path <- here::here("results", paste0("wilcox_", yr1, "_", yr2, ".rds"))
readr::write_rds(wil, wil_path)
