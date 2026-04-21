library(terra)
library(igraph)
library(nabor)
library(here)
library(sf)
library(data.table)
library(dplyr)
library(tidygraph)
library(sfnetworks)
library(psych)
library(future.apply)

source("Road cleaning.R")

plan(multisession, workers = parallelly::availableCores())
if (!dir.exists("results/chunks/WA_para_testing/")) {
  dir.create("results/chunks/WA_para_testing/", recursive = TRUE)
}

t_flow_chunked <- function(year,
                           dest_grid,
                           dest_boundary,
                           task_id,
                           num_chunks,
                           save_dir = "results/chunks/WA_para_testing/") {
  all_roads <- sf::read_sf(here::here("CA_WA_OR_roads.shp"))

  full_graph <- build_road_network(all_roads, directed = FALSE, snap_tolerance = 100) %>%
    activate(nodes) %>%
    mutate(component = group_components()) %>%
    filter(component == names(which.max(table(component))))

  ca_pa_origin <- ca[species == "Ae aegypti" & year == year,
                     .(p = as.integer(any(total > 0))),
                     by = c("longitude", "latitude")]
  pts_o <- terra::project(
    terra::vect(ca_pa_origin, geom = c("longitude", "latitude"), crs = "epsg:4326"),
    "epsg:4087"
  )
  r_o <- terra::mask(
    terra::rasterize(pts_o, CA_grid, field = "p", fun = max),
    ca_boundary
  )

  r_d <- terra::mask(dest_grid, dest_boundary)

  all_cells_d <- which(!is.na(values(r_d)))

  chunk_size <- ceiling(length(all_cells_d) / num_chunks)
  start_idx  <- (task_id - 1) * chunk_size + 1
  end_idx    <- min(task_id * chunk_size, length(all_cells_d))

  my_cells <- all_cells_d[start_idx:end_idx]

  nodes  <- full_graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf()
  coords <- st_coordinates(nodes)
  coords_o <- terra::xyFromCell(r_o, which(values(r_o) == 1))
  origin_node_ids <- nabor::knn(data = coords, query = coords_o, k = 1)$nn.idx
  edge_weights <- full_graph %>% activate(edges) %>% pull(Shape_Leng)

  r_d_wrapped <- terra::wrap(r_d)

  cell_metrics <- future_lapply(my_cells, function(target) {
    r_d_local <- terra::unwrap(r_d_wrapped)

    coord_d <- terra::xyFromCell(r_d_local, target)
    node_index_d <- nabor::knn(
      data  = coords,
      query = matrix(coord_d, ncol = 2),
      k     = 1
    )

    aadt_vals <- sapply(origin_node_ids, function(o_node) {
      path <- igraph::shortest_paths(
        full_graph,
        from    = o_node,
        to      = node_index_d$nn.idx,
        output  = "both",
        weights = edge_weights
      )
      epath <- path$epath[[1]]
      if (length(epath) > 0) {
        edge_aadt <- igraph::edge_attr(full_graph, "AADT", index = epath)
        scaled    <- edge_aadt / 1e6
        c(mean     = mean(scaled, na.rm = TRUE),
          median   = median(scaled, na.rm = TRUE),
          harmonic = psych::harmonic.mean(scaled, na.rm = TRUE, zero = FALSE))
      } else {
        c(mean = 0, median = 0, harmonic = 0)
      }
    })

    # aadt_vals is a 3 x n_origins matrix; sum across origins
    if (is.matrix(aadt_vals)) {
      data.frame(
        mean_aadt     = sum(aadt_vals["mean", ],     na.rm = TRUE),
        median_aadt   = sum(aadt_vals["median", ],   na.rm = TRUE),
        harmonic_aadt = sum(aadt_vals["harmonic", ], na.rm = TRUE)
      )
    } else {
      data.frame(mean_aadt = 0, median_aadt = 0, harmonic_aadt = 0)
    }
  }, future.seed = TRUE)

  res_df <- dplyr::bind_rows(cell_metrics)
  save_name <- paste0(save_dir, "WA_", year, "_chunk_", task_id, ".rds")
  saveRDS(list(cells = my_cells, metrics = res_df), save_name)

  return(paste("Task", task_id, "completed and saved."))
}

task <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
total_chunks <- 50

WA_grid <- terra::rast(paste0(here(), "/Washington5km_raster.tif"))

wa_boundary <- readRDS("wa_boundary.rds") %>%
  terra::vect() %>%
  terra::project(WA_grid)

CA_grid <- terra::rast(here::here("California5km_raster.tif"))
ca <- readRDS(here::here("CA aegypti pseudo.rds"))

ca_boundary <- readRDS(here::here("ca_boundary.rds")) %>%
  terra::vect() %>%
  terra::project(CA_grid)

t_flow_chunked(
  year         = 2024,
  dest_grid    = WA_grid,
  dest_boundary = wa_boundary,
  task_id      = task,
  num_chunks   = total_chunks
)
