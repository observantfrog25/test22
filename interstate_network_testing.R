library(terra)
library(igraph)
library(nabor)
library(here)
library(sf)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyterra)
library(tidygraph)
library(dbscan)
library(sfnetworks)
library(stplanr)
source("Road cleaning.R")
t_flow_interstate <- function(year,
                              dest_grid,
                              dest_boundary,
                              write = FALSE,
                              write_dir = tempdir()) {
  all_roads <- sf::read_sf(here::here("Roads", "CA_WA_OR_snapped_clean.shp"))
  full_graph <- build_road_network(all_roads, directed = FALSE)
  
  ca_pa_origin <- ca[species == "Ae aegypti" &
                       year == year, .(p = as.integer(any(total > 0))), by = c("longitude", "latitude")]
  
  pts_o <- terra::vect(ca_pa_origin,
                       geom = c("longitude", "latitude"),
                       crs = "epsg:4326")
  pts_o <- terra::project(pts_o, "epsg:4087")
  
  r_o <- terra::rasterize(
    pts_o,
    CA_grid,
    field = "p",
    fun = max,
    background = NA
  )
  r_o <- terra::mask(r_o, ca_boundary)
  
  r_d <- dest_grid
  values(r_d) <- 1
  r_d <- terra::mask(r_d, dest_boundary)
  
  cells_o <- which(values(r_o) == 1)
  cells_d <- which(!is.na(values(r_d)))
  
  nodes <- full_graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf()
  coords <- st_coordinates(nodes)
  coords_o <- terra::xyFromCell(r_o, cells_o)
  
  node_indices_o <- nabor::knn(data = coords, query = coords_o, k = 1)
  origin_node_ids <- node_indices_o$nn.idx # A vector of IDs
  
  edge_weights <- full_graph %>% activate(edges) %>% pull(Shape_Leng)
  
  cell_metrics <- lapply(cells_d, function(target) {
    coord_d <- terra::xyFromCell(r_d, target)
    node_index_d <- nabor::knn(data = coords,
                               query = matrix(coord_d, ncol = 2),
                               k = 1)
    target_node_id <- node_index_d$nn.idx
    
    aadt_vals <- sapply(origin_node_ids, function(o_node) {
      path <- igraph::shortest_paths(
        graph = full_graph,
        from = o_node,
        to = target_node_id,
        output = 'both',
        weights = edge_weights
      )
      
      epath <- path$epath[[1]]
      if (length(epath) > 0) {
        edge_aadt <- igraph::edge_attr(full_graph, "AADT", index = epath)
        return(mean(edge_aadt / (10^6), na.rm = TRUE))
      } else {
        return(0)
      }
    })
    
    return(data.frame(
      aadt = sum(aadt_vals, na.rm = TRUE),
      dist_d = node_index_d$nn.dists[1, 1]
    ))
  })
  
  cell_metrics_df <- dplyr::bind_rows(cell_metrics)
  
  r_aadt_d <- r_d
  values(r_aadt_d) <- NA
  r_aadt_d[cells_d] <- cell_metrics_df$aadt
  
  r_dist_d <- r_d
  values(r_dist_d) <- NA
  r_dist_d[cells_d] <- cell_metrics_df$dist_d
  
  out_raster <- c(r_aadt_d, r_dist_d)
  names(out_raster) <- c("AADT", "Dist_to_Road")
  
  if (write)
    terra::writeRaster(out_raster, filename = write_dir, overwrite = TRUE)
  
  return(out_raster)
}

CA_grid <- terra::rast(paste0(here(), "/California5km_raster.tif"))
ca <- readRDS('CA aegypti pseudo.rds')
ca_boundary <- tigris::states(cb = FALSE, resolution = "500k") %>%
  dplyr::filter(NAME == "California") %>%
  terra::vect() %>%
  terra::project(CA_grid)
#
# test <- t_flow_interstate(
#   year = 2013,
#   dest_grid = OR_grid_subset,
#   dest_boundary = or_boundary_subset,
#   write = FALSE
# )
#
# test <- readRDS('OR_2024_chunk_6.rds')
#
# test_raster <- OR_grid * NA
# test_raster[test$cells] <- test$metrics$aadt
#
#
# ggplot() +
#   geom_spatraster(data = test_raster[[1]]) +
#   scale_fill_viridis_c(name = "Sum AADT", na.value = "white") +
#   geom_spatvector(
#     data        = or_boundary,
#     inherit.aes = FALSE,
#     fill        = NA,
#     color       = "black",
#     linewidth   = 0.5
#   ) +
#   theme_minimal() # LET'S FUCKING GOOOOOOOOOOOOOOOOOOOOOOOOOO (works 04/01)

# all chunks are finished for OR... testing that they work
WA_grid <- terra::rast(paste0(here(), "/Washington5km_raster.tif"))
wa_boundary <- readRDS("wa_boundary.rds") %>%
  terra::vect() %>%
  terra::project(WA_grid)

final_aadt_rast <- WA_grid * NA
final_dist_rast <- WA_grid * NA

chunk_files <- list.files(
  here::here("Roads/traffic_flow_extracted/WA_alt"),
  pattern = "WA_2024_chunk_.*\\.rds",
  full.names = TRUE
)

for (f in chunk_files) {
  chunk <- readRDS(f)
  
  
  final_aadt_rast[chunk$cells] <- chunk$metrics$aadt
  final_dist_rast[chunk$cells] <- chunk$metrics$dist_d
}

wa_exposure <- c(final_aadt_rast, final_dist_rast)
names(wa_exposure) <- c("AADT_Exposure", "Dist_to_Road")

# wa_exposure <- terra::crop(wa_exposure, wa_boundary)
# wa_exposure <- terra::mask(wa_exposure, wa_boundary)

wa_roads <- build_road_network(sf::read_sf(paste0(
  here::here(), "/Roads/WA state lines/WA_24_clean.shp"
)), directed = FALSE)

ggplot() +
  geom_spatraster(data = wa_exposure[[1]]) +
  scale_fill_viridis_c(name = "Sum AADT (scaled by a factor of 1e^-6)", na.value = "white") +
  geom_spatvector(
    data        = wa_boundary,
    inherit.aes = FALSE,
    fill        = NA,
    color       = "black",
    linewidth   = 0.5
  ) +
  theme_minimal()

ggsave(
  filename = paste0(here(), "/Plots/AADT examples/WA_aadt.png"),
  bg = "white",
  height = 11,
  width = 13,
  units = "in",
  dpi = 300
)

ggplot() +
  geom_spatraster(data = wa_exposure[[2]], alpha = 0.75) +
  scale_fill_viridis_c(name = "Distance to nearest neighbor (m)", na.value = "white") +
  geom_spatvector(
    data        = wa_boundary,
    inherit.aes = FALSE,
    fill        = NA,
    color       = "black",
    linewidth   = 0.5
  ) +
  geom_sf(
    data = wa_roads %>% activate(edges) %>% as_tibble() %>% st_as_sf(),
    col = 'darkgrey'
  ) +
  geom_sf(
    data = wa_roads %>% activate(nodes) %>% as_tibble() %>% st_as_sf(),
    col = 'firebrick',
    size = 0.15
  ) +
  theme_minimal()


ggsave(
  filename = paste0(here(), "/Plots/AADT examples/WA_knn_dist.png"),
  bg = "white",
  height = 11,
  width = 13,
  units = "in",
  dpi = 300
)

ggplot() +
  geom_spatraster(data = wa_exposure[[1]]) +
  scale_fill_viridis_c(name = "", na.value = "white") +
  geom_spatvector(
    data        = wa_boundary,
    inherit.aes = FALSE,
    fill        = NA,
    color       = "black",
    linewidth   = 0.5
  ) +
  geom_sf(
    data = full_graph %>% activate(edges) %>% as_tibble() %>% st_as_sf(),
    col = 'darkgrey'
  ) +
  geom_sf(
    data = full_graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf(),
    col = 'firebrick',
    size = 0.05
  ) +
  theme_minimal()

ggsave(
  filename = paste0(here(), "/Plots/AADT examples/WA_aadt_with_road.png"),
  bg = "white",
  height = 11,
  width = 13,
  units = "in",
  dpi = 300
)

library(mapview)

m <- mapview(wa_exposure[[1]],
        layer.name = "Exposure",
        col.regions = viridis::viridis(100)) +
  mapview(
    full_graph %>% activate(edges) %>% as_tibble() %>% st_as_sf(),
    color = "darkgrey",
    legend = FALSE
  ) +
  mapview(
    full_graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf(),
    color = "firebrick",
    cex = 0.45,
    legend = FALSE
  )
mapshot(m, url = "WA_AADT.html")


# testing out why AADT differs in clustered cells
all_roads <- sf::read_sf(here::here("Roads", "CA_WA_OR_avg_aadt_clean.shp"))
full_graph <- build_road_network(all_roads, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(component = group_components()) %>%
  filter(component == names(which.max(table(component))))

ca_pa_origin <- ca[species == "Ae aegypti" &
                     year == year, .(p = as.integer(any(total > 0))), by = c("longitude", "latitude")]

pts_o <- terra::vect(ca_pa_origin,
                     geom = c("longitude", "latitude"),
                     crs = "epsg:4326")
pts_o <- terra::project(pts_o, "epsg:4087")

r_o <- terra::rasterize(pts_o,
                        CA_grid,
                        field = "p",
                        fun = max,
                        background = NA)
r_o <- terra::mask(r_o, ca_boundary)

r_d <- WA_grid
values(r_d) <- 1
r_d <- terra::mask(r_d, wa_boundary)

cells_o <- which(values(r_o) == 1)
cells_d <- r_o[13159]

nodes <- full_graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf()
coords <- st_coordinates(nodes)
coords_o <- terra::xyFromCell(r_o, cells_o)

node_indices_o <- nabor::knn(data = coords, query = coords_o, k = 1)
origin_node_ids <- node_indices_o$nn.idx

edge_weights <- full_graph %>% activate(edges) %>% pull(Shape_Leng)

coord_d <- terra::xyFromCell(r_d, 13336)
node_index_d <- nabor::knn(data = coords,
                           query = matrix(coord_d, ncol = 2),
                           k = 1)

target_node_id <- node_index_d$nn.idx

aadt_vals <- sapply(origin_node_ids, function(o_node) {
  path <- igraph::shortest_paths(
    graph = full_graph,
    from = o_node,
    to = target_node_id,
    output = 'both',
    weights = edge_weights
  )
  
  epath <- path$epath[[1]]
  if (length(epath) > 0) {
    edge_aadt <- igraph::edge_attr(full_graph, "AADT", index = epath)
    return(mean(edge_aadt / (10^6), na.rm = TRUE))
  } else {
    return(0)
  }
})

data.frame(
  aadt = sum(aadt_vals, na.rm = TRUE),
  dist_d = node_index_d$nn.dists[1, 1]
)

cell_metrics_df <- dplyr::bind_rows(cell_metrics)

r_aadt_d <- r_d
values(r_aadt_d) <- NA
r_aadt_d[cells_d] <- cell_metrics_df$aadt

r_dist_d <- r_d
values(r_dist_d) <- NA
r_dist_d[cells_d] <- cell_metrics_df$dist_d

out_raster <- c(r_aadt_d, r_dist_d)
names(out_raster) <- c("AADT", "Dist_to_Road")


all_roads %<%
  transform(AADT_scaled = AADT/(10^6)) |> 
  mapview(zcol = "AADT_scaled", 
          col.regions = mapviewColors("viridis"), 
          lwd = 3)

#### TESTING SEATTLE NODES #####
seattle_pt <- st_point(c(-122.3321, 47.6062)) %>% 
  st_sfc(crs = 4326) %>% 
  st_transform(4087)

all_roads <- sf::read_sf(here::here("Roads", "merge_snap.shp"))

full_graph <- build_road_network(all_roads, directed = FALSE)
# # 
# mapview(full_graph %>%
#           activate(edges) %>%
#           as_tibble() %>%
#           st_as_sf(), lwd = 1, color = "black") +
#   mapview(full_graph %>%
#   activate(nodes) %>%
#   as_tibble() %>%
#   st_as_sf(), color = "firebrick", cex = 0.45,
#   legend = FALSE)

full_graph <- full_graph %>%
  activate(nodes) %>%
  mutate(component = group_components()) %>%
  filter(component == names(which.max(table(component))))


# 
# mapview(full_graph %>%
#           activate(edges) %>%
#           as_tibble() %>%
#           st_as_sf(), lwd = 1, color = "black") +
#   mapview(st_transform(st_as_sf(seattle_pt), 4326), col.regions = "forestgreen", layer.name = "Seattle point") +
#   mapview(
#     full_graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf(),
#     color = "firebrick",
#     cex = 0.45,
#     legend = FALSE
#   )

test_year <- 2013
ca_pa_origin <- ca[species == "Ae aegypti" & year == test_year, 
                   .(p = as.integer(any(total > 0))), 
                   by = c("longitude", "latitude")]

pts_o <- terra::vect(ca_pa_origin[p == 1], geom = c("longitude", "latitude"), crs = "epsg:4326")
pts_o <- terra::project(pts_o, "epsg:4087")
coords_o <- terra::crds(pts_o)

nodes_sf <- full_graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf()
node_coords <- st_coordinates(nodes_sf)


node_indices_o <- nabor::knn(data = node_coords, query = coords_o, k = 1)$nn.idx

node_index_d <- nabor::knn(data = node_coords, query = st_coordinates(seattle_pt), k = 1)$nn.idx

snapped_nodes_o <- nodes_sf %>% 
  slice(as.vector(node_indices_o)) %>%
  mutate(label = "Snapped Origin Node")

snapped_node_d <- nodes_sf %>% 
  slice(as.vector(node_index_d)) %>%
  mutate(label = "Snapped Destination Node")

edge_weights <- full_graph %>% activate(edges) %>% pull(Shape_Leng)

valid_paths <- lapply(node_indices_o, function(o_node) {
  path <- igraph::shortest_paths(
    graph = full_graph,
    from = o_node,
    to = node_index_d,
    weights = edge_weights,
    output = 'both'
  )
  
  edge_ids <- as.numeric(path$epath[[1]])
  
  if (length(edge_ids) > 0) {
    # Convert to tibble FIRST, then to sf
    path_sf <- full_graph %>% 
      activate(edges) %>% 
      slice(edge_ids) %>% 
      as_tibble() %>%    # <--- This is the missing link
      st_as_sf()
    
    return(path_sf)
  }
  return(NULL)
})

all_paths_sf <- do.call(rbind, valid_paths)

mapview(all_paths_sf, color = "firebrick", lwd = 2, layer.name = "Shortest Paths") +
  mapview(pts_o, col.regions = "dodgerblue", layer.name = "CA Original Points") +
  mapview(snapped_nodes_o, col.regions = "cyan", cex = 3, layer.name = "Snapped Road Nodes (CA)") +
  mapview(seattle_pt, col.regions = "forestgreen", layer.name = "Seattle Destination") +
  mapview(snapped_node_d, col.regions = "limegreen", cex = 3, layer.name = "Snapped Road Node (SEA)")
