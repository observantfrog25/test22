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
source("Road cleaning.R")

##### DON'T RUN ###############
# ca_list <- list(
#   "CA_13",
#   'CA_14',
#   "CA_15",
#   "CA_17",
#   "CA_18",
#   "CA_19",
#   "CA_20",
#   "CA_22",
#   "CA_23",
#   "CA_24"
# )
#
# reproj_fun <- function(abbr) {
#   temp <- sf::read_sf(paste0(
#     here(),
#     "/Roads/HPMS processed/",
#     abbr,
#     "/",
#     abbr,
#     "_processed.shp"
#   ))
#   temp <- sf::st_transform(temp, crs = "epsg:4087")
#
#   sf::write_sf(
#     temp,
#     paste0(
#       here(),
#       "/Roads/HPMS processed/",
#       abbr,
#       "_processed.shp"
#     )
#   )
# }
#
# for (abbr in ca_list) {
#   reproj_fun(abbr)
# }

# 03/30/26 we've decided to go ahead and average out the AADT between years since it's 
# similar
# ca_list <- list(
#   'CA_14',
#   "CA_15",
#   "CA_17",
#   "CA_18",
#   "CA_19",
#   "CA_20",
#   "CA_22",
#   "CA_23"
# )
# ca_major_list <- list()
# 
# for (net in ca_list){
#   ca_major_list[[net]] <- sf::read_sf(paste0(here::here(), "/Roads/CA state lines/", net, "_clean.shp"))
# }
# 
# ca_major_all <- dplyr::bind_rows(ca_major_list)
# 
# avg_aadt <- ca_major_all %>%
#   dplyr::group_by(Direction, Shape_Leng, geometry) %>%
#   dplyr::summarise(
#     AADT = mean(AADT, na.rm = TRUE)
#   ) %>%
#   dplyr::ungroup()
# 
# sf::write_sf(avg_aadt, paste0(here::here(), "/Roads/CA state lines/CA_avg_aadt_clean.shp"))

########## NETWORK CREATON ###########
# we'll eventually want to generalize this further when we determine how we want to do OR and WA
CA_grid <- terra::rast(paste0(here(), "/California5km_raster.tif"))
ca <- readRDS('CA aegypti pseudo.rds')
ca_boundary <- tigris::states(cb = FALSE, resolution = "500k") %>%
    dplyr::filter(NAME == "California") %>%
  terra::vect() %>%
  terra::project(CA_grid)

ca_major <- sf::read_sf(paste0(here::here(), "/Roads/CA state lines/CA_avg_aadt_clean.shp"))

ca_graph <- build_road_network(ca_major, directed = FALSE)

t_flow_fun <- function(years,
                       write = TRUE,
                       write_dir = tempdir()) {
  
  if (length(years) != 2) {
    stop("Years must be of length 2 :(")
  }
  
  ca_pa <- ca[species == "Ae aegypti" &
                year %in% years, .(p = as.integer(any(total > 0))), by = c("longitude", "latitude", "year")]
  
  pts <- terra::vect(ca_pa,
                     geom = c("longitude", "latitude"),
                     crs = "epsg:4326")
  pts <- terra::project(pts, CA_grid) # project it back 4087 (meters) to match everything else 
  
  pts_o <- pts[pts$year == years[1], ] # grab origin nodes (year = i - 1);
  # so if we want to model 2014 presence probability, pts_o is 2013 presence 
  pts_d <- pts[pts$year == years[2], ] # grab destination nodes (year = i)
  
  r_o <- terra::rasterize(
    pts_o,
    CA_grid,
    field = "p",
    fun = max,
    background = NA
  )
  
  r_d <- terra::rasterize(
    pts_d,
    CA_grid,
    field = "p",
    fun = max,
    background = NA
  )
  
  r_o <- terra::mask(r_o, ca_boundary)
  r_d <- terra::mask(r_d, ca_boundary)
  
  # get cell indices for presence year i - 1 and all collections in year i
  cells_o <- which(values(r_o) == 1)
  cells_d <- which(values(r_d) %in% c(0, 1))
  
  edges <- expand.grid(from = cells_o, to = cells_d) %>%
    subset(from != to) # remove self-loops
  
  nodes <- ca_graph %>%
    activate(nodes) %>%
    as_tibble() %>%
    st_as_sf()
  
  coords <- nodes %>%
    st_coordinates()
  
  cell_metrics <- lapply(cells_d, function(target) {
    connected <- edges$from[edges$to == target]
    
    if (length(connected) == 0) {
      return(data.frame(aadt = NA, dist_d = NA))
    }
    
    coords_o <- terra::xyFromCell(r_o, connected)
    coord_d  <- terra::xyFromCell(r_d, target)
    
    node_indices_o <- knn(data = coords, 
                          query = coords_o, 
                          k = 1) 
    node_index_d   <- knn(data = coords, 
                          query = matrix(coord_d, ncol = 2), 
                          k = 1)
    node_d <- nodes[node_index_d$nn.idx, ]
    
    aadt_vals <- sapply(seq_along(node_indices_o$nn.idx), function(i) {
      idx <- node_indices_o$nn.idx[i] 
      
      node_o <- nodes[idx, ]
      path <- shortest_paths(
        graph = ca_graph,
        from = node_o$nodeID,
        to = node_d$nodeID,
        output = 'both',
        weights = ca_graph %>% activate(edges) %>% pull(Shape_Leng)
      ) 
      
      pg <- ca_graph %>%
        subgraph_from_edges(eids = path$epath %>% unlist()) %>%
        as_tbl_graph()
      
      mean_aadt <- pg %>%
        activate(edges) %>%
        as_tibble() %>%
        summarise(mean_aadt = mean(AADT, na.rm = TRUE)) %>%
        pull(mean_aadt) 
      
      return(mean_aadt)
    })
    
    # AS OF 03/31/26 NO distance penalty
    raw_exposure <- sum(aadt_vals, na.rm = TRUE)
    
    # destination distance
    dist_d <- node_index_d$nn.dists[1, 1]
    
    return(data.frame(aadt = raw_exposure, dist_d = dist_d))
  })
  
  cell_metrics_df <- dplyr::bind_rows(cell_metrics)
  
  # create raster for AADT
  r_aadt_d <- CA_grid 
  values(r_aadt_d) <- NA 
  r_aadt_d[cells_d] <- cell_metrics_df$aadt
  r_aadt_d <- terra::mask(r_aadt_d, ca_boundary)
  
  # create raster for Destination Distances
  r_dist_d <- CA_grid
  values(r_dist_d) <- NA
  r_dist_d[cells_d] <- cell_metrics_df$dist_d
  r_dist_d <- terra::mask(r_dist_d, ca_boundary)
  
  # stack them into a multi-layer SpatRaster
  out_raster <- c(r_aadt_d, r_dist_d)
  names(out_raster) <- c("AADT", "Dist_to_Road")
  
  if (write) {
    terra::writeRaster(out_raster, filename = write_dir)
  }
  
  return(out_raster)
}
