library(dplyr)
library(data.table)
highway_extract <- function(location, write_dir, return = TRUE) {
  a <- osmdata::getbb(location) %>%
    osmdata::opq() %>%
    osmdata::add_osm_feature(key = "highway", value = "primary")
  
  major <- osmdata::osmdata_sf(a)
  major <- major$osm_lines %>%
    st_transform(crs = "epsg:4087")
  
  sf::write_sf(major, write_dir)
  
  if (return) {
    return(major)
  }
}

aadt_append <- function(aadt_path, write_dir, ca_lines = ca_major) {
  aadt <- data.table::setDT(readxl::read_excel(aadt_path))
  aadt$ROUTE <- as.integer(aadt$ROUTE)
  aadt$AHEAD_AADT <- as.integer(aadt$AHEAD_AADT)
  aadt$BACK_AADT <- as.integer(aadt$BACK_AADT)
  aadt$PM <- as.double(aadt$PM)
  
  data.table::setnames(aadt,
                       old = c("ROUTE", "COUNTY"),
                       new = c("RTE", "CNTY"))
  
  aadt_clean <- aadt[!is.na(PM) &
                       (!is.na(BACK_AADT) | !is.na(AHEAD_AADT)), .(RTE, CNTY, PM, BACK_AADT, AHEAD_AADT)]
  
  lines_df <- ca_lines %>%
    st_drop_geometry() %>%
    as.data.table()
  
  matched <- aadt_clean[lines_df, on = .(RTE = Route, CNTY = County, PM >= bPM, PM <= ePM), .(
    OBJECTID,
    Route = RTE,
    County = CNTY,
    bPM = i.bPM,
    ePM = i.ePM,
    BACK_AADT,
    AHEAD_AADT
  ), nomatch = 0]
  
  collapsed <- matched[, .(
    BACK_AADT  = mean(BACK_AADT, na.rm = TRUE),
    AHEAD_AADT = mean(AHEAD_AADT, na.rm = TRUE)
  ), by = .(OBJECTID, Route, County, bPM, ePM)]
  
  
  aadt_long <- rbindlist(list(collapsed[, .(
    OBJECTID,
    Route,
    County,
    bPM,
    ePM,
    AADT      = BACK_AADT,
    Direction = fifelse(Route %% 2 == 0, "WB", "SB")
  )], collapsed[, .(
    OBJECTID,
    Route,
    County,
    bPM,
    ePM,
    AADT      = AHEAD_AADT,
    Direction = fifelse(Route %% 2 == 0, "EB", "NB")
  )]))
  
  ca_major_aadt <- ca_lines %>%
    left_join(aadt_long,
              by = c("OBJECTID", "Route", "County" = "County", "bPM", "ePM", "Direction"))
  
  ca_major_aadt$AADT <- ifelse(is.na(ca_major_aadt$AADT), 0, ca_major_aadt$AADT)
  
  ca_major_aadt <- st_transform(ca_major_aadt, "epsg:4087")
  
  sf::write_sf(ca_major_aadt, write_dir)
}

highway_append <- function(NHS_shp,
                           # road shapefile,
                           aadt_shp # aadt data
                           ) {
                           NHS_shp <- sf::st_transform(NHS_shp, "epsg:4087")[c("Route", "County", "bPM", "ePM", "Direction", "Shape_Leng")]
                           
                           aadt_shp <- sf::st_zm(aadt_shp)
                           
                           aadt_shp <- aadt_shp %>%
                             sf::st_drop_geometry()
                           
                           data.table::setDT(aadt_shp)
                           
                           aadt_shp <- aadt_shp[, .(aadt = sum(aadt, na.rm = TRUE)), by = route_id]
                           
                           dplyr::left_join(NHS_shp, aadt_shp, by = "route_id")
                           }

sf_to_tidygraph <- function(x, directed = TRUE, snap_tolerance = 1) {
  edges <- x %>%
    st_cast("LINESTRING") %>%
    mutate(edgeID = row_number())

  coords <- edges %>%
    st_coordinates() %>%
    as_tibble()

  id_col <- if ("L2" %in% names(coords)) "L2" else "L1"
  coords <- coords %>% rename(edgeID = !!sym(id_col))

  endpoints <- coords %>%
    group_by(edgeID) %>%
    slice(c(1, n())) %>%
    ungroup() %>%
    mutate(start_end = rep(c('start', 'end'), times = n() / 2))

  endpoints <- endpoints %>%
    mutate(
      X_snap = round(X / snap_tolerance) * snap_tolerance,
      Y_snap = round(Y / snap_tolerance) * snap_tolerance,
      xy = paste(X_snap, Y_snap)
    ) %>%
    group_by(xy) %>%
    mutate(nodeID = cur_group_id()) %>%
    ungroup()

  source_nodes <- endpoints %>%
    filter(start_end == 'start') %>%
    pull(nodeID)

  target_nodes <- endpoints %>%
    filter(start_end == 'end') %>%
    pull(nodeID)

  edges <- edges %>%
    mutate(from = source_nodes, to = target_nodes)

  nodes <- endpoints %>%
    distinct(nodeID, .keep_all = TRUE) %>%
    select(-c(edgeID, start_end, xy, X_snap, Y_snap)) %>%
    st_as_sf(coords = c('X', 'Y')) %>%
    st_set_crs(st_crs(edges))

  tbl_graph(nodes = nodes,
            edges = as_tibble(edges),
            directed = directed)
}

#' Build a topologically-clean road network using sfnetworks.
#' Snaps coordinates to a grid (snap_tolerance in CRS units, e.g. meters for
#' EPSG:4087), subdivides edges at crossings, and removes pseudo-nodes.
#' Returns a tbl_graph compatible with tidygraph/igraph operations.
build_road_network <- function(x, directed = FALSE, snap_tolerance = 1) {
  x <- st_cast(x, "LINESTRING")

  # precision = 1/snap_tolerance: round(coord * precision) / precision
  # e.g. snap_tolerance = 1 m  -> precision = 1 -> round to nearest meter
  x <- st_set_precision(x, 1 / snap_tolerance)
  x <- st_make_valid(x)

  net <- sfnetworks::as_sfnetwork(x, directed = directed)

  net <- net %>%
    tidygraph::convert(sfnetworks::to_spatial_subdivision) %>%
    tidygraph::convert(sfnetworks::to_spatial_smooth)

  net <- net %>%
    activate(edges) %>%
    mutate(Shape_Leng = as.numeric(sfnetworks::edge_length()))

  as_tbl_graph(net)
}