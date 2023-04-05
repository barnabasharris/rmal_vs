
library(trajr)
library(terra)
library(parallel)
library(tidyr)
require(purrr)
require(glue)
require(geosphere)
require(data.table)
require(future.apply)
require(future)
require(furrr)
require(curl)
require(purrr)
require(sf)
require(dplyr)
require(stringr)
require(whitebox)
require(stars)
library(rvest)
library(raster)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# sfPolygon <- landscapeMax
extractSRTM <- function(srtmLoc="/media/mal/srtm",
                        sfPolygon, 
                        aggregateBy=1) {
  
  # get relevant SRTM data
  srtmIndex <- st_read(file.path(srtmLoc,'index.shp'))
  
  landscapeMax.union <- 
    sfPolygon %>% 
    st_union() %>% 
    st_sf %>% 
    st_transform(4326) %>% 
    st_filter(srtmIndex, .)
  
  if (nrow(landscapeMax.union)==0) return(NULL)
  
  srtm <- map(landscapeMax.union$location,rast)
  names(srtm) <- NULL
  if (length(srtm) > 1) {
    mos <- do.call(terra::merge, srtm)
  } else mos <- srtm[[1]]
  
  # get elevation, first maximum possible elevation needed
  landscapeMaxSp.wgs <- vect( sfPolygon %>% 
                                st_union() %>% 
                                st_sf %>% 
                                st_transform(4326) )
  
  mos.m <- trim(mask(mos, landscapeMaxSp.wgs))
  
  mos.m.gb <- project(mos.m, 'EPSG:27700')
  
  if (aggregateBy > 1) {
    landscapeMaxElevBuff.t <- 
      terra::aggregate(mos.m.gb,
                       fact = aggregateBy,
                       fun='mean')
  } else return(mos.m.gb)
  return(landscapeMaxElevBuff.t)
}


 sliceTerrain <- 
  function(secDat = sectionData,
           landscapeMetrics = 
             c('slope', 'TPI', 'TRI', 'roughness'),
           metricFocalDistance=500,
           viewdist=10000,
           rdsDataDir = '/media/mal/working_files/rmal/rdata/rds',
           srtmDir = '/media/mal/srtm',
           wbtLoc = '/home/barneyharris/git/whitebox-tools/target/release/whitebox_tools',
           processedSrtmDir = '/media/mal/processed_srtm',
           cumVizLoc = '/media/mal/vis_raster/britain_cum_vis.tif',
           tmpDir = '/tmp',
           clusterObj=NULL,
           refresh=F) {
  
    uid <- secDat$uid
    cl <- clusterObj
    ##┣ viewshed function ----
    # x <- ewPetLocs[[1]]
    wbt_cl_viewshed <- function(x, pols.r.Loc) {
      library(terra)
      library(glue)
      library(whitebox)
      library(purrr)
      wbt_init(exe_path = wbtLoc)
      
      print(paste0('processing',x,'...'))
      
      map_df_progress <- function(.x, .f, ..., .id = NULL) {
        .f <- purrr::as_mapper(.f, ...)
        pb <- progress::progress_bar$new(total = length(.x), format = " [:bar] :current/:total (:percent) eta: :eta", force = TRUE)
        
        f <- function(...) {
          pb$tick()
          .f(...)
        }
        purrr::map_df(.x, f, ..., .id = .id)
      }
      
      vpCoords <- vect(x)
      # row_num <- 1
      vpVsStats <- 1:nrow(vpCoords) %>%
        # vpVsStats <- 1:200 %>%
        map_df_progress(.f = function(row_num) {
          vpRaw <- vpCoords[row_num,]
          vpRawLoc <- glue('{tmpDir}/uid{uid}_cell_number_{vpRaw$cellnum}_vpRaw.shp')
          pat <- paste0(substr(basename(vpRawLoc),1,nchar(basename(vpRawLoc))-4),'*')
          if (file.exists(vpRawLoc)) {
            file.remove(
              list.files(tmpDir,
                         pattern=pat,
                         full.names = T)
            )
          }
          writeVector(vpRaw,vpRawLoc,
                      overwrite = T)
          
          visIndexOutput <- glue('{tmpDir}/viz_out_uid{uid}_cell_number_{vpRaw$cellnum}.tif')
          
          wbt_viewshed(dem = rasVpDemLoc, stations = vpRawLoc, height=1.6,
                       output = visIndexOutput)
          
          pols.r <- rast(pols.r.Loc)
          r <- rast(visIndexOutput)
          r.stats <- zonal(r, pols.r,
                           fun = sum,
                           na.rm=T)
          names(r.stats)[2] <- 'sum'
          r.stats$cell_number <- vpRaw$cellnum
          r.stats$e_uid <- uid
          return(r.stats)
        })
      return(vpVsStats)
    }
    
    ##┣ raster process ----
    # density of 0.1 gives one point per 10 metres (matches SRTM res)
    warning('processing buffers and SRTM rasters...')
    slicePoints <- secDat$sliceLines %>% 
      st_transform(27700) %>% 
      st_line_sample(density=0.1) %>% 
      st_sf %>% 
      mutate(line_id = 1:nrow(.)) %>% 
      st_cast('POINT') %>% 
      group_by(line_id) %>% 
      mutate(p_id = 1:n()) %>% 
      ungroup()
    
    sliceLength <- secDat$sliceLines[1,] %>% 
      st_transform(27700) %>% 
      st_length %>% 
      round %>% 
      as.numeric
    
    # buff required for landscape terrain metrics
    buffDistance <- 
      metricFocalDistance + 
      sliceLength / 2
    
    landscapeMax <- secDat$sliceLines %>% 
      st_transform(27700) %>% 
      st_buffer(buffDistance) %>% 
      st_union()
    
    # sfPolygon <- landscapeMax
    
    landscapeMaxElevBuff.t <- 
      extractSRTM(sfPolygon = landscapeMax,
                  aggregateBy = 1,
                  srtmLoc = srtmDir)
    
    if (is.null(landscapeMaxElevBuff.t)) {
      print('no SRTM coverage!')
      return(NULL)
    }
    
    # buff required for landscape terrain metrics
    vpBuffDistance <- 
      metricFocalDistance + 
      viewdist
    
    landscapeMaxVp <- secDat$sliceLines %>% 
      st_transform(27700) %>% 
      st_buffer(vpBuffDistance) %>% 
      st_union()
    
    # sfPolygon <- landscapeMax
    landscapeMaxElevBuffVp.t <- 
      extractSRTM(sfPolygon = landscapeMaxVp,
                  aggregateBy = 1,
                  srtmLoc = srtmDir)
    
    
    ##┣ create zonal polygons ----
    warning('creeating zonal polygons and rasters...')
    # function for generating buffer zone split by ew line
    # select nearest vertices on bearing of either end of line
    {
      a <- secDat$ew_points[1,] %>% 
        st_transform(4326)
      b <- secDat$ew_points[nrow(secDat$ew_points)-1,] %>% 
        st_transform(4326)
      
      # construct line projecting from end of ew line
      a.dest <- 
        geosphere::destPoint(st_coordinates(a),
                             a$bearings + 180,
                             vpBuffDistance + 5000) %>%
        as.data.frame() %>% 
        st_as_sf(coords=c('lon','lat'),crs=4326) %>% 
        st_transform(27700)
      
      b.dest <- 
        geosphere::destPoint(st_coordinates(b),
                             b$bearings,
                             vpBuffDistance + 5000) %>%
        
        as.data.frame() %>% 
        st_as_sf(coords=c('lon','lat'),crs=4326) %>% 
        st_transform(27700)
      
      
      # plot(landscapeMaxElevBuff.t)
      # plot(secDat$joinedLinear$geometry,add=T)
      # plot(a.dest,add=T)
      # plot(b.dest,add=T)
      
      # construct split line using original line + projected points
      ewSplitLine <- rbind(
        st_coordinates(a.dest),
        st_coordinates(secDat$ew_points),
        st_coordinates(b.dest)
      ) %>% 
        as.data.frame() %>% 
        st_as_sf(coords=c('X','Y')) %>% 
        st_sf %>% 
        summarise(geometry = st_combine(geometry)) %>% 
        st_cast('LINESTRING') %>% 
        st_set_crs(st_crs(landscapeMax)) %>% 
        st_make_valid() %>% 
        st_as_sfc
      
      # plot(landscapeMaxVp)
      # plot(ewSplitLine,add=T)
      # split bufferp polygon with line
      pols <- 
        lwgeom::st_split(landscapeMaxVp,
                         ewSplitLine) %>% 
        st_collection_extract(c("POLYGON")) %>% 
        st_sf %>% 
        mutate(area = as.numeric(st_area(geometry))) %>% 
        slice_max(order_by = area, n=2) %>% 
        mutate(id = 1:2)
      
      pols.r <- pols %>% 
        st_sf %>% 
        vect() %>% 
        rasterize(x=.,
                  y=landscapeMaxElevBuffVp.t,
                  field="id")
      
      # wrap zonal poly for parallel use
      pols.r.wrapped <- wrap(pols.r)
      # write pols.r for same use
      pols.r.Loc <- glue('{tmpDir}/ew{uid}_zonal_pols.tif',
                         overwrite=T)
      writeRaster(pols.r,pols.r.Loc)
    }
    # plot(pols[1,]$geometry,add=T,
    #      col='red')
    # plot(pols[2,]$geometry,add=T,
    #      col='blue')
    
    ##┣ path analysis -----
    warning('analysing path and surrounding terrain...')
    trj <- secDat$ew_points %>% 
      st_sf %>% 
      st_set_geometry('geometry') %>% 
      st_coordinates(.) %>% 
      as.data.frame() %>% 
      TrajFromCoords(., spatialUnits = 'metres')
    
    traj_angles <- 
      c(NA,rad2deg(TrajAngles(trj)),NA)
    
    ew.trj <- 
      data.frame(
        e_uid = secDat$uid,
        trjs_length = TrajLength(trj),
        trjs_straightness = TrajStraightness(trj),
        trjs_sinu = TrajSinuosity2(trj),
        trjs_mean_angle_change = mean(abs(traj_angles),na.rm=T),
        trjs_sd_angle_change = sd(abs(traj_angles),na.rm=T)
      )
    
    ##┣ landscape metrics ----
    terrainMetrics <-
      terrain(x = landscapeMaxElevBuff.t,
              v=landscapeMetrics,
              neighbors = 8)
    
    neighCells <- 
      ceiling(
        metricFocalDistance / 
          mean(terra::res(landscapeMaxElevBuff.t))
      )
    
    # check if odd
    if (neighCells %% 2 == 0) neighCells <- neighCells + 1
    
    centreCell <- ceiling(neighCells^2/2)
    focalTPI <- function(x) {
      if (all(!is.na(x))) {
        x[centreCell] - 
          mean(x[-centreCell])
      } else NA
    }
    warning('calculating custom focal TPI...')
    TPI <- focal(landscapeMaxElevBuff.t, 
                 w=neighCells,
                 fun=focalTPI)
    
    terrainMetrics$TPI_customN <- 
      TPI
    
    # tif based metrics
    # write smaller dem
    rasMetricDemLoc <- glue('{tmpDir}/metric_dem_euid_{secDat$uid}.tif')
    writeRaster(landscapeMaxElevBuff.t,rasMetricDemLoc,overwrite=T)
    # write viewshed dem
    rasVpDemLoc <- glue('{tmpDir}/vp_dem_euid_{secDat$uid}.tif')
    writeRaster(landscapeMaxElevBuffVp.t,rasVpDemLoc,overwrite=T)
    
    # set up WBT
    wbt_init(exe_path = wbtLoc)
    
    # geomorphons
    geoSearch <- 50
    geoOut <- glue('{processedSrtmDir}/geomorphons/ew{secDat$uid}_geomorphons_search{geoSearch}.tif')
    wbt_geomorphons(dem = rasMetricDemLoc,
                    output = geoOut,
                    search = geoSearch)
    # import geomorph maps and extract along ew
    geoMorph <- rast(geoOut)
    ew_points_vect <- vect(secDat$ew_points)
    
    # define summary stats
    # analyse start, end, middle, and entirity of earthwork in terms of geomorph
    startIndices <- 1:ceiling(nrow(ew_points_vect)*0.1)
    endIndices <- nrow(ew_points_vect):(nrow(ew_points_vect)-ceiling(nrow(ew_points_vect)*0.1))
    middleIndices <- (ceiling(nrow(ew_points_vect)*0.45)):(ceiling(nrow(ew_points_vect)*0.55))
    allIndices <- list(end_a = startIndices,
                       end_b = endIndices,
                       middle = middleIndices,
                       all = 1:nrow(ew_points_vect))
    geoIndices <- allIndices %>% 
      map(.f = function(x) {
        extract(geoMorph,ew_points_vect[x,])
      }) |> rbindlist(idcol = 'length_group') %>% 
      as.data.frame() %>% 
      rename('geomorph' = 3)
    
    geoIndices.sum <- 
      geoIndices %>% 
      group_by(length_group) %>% 
      summarise(
        count_1 = length(geomorph[geomorph==1]),
        count_2 = length(geomorph[geomorph==2]),
        count_3 = length(geomorph[geomorph==3]),
        count_4 = length(geomorph[geomorph==4]),
        count_5 = length(geomorph[geomorph==5]),
        count_6 = length(geomorph[geomorph==6]),
        count_7 = length(geomorph[geomorph==7]),
        count_8 = length(geomorph[geomorph==8]),
        count_9 = length(geomorph[geomorph==9]),
        count_10 = length(geomorph[geomorph==10])
      ) %>% 
      ungroup() %>% 
      mutate(modal_geomorph = names(.[,2:11])[max.col(.[,2:11])]) %>%
      rowwise(length_group) %>% 
      mutate(total_count = sum(c_across(where(is.numeric))),
             total_unique = sum(c_across(starts_with('count_')) > 0, na.rm = TRUE)
      ) %>% 
      mutate(across(count_1:count_10, ~ round((.x / total_count)*100,2),
                    .names = "per_{.col}"))
    
    # elevation
    terrainMetrics$elev <- landscapeMaxElevBuff.t
    
    # cum visibility
    vp <- terra::rast(cumVizLoc)
    
    vp.resamp <- resample(vp,landscapeMaxElevBuff.t)
    vp.masked <- mask(vp.resamp,landscapeMaxElevBuff.t)
    
    terrainMetrics$cumviz <- vp.masked
    
    ##┣ viewshed calculations ----
  
    workers = length(clusterObj)
    ###•┣ viewshed data preparation  ----
    warning('preparing data for viewshed analyses...')
    # viewsheds, first do along path of ew
    # multiple ew locations per ras cell, find unique cells to avoid
    # needless repeating of viewshed call
    {
      ewLocsVect <- vect(secDat$ew_points)
      ewCoords <- terra::crds(ewLocsVect)
      cellsNums <- cellFromXY(landscapeMaxElevBuffVp.t,ewCoords)
      ewLocsVect$cell_number <- cellsNums
      
      # make unique =
      cellsNums.unique <- unique(cellsNums)
      uniqueCrds <- xyFromCell(landscapeMaxElevBuffVp.t,cellsNums.unique)
      vpCoords <- vect(uniqueCrds, crs='EPSG:27700')
      vpCoords$cellnum <- cellsNums.unique
      
      # split points into 4 for parallel processing
      q <- ceiling(nrow(vpCoords) / workers)
      vecChunks <- split(1:nrow(vpCoords), ceiling(seq_along(1:nrow(vpCoords))/q))
      ewLocs <- names(vecChunks) %>% 
        map(.f = function(x) {
          vpLoc <- glue('{tmpDir}/uid{uid}_ew_locs_part{x}.shp')
          writeVector(vpCoords[vecChunks[[x]],], vpLoc,
                      overwrite=T)
          vpLoc
        })
    }
    
    # now peturbated locations vs
    # multiple pet locations per ras cell, find unique cells to avoid
    # needless repeating of viewsehd script
    {
      ewPetCoords <- st_coordinates(slicePoints)
      
      cellsNums <- cellFromXY(landscapeMaxElevBuffVp.t,ewPetCoords)
      slicePoints$cell_number <- cellsNums
      
      cellsNums.unique <- unique(cellsNums)
      uniqueCrds <- xyFromCell(landscapeMaxElevBuffVp.t,cellsNums.unique)
      vpCoords <- vect(uniqueCrds, crs='EPSG:27700')
      vpCoords$cellnum <- cellsNums.unique
      
      # split points into 4 for parralel processing
      q <- ceiling(nrow(vpCoords) / workers)
      vecChunks <- split(1:nrow(vpCoords), ceiling(seq_along(1:nrow(vpCoords))/q))
      ewPetLocs <- names(vecChunks) %>% 
        map(.f = function(x) {
          vpLoc <- glue('{tmpDir}/uid{uid}_ew_pet_locs_part{x}.shp')
          writeVector(vpCoords[vecChunks[[x]],], vpLoc,
                      overwrite=T)
          vpLoc
        })
    }
    
    ###•┣ summed vs map  ----
    
    if (!file.exists(glue('{processedSrtmDir}/vs/ew{uid}_summed_vs_rad{viewdist}.tif'))) {
      warning('generating summed binary viewshed map...')  
      # first produce summed binary raster map showing viewshed along whole
      # earthwork
      
      clusterExport(cl, c('rasVpDemLoc','uid', 'wbtLoc'),envir=environment())
      st <- Sys.time()
      clusterApply(cl, ewLocs, function(x) {
        library(whitebox)
        library(stringr)
        wbt_init(exe_path = wbtLoc)
        
        summedVsOutput <- str_replace(x,'.shp','.tif')
        
        wbt_viewshed(dem = rasVpDemLoc, stations = x, height=1.6,
                     output = summedVsOutput)
        
      })
      Sys.time() - st

      # merge resulting parts
      summedVsOutputs <- ewLocs %>% 
        map(str_replace,'.shp','.tif')
      # load and sum
      r <- rast(unlist(summedVsOutputs))
      r.sum <- sum(r)
      r.sum[r.sum==0] <- NaN
      r.sum.mean <- focal(r.sum,fun=mean,na.rm=F)
      
      writeRaster(r.sum.mean,glue('{processedSrtmDir}/vs/ew{uid}_summed_vs_rad{viewdist}.tif'),
                  overwrite=T)
    } else {
      r.sum.mean <- rast(glue('{processedSrtmDir}/vs/ew{uid}_summed_vs_rad{viewdist}.tif'))
    }
    
    summedVsStats <- 
      terra::zonal(r.sum.mean,pols.r,
                   fun = sum,
                   na.rm=T)
    
    ###•┣ ew + peturbated location vs stats  ----
    if (
      !file.exists(glue('{rdsDataDir}/ew{uid}_vs_ew_pet_stats.RDS')) |
      refresh==T
    ) {
      warning('measuring viewshed along peturbated earthwork points...')
      clusterExport(cl, c('rasVpDemLoc','uid','pols.r.wrapped'),
                    envir = environment())
      st <- Sys.time()
      vsZonePetStats <- 
        clusterApply(cl, ewPetLocs, wbt_cl_viewshed, pols.r.Loc)
      print(Sys.time() - st)
      # run parallel !
      ewPetLocVsStats <- bind_rows(vsZonePetStats) %>% 
        tidyr::pivot_wider(id_cols = c('cell_number','e_uid'),
                           names_from = id,
                           values_from = sum,
                           names_prefix = 'zone_')
      
      # join to original line ids etc.
      ewPetLocsVect.j <- slicePoints %>% 
        left_join(ewPetLocVsStats,
                  by='cell_number') %>% 
        vect()

      saveRDS(ewPetLocsVect.j,
              file=glue('{rdsDataDir}/ew{uid}_vs_ew_pet_stats.RDS'))
    } else {
      ewPetLocsVect.j <- readRDS(glue('{rdsDataDir}/ew{uid}_vs_ew_pet_stats.RDS'))
    }
    
    warning('finished!')
    # > Time difference of 3.72654 mins (workers = 6)
    # > Time difference of 5.952829 mins (uid = 208, workers = 4)
    
    
    ewPetLocsVect.z1 <- rasterize(ewPetLocsVect.j, terrainMetrics, field='zone_1')
    ewPetLocsVect.z2 <- rasterize(ewPetLocsVect.j, terrainMetrics, field='zone_2')
    
    terrainMetrics$vs_zone1 <- ewPetLocsVect.z1
    terrainMetrics$vs_zone2 <- ewPetLocsVect.z2
    
    ##┣ extract and aggregate data -----
    terrainData <- st_as_stars(terrainMetrics) %>% 
      st_set_crs(st_crs(slicePoints))
    
    slicePoints.wgs <- slicePoints %>% 
      st_transform(st_crs(terrainData)) %>% 
      st_sf
    
    print('extracting terrain data at slice points...')
    terPet <- st_extract(terrainData,slicePoints) %>% 
      st_as_sf %>% st_sf %>% 
      st_join(slicePoints, st_equals) %>% 
      mutate(ew_class = 
               ifelse(p_id == 25, 'ew','slice'),
             wt = abs(p_id - 25) + 1)
    
    allTerrain <- terPet %>% 
      st_drop_geometry() %>% 
      tidyr::pivot_longer(cols = !!names(terrainMetrics))
    
    terrainDat <- 
      list(terPet = terPet,
           allTerrain = allTerrain,
           zonePols = pols,
           summedVsStats = summedVsStats,
           geoMorphStats = geoIndices.sum,
           trajAngles = traj_angles)
    return(terrainDat)
  }

rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

