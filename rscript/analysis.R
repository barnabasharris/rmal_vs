dataDir = '/home/tcrnbgh/secdat_initials'
source('rscript/functions.R')

print('getting MPI cluster...')
cl <- snow::getMPIcluster()
print('done!')

# Display info about each process in the cluster
print(clusterCall(cl, function() Sys.info()))

# load secdata
initials <- list.files(file.path(dataDir,'rds'),
                       pattern = 'sectionData_initial*',
                       full.names = T)
  
x <- initials[100]
  
l <- readRDS(x)

terrainDat <- 
  sliceTerrain(secDat = l,
               srtmDir = '/home/tcrnbgh/srtm',
               wbtLoc = '/home/tcrnbgh/WBT/whitebox_tools',
               rdsDataDir = '/home/tcrnbgh/Scratch/rds',
               processedSrtmDir = '/home/tcrnbgh/Scratch/processed_srtm',
               cumVizLoc = '/home/tcrnbgh/britain_cum_vis.tif',
               tmpDir = '/home/tcrnbgh/Scratch/tmp',
               refresh=T,
               clusterObj = cl)
  
# initials %>% 
#   walk(.f = function(x) {
#     print(x)
#     l <- readRDS(x)
#     # print(l)
#     uid <- str_extract(x,'\\d+')
#     l$terrainDat$terPet$e_uid <- uid
#     l$terrainDat$allTerrain$e_uid <- uid
#     terrainDat <- sliceTerrain(secDat = l,
#                                refresh=T)
#     l$terrainDat <- terrainDat
#     saveRDS(l, file=str_replace(x,'initial','terrain_w_vs'))
#     # return(l)
#   })

stopCluster(cl)
mpi.quit()