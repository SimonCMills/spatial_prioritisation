// script to generate a series of rasters at a 1*1km resolution as  input 
// layers to the spatial prioritisation exercise. 
// Information we need to generate:
// (1) average elevation
// (2) percentage "forest" treecover (area covered by pixels with >50% 
// canopy cover that is contiguous with many other like pixels)
// (3) percentage "not forest" treecover (i.e. treecover that is not 
// contiguous with large areas of forest)
// 
// ultimately want a series of rasters that all: 
// (1) span the eastern cordillera of Colombia
// (2) share the same resolution and projection (1km res)
// (3) are clipped to be above 1000 m.a.s.l. (on reduced res ele raster)


// var urban_area = ee.Image("DLR/WSF/WSF2015/v1"),
//     imageCollection = ee.ImageCollection("CSP/HM/GlobalHumanModification"),
//     imageCollection2 = ee.ImageCollection("JRC/GHSL/P2016/SMOD_POP_GLOBE_V1"),
//     hansen = ee.Image("UMD/hansen/global_forest_change_2021_v1_9"),
//     ele = ee.Image("CGIAR/SRTM90_V4"),
//     mountains = ee.FeatureCollection("users/scmills/GMBAMountainInventory_v1_2-World"),
//     countries = ee.FeatureCollection("USDOS/LSIB/2013"),
//     inland_water = ee.ImageCollection("GLCF/GLS_WATER"),
//     copernicus = ee.ImageCollection("COPERNICUS/Landcover/100m/Proba-V-C3/Global");
    
// layer: eastern cordillera of Colombia
var EC = mountains.filter("Name == 'Cordillera Oriental Colombia Venezuela'");

// layers: treecover and water
var tc = hansen.select("treecover2000");

// GLCF is covered in cloud.. so use copernicus
var land = copernicus
  .select("discrete_classification")
  .toBands()
  .select("2019_discrete_classification")
  .eq(80).neq(1).clip(EC);

// clip elevation, urban areas, and treecover (to add: other infrastructure) 
// rasters to the EC 
var ele_EC = ele.clip(EC);

// drop resolution of elevation raster to 1*1km 
// var ele_EC_1km = ele_EC
//     // Force the next reprojection to aggregate instead of resampling.
//     .reduceResolution({
//       reducer: ee.Reducer.mean(),
//       maxPixels: 1024
//     })
//     // Request the data at the scale and projection of the MODIS image.
//     .reproject({
//       crs: hansen.projection(), 
//       scale: 1000
//     });
//------------------------------------------------------------------------------ 
// drop resolution of EC elevation layer to 177: this corresponds to the area 
// sampled (100^2*pi). Also reproject to Hansen. 
var ele_EC_177m = ele_EC
    // Force the next reprojection to aggregate instead of resampling.
    .reduceResolution({
      reducer: ee.Reducer.mean(),
    })
    // Request the data at the scale and projection of the MODIS image.
    .reproject({
      crs: hansen.projection(), 
      scale: 177
    });
    
    
// get cells above 1000m in the reduced resolution EC raster
// previously ele_EC2
var ele_processed = ele_EC_177m.updateMask(ele_EC_177m.gte(1000));

// need to get urban area into a [1,0] format 
var ua_EC = urban_area.eq(255).unmask(0).clip(EC);

//// urban area processing ////
// refine the urban area raster so that small "not urban" sections that lie
// entirely within "urban" areas are reclassified as "urban". 

// resample urban areas down to 50 metre resolution and classify any cell with
// greater than 50% urban cover as urban. 
var ua_resampled = ua_EC
  .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 1024
    })
  .reproject(
    {crs:ele_EC.projection(), 
    scale:50})
    .gt(.5);

// a cluster is found if it is separated from at least 99 other like-pixels
// and assigned a unique value
var ua_EC_clust = ua_resampled.connectedComponents({
  connectedness: ee.Kernel.plus(1),
  maxSize: 100
})
// this has identified clusters of (1) urban areas surrounded by "not urban"
// areas and (2) "not urban areas" surrounded by urban areas, 
// the trick now is to extract just the "not urban areas" (i.e. areas 
// that are not urban and lie internal to urban areas)

// I don't *think* any clusters get assigned a value of 0, but to be sure I 
// have set a potential 0-valued cluster to 1. 
var ua_EC_clust_0safe = ua_EC_clust.where(ua_EC_clust.eq(0), 1)

// the first line sets all clusters that are clusters of "urban" cells to 0 
// and then queries for all non-0 clusters (i.e. gets clusters of "not urban"
// cells, which by definition are gaps internal to urban areas). 
// last line just selects the labels band that got created in this process
var ua_EC_internal_gaps = ua_EC_clust_0safe
  .where(ua_resampled.eq(1), 0)
  .neq(0)
  .select("labels")

// lastly, update values in the resampled raster that have been identified
// as internal clusters above. 
var ua_resampled_filled = ua_resampled.where(ua_EC_internal_gaps.eq(1), 1); 

// note: could additionally buffer, but I think this is overkill given 
// the above. 
// var buffer = urban_area.eq(1).focal_max({
// radius: 100,
// units:'meters',
// iterations: 1
// })

//// TO DO: add additional infrastructure layers ////

//// FOREST AREA PROCESSING ---------------------------------------------------
// clip 2000 treecover to Eastern Cordillera
var tc_EC = tc.clip(EC)
// map to forest and mask by the 177 m resolution >1000 m elevation layer. 
var forest = tc_EC.gte(50).updateMask(ele_processed);

/* 
We don't want to classify small isolated fragments of forest as forest
proper. Follow the same method as with urban areas to identify fragments of 
forest that fall below a particular size.

A cluster is found if it is separated from at least 99 other like-pixels and 
assigned a unique value
*/
var forest_EC_clust = forest.connectedComponents({
  connectedness: ee.Kernel.plus(1),
  maxSize: 100
});

/* 
This has identified clusters of (1) forest areas surrounded by "not forest"
areas and (2) "not forest areas" surrounded by forest areas, the trick now is 
to extract just the "forest areas" (i.e. areas that are not urban and lie 
internal to urban areas)
*/

// as above
var forest_EC_clust_0safe = forest_EC_clust.where(forest_EC_clust.eq(0), 1);

// As before, but now getting clusters of forest cells, clusters that are 
// forest, surrounded by non-forest. 
var forest_fragments = forest_EC_clust
  .where(forest.eq(0), 0)
  .neq(0)
  .select("labels")

/* 
Lastly, set any forest cells that are in fragments to 0 (not forest) and 
mask out any cells that are within urban areas or are water (so that they don't
contribute to the average).

var forest_contig is therefore contiguous forest, masked by urban areas and 
waterbodies (note that urban areas have had small "non-urban" areas within them
set to be urban).
*/

var forest_contig = forest
  .where(forest_fragments.eq(1), 0)
  .updateMask(ua_resampled_filled.eq(0))
  .updateMask(land);
  
// reproject to 177 m resolution 
var forest_contig_177m = forest_contig
    // Force the next reprojection to aggregate instead of resampling.
    .reduceResolution({
      reducer: ee.Reducer.mean()
    })
    // Request the data at the scale and projection of the MODIS image.
    .reproject({
      crs: forest_contig.projection(),
      scale:177
    })
    
// store this mask layer too; and get the proportion of each 177 m cell that 
// is unmasked
var unmasked = forest_contig.mask();
  
var prop_unmasked_177m = unmasked
    // Force the next reprojection to aggregate instead of resampling.
    .reduceResolution({
      reducer: ee.Reducer.mean()
    })
    // Request the data at the scale and projection of the MODIS image.
    .reproject({
      crs: forest_contig.projection(),
      scale:177
    });

// mask treecover layer on urban areas and waterbodies (i.e. the forest_contig
// mask, and additionally mask out contiguous forest 
var tc_masked = tc_EC
  .updateMask(forest_contig.mask())
  .updateMask(forest_contig.eq(0));
  
var tc_masked_177m = tc_masked
    // Force the next reprojection to aggregate instead of resampling.
    .reduceResolution({
      reducer: ee.Reducer.mean()
    })
    // Request the data at the scale and projection of the MODIS image.
    .reproject({
      crs: forest_contig.projection(),
      scale:177
    });
    
  
 
//Map.addLayer(tc_EC, {}, "tc"); 
//Map.addLayer(forest.updateMask(forest), {}, "forest");
Map.addLayer(forest_contig.updateMask(forest_contig), {}, "forest_contig");
Map.addLayer(forest_contig_177m, {}, "forest_contig_177");

Map.addLayer(unmasked, {}, "unmasked");
Map.addLayer(prop_unmasked_177m, {}, "unmasked_177");

Map.addLayer(tc_masked, {}, "treecover2000");
Map.addLayer(tc_masked_177m, {}, "treecover2000_177");


//// EXPORT ////
Export.image.toDrive({
    image: forest_contig_177m,
    description: "pct_contig_forest_177m_res",
    folder: "gee_exports",
    scale: 177,
    region: EC
  });

Export.image.toDrive({
    image: tc_masked_177m,
    description: "pct_treecover_in_pasture_177m_res",
    folder: "gee_exports",
    scale: 177,  
    region: EC
  });
  
Export.image.toDrive({
    image: prop_unmasked_177m,
    description: "proportion_unmasked_177m",
    folder: "gee_exports",
    scale: 177,  
    region: EC
  });
  
Export.image.toDrive({
    image: ele_EC_177m,
    description: "elev_ALOS_177m_res",
    folder: "gee_exports",
    scale: 177,
    region: EC
  });
