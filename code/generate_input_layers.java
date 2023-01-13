///
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

// get treecover
var tc = hansen.select("treecover2000");

// select eastern cordillera of Colombia
var EC = mountains.filter("Name == 'Cordillera Oriental Colombia Venezuela'")

// clip elevation, urban areas, and treecover (to add: other infrastructure)
// rasters to the EC
var ele_EC = ele.clip(EC);
// need to get urban area into a [1,0] format
var ua_EC = urban_area.eq(255).unmask(0).clip(EC);
var tc_EC = tc.clip(EC);
// tc to forest layer
var forest = tc_EC.gte(50);

// drop resolution of elevation raster to 1*1km
var ele_EC_1km = ele_EC
    // Force the next reprojection to aggregate instead of resampling.
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 1024
    })
    // Request the data at the scale and projection of the MODIS image.
    .reproject({
      crs: ele_EC.projection(),
      scale: 1000
    });

// get cells above 1000m in the reduced resolution EC raster
var ele_EC2 = ele_EC_1km.updateMask(ele_EC_1km.gte(1000));


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

//// FOREST AREA PROCESSING ////
// We don't want to classify small isolated fragments of forest as forest
// proper. Follow the same method as with urban areas to identify fragments
// of forest that fall below a particular size.

// a cluster is found if it is separated from at least 99 other like-pixels
// and assigned a unique value
var forest_EC_clust = forest.connectedComponents({
  connectedness: ee.Kernel.plus(1),
  maxSize: 100
});
// this has identified clusters of (1) forest areas surrounded by "not forest"
// areas and (2) "not forest areas" surrounded by forest areas,
// the trick now is to extract just the "forest areas" (i.e. areas
// that are not urban and lie internal to urban areas)

// as before
var forest_EC_clust_0safe = forest_EC_clust.where(forest_EC_clust.eq(0), 1);

// As before, but now getting clusters of forest cells, clusters that are
// forest, surrounded by non-forest.
var forest_fragments = forest_EC_clust
  .where(forest.eq(0), 0)
  .neq(0)
  .select("labels")

// lastly, set any forest cells that are in fragments to 0 (not forest) and
// mask out any "forest cells" that are within urban areas (in principle
// shouldn't happen much)
var forest_contig = forest
  .where(forest_fragments.eq(1), 0)
  .updateMask(ua_resampled_filled.eq(0));


// mask treecover layer on urban areas and forest
var tc_masked = tc
  .updateMask(ua_resampled_filled.eq(0))
  .updateMask(forest_contig.eq(0));



//Map.addLayer(buffer, {palette: "red"});
//Map.addLayer(ele_EC2, {min: 0, max: 5000}, "elevation");

//Map.addLayer(ua_EC_filled, {}, "filled");

// Map.addLayer(ua_resampled_filled);
// Map.addLayer(ua_EC_filled);
///Map.addLayer(forest_masked, {}, "forest_masked")
Map.addLayer(forest_contig, {}, "forest_contig")
Map.addLayer(tc_masked, {}, "tc_masked")
Map.addLayer(ua_resampled_filled, {}, "ua")
// Map.addLayer(ua_resampled_filled2, {}, "ua2")



// Map.addLayer(imageCollection2, {min:0, max:3, palette:["red", "blue", "green", "yellow"]})