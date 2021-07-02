// Required Data Inputs 
// ===================
// * USGS/NASA's Landsat 4 surface reflectance tier 1 dataset (August 1982 - December 1993)
// * USGS/NASA's Landsat 5 surface reflectance tier 1 dataset (January 1, 1984 - May 5, 2012)
// * USGS/NASA's Landsat 7 surface reflectance tier 1 dataset (January 1, 1999 - December 31, 2019)
// * USGS/NASA's Landsat 8 surface reflectance tier 1 dataset (April 11, 2013 - December 31, 2019)
// * Study Area Polygon

var countries = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017");
var country = 'Egypt';

var ALGO = "MEAN"; // PERCENTILE75  // PERCENTILE65 // PERCENTILE60 // MEDIAN
var cloudCoveragePercentage = 80;

var studyArea = countries.filter(ee.Filter.eq('country_na',country ))

if (country === 'Egypt'){
  studyArea = ee.FeatureCollection("users/derickongeri/EgyptGrids")//.filter(ee.Filter.eq('id',3))
  var tileList = studyArea.aggregate_array('id').distinct();
  tileList = tileList.getInfo();
  print(tileList)
}

Map.centerObject(studyArea);



var trend_baselinePeriod = ee.List.sequence(2000, 2015, 1).getInfo(); //baseline years for reporting trend 2000-2015 (16 years)
var trend_reportingPeriod = ee.List.sequence(2005, 2020, 1).getInfo(); //reporting period for trend 2005 - 2020 (16 years)

var state_baselinePeriod = ee.List.sequence(2000, 2012, 1).getInfo(); // first 13year for reporting state 2000 - 2012
var state_baselinePeriod_final = ee.List.sequence(2013, 2015, 1).getInfo(); // final 3 yars of the reporing period 2013 - 2015

var state_reportingPeriod = ee.List.sequence(2005, 2017, 1).getInfo(); // first 13year of the reporting period
var state_reportingPeriod_final = ee.List.sequence(2018, 2020, 1).getInfo(); // last 3 years of the roporting period

function yearlyNDVI(year){
    var start_date = year+ '-01-01';
    var end_date   = year+ '-12-31';
  
    //--------------------------------------------------------------------
    //       Landsat 4, 5, 7 cloudmask
    //--------------------------------------------------------------------
    
        // If the cloud bit (5) is set and the cloud confidence (7) is high
        // or the cloud shadow bit is set (3), then it's a bad pixel.
      var cloudMaskL7 = function(image) {
      var qa = image.select('pixel_qa');
      var cloud = qa.bitwiseAnd(1 << 5)
                      .and(qa.bitwiseAnd(1 << 7))
                      .or(qa.bitwiseAnd(1 << 3));
      
        // Remove edge pixels that don't occur in all bands
      //var mask2 = image.mask().reduce(ee.Reducer.min())//.focal_min(300,'square','meters').eq(0);
      //var mask2 = image.select('B4').reduce(ee.Reducer.min()).gt(0)//.focal_min(500,'square','meters');
      // Remove edge pixels that don't occur in all bands
      var mask3 =  
                  (image.select('B3').gt(100))
                  .and(image.select('B4').gt(100))
                  
    
                  .and(image.select('B4').lt(10000))
                  .and(image.select('B3').lt(10000))

      return image.updateMask(cloud.not()).updateMask(mask3)//.updateMask(mask2)//.clip(image.geometry().buffer(-5000))//.or(mask3));
    };

      var cloudMaskL45 = function(image) {
      var qa = image.select('pixel_qa');
      var cloud = qa.bitwiseAnd(1 << 5)
                      .and(qa.bitwiseAnd(1 << 7))
                      .or(qa.bitwiseAnd(1 << 3));
      
      // Remove edge pixels that don't occur in all bands
      //var mask2 = image.mask().reduce(ee.Reducer.min());
        var mask2 =  
                  (image.select('B3').gt(100))
                  .and(image.select('B4').gt(100))
                  
    
                  .and(image.select('B4').lt(10000))
                  .and(image.select('B3').lt(10000))
      
      return (image.updateMask(cloud.not()).updateMask(mask2))//.clip(image.geometry().buffer(-5000))//.updateMask(mask2);
    };
    
    //--------------------------------------------------------------------
    //         Landsat 8 cloudmask
    //--------------------------------------------------------------------
    
        // Bits 3 and 5 are cloud shadow and cloud, respectively.
    function maskL8sr(image) {
      var cloudShadowBitMask = (1 << 3);
      var cloudsBitMask = (1 << 5);
      
        // Get the pixel QA band.
      var qa = image.select('pixel_qa');
    
        // Both flags should be set to zero, indicating clear conditions.
      var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                     .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
      var mask2 =  
                  
                  (image.select('B5').gt(100))
                  .and(image.select('B4').gt(100))
                  
    
                  .and(image.select('B5').lt(10000))
                  .and(image.select('B4').lt(10000))
                  
       //var mask2 = image.mask().reduce(ee.Reducer.min()).focal_min(500,'square','meters');
      //return image
      return image.updateMask(mask).updateMask(mask2)//.clip(image.geometry().buffer(-5000));
    }
    
        // Apply Cloudmask to L4.5.7
    var L4 = ee.ImageCollection("LANDSAT/LT04/C01/T1_SR")
                      .filterDate(start_date, end_date)
                      .filter(ee.Filter.lessThan('CLOUD_COVER_LAND', cloudCoveragePercentage))
                      .filterBounds(studyArea)
                      .map(cloudMaskL45)
                      .select(['B3', 'B4'], ['RED', 'NIR']);
    var L5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
                      .filterDate(start_date, end_date)
                      .filter(ee.Filter.lessThan('CLOUD_COVER_LAND', cloudCoveragePercentage))
                      .filterBounds(studyArea)
                      .map(cloudMaskL45)
                      .select(['B3', 'B4'], ['RED', 'NIR']);
    var L7a = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
                      .filterDate('1999-01-01', '2003-04-01')
                      .filterDate(start_date, end_date)
                      .filter(ee.Filter.lessThan('CLOUD_COVER_LAND', 100))
                      .filterBounds(studyArea)
                      .map(cloudMaskL7)
                      .select(['B3', 'B4'], ['RED', 'NIR']);
    var L7b = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
                      .filterDate('2012-01-01', '2013-12-31')
                      .filterDate(start_date, end_date)
                      .filter(ee.Filter.lessThan('CLOUD_COVER_LAND', 100))
                      .filterBounds(studyArea)
                      .map(cloudMaskL7)
                      .select(['B3', 'B4'], ['RED', 'NIR']);
                      
    var L7 = L7a.merge(L7b);
    
    var L8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
                      .filterDate(start_date, end_date)
                      .filter(ee.Filter.lessThan('CLOUD_COVER', cloudCoveragePercentage))
                      .filterBounds(studyArea)
                      //.filterBounds(AOI)
                      .map(maskL8sr)
                      .select(['B4', 'B5'], ['RED', 'NIR']);
    
    //--------------------------------------------------------------------
    // Merge Landsat 4, 5, 8 imagery collections and filter all by date/place
    //--------------------------------------------------------------------
    
    var L4578 = L4.merge(L5).merge(L7).merge(L8);
    
    //--------------------------------------------------------------------
    //                  Create NDVI Collection 
    //--------------------------------------------------------------------
    
    var NDVI = function(image) {
      return image.normalizedDifference(['NIR', 'RED']).rename('NDVI'+year);
      //return image.addBands(ndvi);
    };
    
    if (ALGO=='MEDIAN'){
        var suffix = 'median'; 
        var annualNDVI = L4578.map(NDVI).median().clip(studyArea).rename(suffix+'NDVI');
    }
    
    if (ALGO=='MEAN'){
        var suffix = 'mean'; 
        var annualNDVI = L4578.map(NDVI).mean().clip(studyArea).rename(suffix+'NDVI');
    }
    
    if (ALGO=='PERCENTILE75'){
        var suffix = '75pc'; 
        var annualNDVI = L4578.map(NDVI).reduce(ee.Reducer.percentile([75])).clip(studyArea).rename(suffix+'NDVI'+'_'+year);
    }
    
    if (ALGO=='PERCENTILE65'){
        var suffix = '65pc'; 
        var annualNDVI = L4578.map(NDVI).reduce(ee.Reducer.percentile([65])).clip(studyArea).rename(suffix+'NDVI'+'_'+year);
    }
    
    //mask out permanent water bodies
    var surfaceWater_dataset = ee.ImageCollection('JRC/GSW1_3/YearlyHistory').filterMetadata('year', 'equals', year);
    var surfaceWaterYearly = ee.Image(surfaceWater_dataset.first());
    var swater_mask = surfaceWaterYearly.updateMask(surfaceWaterYearly.eq(3));
    var annualNDVI = annualNDVI.where(swater_mask,0)
    
  return annualNDVI
}

//pass computed yearlyNDVI to collection
function ndviCollection(listYears){
  var ImgColl_NDVI = ee.ImageCollection.fromImages(listYears.map(yearlyNDVI))
  return ImgColl_NDVI
}

//coputing NDVI trends for reporting and baseline periods.
var trendBaseline_NDVIcollection = ndviCollection(trend_baselinePeriod);
var trendReporting_NDVIcollection = ndviCollection(trend_reportingPeriod);

//Computing NDVI state for baseline and reporting periods
var stateBline_NDVIcollection = ndviCollection(state_baselinePeriod);
var stateBline_NDVIcollection_final = ndviCollection(state_baselinePeriod_final);

var stateReporting_NDVIcollection = ndviCollection(state_baselinePeriod);
var stateReporting_NDVIcollection_final = ndviCollection(state_baselinePeriod_final);

//=============================================================================================================
// Function to compute state 
function productivityState(coll1, coll2){
  //compute NDVI mean and standard deviation for baseline period
  var coll_1_mean = coll1.mean().rename('mean');
  var coll_1_stdDev = coll1.reduce(ee.Reducer.stdDev()).rename('stdDev');
  var coll_2_mean = coll2.mean().rename('mean');
 //compute z-statistics 
  var zstats = coll_2_mean.expression(
    '( X - U )/(SD/3**0.5)',{
    'X':coll_2_mean.select('mean'),
    'U':coll_1_mean.select('mean'),
    'SD':coll_1_stdDev.select('stdDev'),
  }).rename('zstats')
  
  return zstats
}

//===============================================================================================================

var palette = ['red', 'white', 'green'];
var productivityState_baseline = productivityState(stateBline_NDVIcollection, stateBline_NDVIcollection_final);
var productivityState_reporting = productivityState(stateReporting_NDVIcollection, stateReporting_NDVIcollection_final)

//var trendBaselinePeriod = trendBaseline_NDVIcollection.reduce(ee.Reducer.kendallsCorrelation());
//var trendReportingPeriod = trendReporting_NDVIcollection.reduce(ee.Reducer.kendallsCorrelation());
//print(trendReportingPeriod)

// ============================================================================================================================
// adding the layers to map

// Map.addLayer(productivityState_baseline, {palette:palette}, 'Baseline State');
// Map.addLayer(productivityState_reporting, {palette:palette}, 'reporting State');
// Map.addLayer( trendBaselinePeriod.select('meanNDVI_p-value'), {palette:palette}, 'Baseline p-value')
// Map.addLayer( trendReportingPeriod.select('meanNDVI_p-value'), {palette:palette}, 'Reporting p-value')
// Map.addLayer( trendBaselinePeriod.select('meanNDVI_tau'), {palette:palette}, 'Baseline tau')
// Map.addLayer( trendReportingPeriod.select('meanNDVI_tau'), {palette:palette}, 'Reporting tau')

// =========================================================================================================================
/*
// Seperate result into 5 classes
var thresholds = ee.Image([-1.96, -1.28, 1.28, 1.96, 3]);
var classified = zstats.lt(thresholds).reduce('sum').toInt();

//Map.addLayer(classified, {}, 'state classified');

var ndvi_visualization = {
  min: -0.22789797020331423, 
  max: 0.6575894075894075,
  palette: 'FFFFFF, CE7E45, DF923D, F1B555, FCD163, 99B718, 74A901, 66A000, 529400,' +
    '3E8601, 207401, 056201, 004C00, 023B01, 012E01, 011D01, 011301'
};*/
//======================================MannKendall test=========================================
function mannKendall(imgCollection){
  var afterFilter = ee.Filter.lessThan({
  leftField: 'system:index',
  rightField: 'system:index'
  });
  
  var joined = ee.ImageCollection(ee.Join.saveAll('after').apply({
  primary: imgCollection,
  secondary: imgCollection,
  condition: afterFilter
  }));
  
  var sign = function(i, j) { // i and j are images
  return ee.Image(j).neq(i) // Zero case
      .multiply(ee.Image(j).subtract(i));
  };
  
  var kendall = ee.ImageCollection(joined.map(function(current) {
  var afterCollection = ee.ImageCollection.fromImages(current.get('after'));
  return afterCollection.map(function(image) {
    // The unmask is to prevent accumulation of masked pixels that
    // result from the undefined case of when either current or image
    // is masked.  It won't affect the sum, since it's unmasked to zero.
  return ee.Image(sign(current, image))
  });
  // Set parallelScale to avoid User memory limit exceeded.
  }).flatten()).reduce('sum', 2);
  
  var threshold_z_score = [5, 1.96, 1.28, -1.28, -1.96]
  var classiffied_Kendall = kendall.lt(threshold_z_score).reduce('sum').toInt()
  
  return classiffied_Kendall
}

var trendBaselinePeriod = mannKendall(trendBaseline_NDVIcollection);

//var reportingKendall = mannKendall(collReporting);
Map.addLayer(trendBaselinePeriod, {min:0, max:5, palette:['1b8607','b8ff68','ffffff','ffc443','ff1919']}, 'Baseline p-value');
/*

//==========================================================================================

function kendall_test(imageCollection){
  var TimeSeriesList = imageCollection.toList(13)
  var NumberOfItems = TimeSeriesList.length().getInfo()
  var  ConcordantArray = ee.List([])
  var  DiscordantArray = ee.List([])
    for (var k = 0; k < 14; k++){
      var CurrentImage = ee.Image(TimeSeriesList.get(k))
        for (var l = k+1; l < NumberOfItems-1; l++){
          var  nextImage = ee.Image(TimeSeriesList.get(l))
          var  Concordant = CurrentImage.lt(nextImage)
          var ConcordantArray = ConcordantArray.add(Concordant)
          var  Discordant = CurrentImage.gt(nextImage)
          var DiscordantArray = DiscordantArray.add(Discordant)}
    }
  var  ConcordantSum = ee.ImageCollection(ConcordantArray).sum().rename('csum')
  var  DiscordantSum = ee.ImageCollection(DiscordantArray).sum().rename('dsum')
  var  MKSstat = ConcordantSum.subtract(DiscordantSum)//.addBands([ConcordantSum, DiscordantSum])
    
  return MKSstat
}

var kendall_stats = kendall_test(ndviCollection);
// print(kendall_stats)


// Stretch this as necessary.

ui.root.clear()

var map1 = ui.Map()
var map2 = ui.Map()
var map3 = ui.Map()
var map4 = ui.Map()

map1.addLayer(trend.select('meanNDVI_tau'), {min:-1, max:1, palette: palette}, 'trend'); //palette: palette
map2.addLayer(kendall_stats, {min:-1, max:1, palette: palette}, 'MKtrend');
map1.addLayer(baseKendall, {min:-1, max:1,palette: palette}, 'kendall'); //palette: palette
map2.addLayer(reportingKendall, {min:-1, max:1, palette: palette}, 'kendall reporting')

var maps = [map1, map2, map3, map4]

var linker = ui.Map.Linker(maps);

// Create a grid of maps.
var mapGrid = ui.Panel(
    [
      ui.Panel([maps[0], maps[1]], null, {stretch: 'both'}),
      ui.Panel([maps[2], maps[3]], null, {stretch: 'both'})
    ],
    ui.Panel.Layout.Flow('horizontal'), {stretch: 'both'});

ui.root.widgets().reset([mapGrid]);
ui.root.setLayout(ui.Panel.Layout.Flow('horizontal'));

// var trajectory = collReporting.reduce(ee.Reducer.kendallsCorrelation());
// Map.addLayer(trajectory.select('meanNDVI_tau'), {palette: palette}, 'trendReporting'); //palette: palette
// Map.addLayer(kendall, {palette: palette}, 'kendall'); //palette: palette
// Map.addLayer(kendallReporting, {palette:palette}, 'kendallReporting')

// //   =================      ADD Sen's Slope     ============== 

// var slope = function(i, j) { // i and j are images
//   return ee.Image(j).subtract(i)
//       .divide(ee.Image(j).metadata('index').difference(ee.Image(i).metadata('index'), 'years'))
//       .rename('slope')
//       .float();
// };

// var slopes = ee.ImageCollection(joined.map(function(current) {
//   var afterCollection = ee.ImageCollection.fromImages(current.get('after'));
//   return afterCollection.map(function(image) {
//       return ee.Image(slope(current, image));
//   });
// }).flatten());

// var sensSlope = slopes.reduce(ee.Reducer.median(), 2); // Set parallelScale.

// Map.addLayer(sensSlope, {palette: palette}, 'sensSlope');


// ==============    Get Sen's Intercept    ==================

// var epochDate = ee.Date('1970-01-01');
// var sensIntercept = ndviCollection.map(function(image) {
//   var epochDays = image.date().difference(epochDate, 'days').float();
//   return image.subtract(sensSlope.multiply(epochDays)).float();
// }).reduce(ee.Reducer.median(), 2);

// Map.addLayer(sensIntercept, {}, 'sensIntercept');


// // ==============    Kendall Variance   ==================

// // Values that are in a group (ties).  Set all else to zero.
// var groups = ndviCollection.map(function(i) {
//   var matches = ndviCollection.map(function(j) {
//     return i.eq(j); // i and j are images.
//   }).sum();
//   return i.multiply(matches.gt(1));
// });

// // Compute tie group sizes in a sequence.  The first group is discarded.
// var group = function(array) {
//   var length = array.arrayLength(0);
//   // Array of indices.  These are 1-indexed.
//   var indices = ee.Image([1])
//       .arrayRepeat(0, length)
//       .arrayAccum(0, ee.Reducer.sum())
//       .toArray(1);
//   var sorted = array.arraySort();
//   var left = sorted.arraySlice(0, 1);
//   var right = sorted.arraySlice(0, 0, -1);
//   // Indices of the end of runs.
//   var mask = left.neq(right)
//   // Always keep the last index, the end of the sequence.
//       .arrayCat(ee.Image(ee.Array([[1]])), 0);
//   var runIndices = indices.arrayMask(mask);
//   // Subtract the indices to get run lengths.
//   var groupSizes = runIndices.arraySlice(0, 1)
//       .subtract(runIndices.arraySlice(0, 0, -1));
//   return groupSizes;
// };

// // See equation 2.6 in Sen (1968).
// var factors = function(image) {
//   return image.expression('b() * (b() - 1) * (b() * 2 + 5)');
// };

// var groupSizes = group(groups.toArray());
// var groupFactors = factors(groupSizes);
// var groupFactorSum = groupFactors.arrayReduce('sum', [0])
//       .arrayGet([0, 0]);

// var count = joined.count();

// var kendallVariance = factors(count)
//     .subtract(groupFactorSum)
//     .divide(18)
//     .float();
// //Map.addLayer(kendallVariance, {}, 'kendallVariance');

// // ===========    Compute Z-statistics  ============


// var zero = kendall.multiply(kendall.eq(0));
// var pos = kendall.multiply(kendall.gt(0)).subtract(1);
// var neg = kendall.multiply(kendall.lt(0)).add(1);

// var z = zero
//     .add(pos.divide(kendallVariance.sqrt()))
//     .add(neg.divide(kendallVariance.sqrt()));
// //Map.addLayer(z, {min: -2, max: 2}, 'z');

// // https://en.wikipedia.org/wiki/Error_function#Cumulative_distribution_function
// function eeCdf(z) {
//   return ee.Image(0.5)
//       .multiply(ee.Image(1).add(ee.Image(z).divide(ee.Image(2).sqrt()).erf()));
// }

// function invCdf(p) {
//   return ee.Image(2).sqrt()
//       .multiply(ee.Image(p).multiply(2).subtract(1).erfInv());
// }
    
// //           ==============  Compute P-values ============

// var p = ee.Image(1).subtract(eeCdf(z.abs()));
// //Map.addLayer(p, {min: 0, max: 1}, 'p');

// // Pixels that can have the null hypothesis (there is no trend) rejected.
// // Specifically, if the true trend is zero, there would be less than 5%
// // chance of randomly obtaining the observed result (that there is a trend).
// //Map.addLayer(p, {}, 'significant trends');



//------------------------------------------------------------------------------------------
//      Export as GeoTIFF
//------------------------------------------------------------------------------------------

// if (country === 'Tunisia'){

//         Export.image.toDrive({
//         image: annualNDVI,
//         description: country + '_NDVI_' + suffix + '_' + Year,
//         scale: 30,
//         region: studyArea,
//         maxPixels:  1e13,
//         fileFormat: 'GeoTIFF',
//         folder:'GEE_classification',
//         formatOptions: {
//           cloudOptimized: true
//             },
//           skipEmptyTiles: true
//           });
// }

// if (country === 'Egypt'){
        
//         tileList.map(function(tile){
//               var tmpAOI = studyArea.filter(ee.Filter.eq('id', tile));
//               //print(tile)
//               Map.addLayer(tmpAOI,{},'ID',1);

//               Export.image.toDrive({
//                     image: annualNDVI,
//                     description: country +'_'+tile+'_NDVI_' + suffix + '_' + Year,
//                     scale: 30,
//                     region: tmpAOI,
//                     maxPixels:  1e13,
//                     fileFormat: 'GeoTIFF',
//                     folder:'GEE_classification',
//                     formatOptions: {
//                       cloudOptimized: true
//                         },
//                       skipEmptyTiles: true
//               })
//               //return 0
//             })
        
        
        
// }
  
  //Map.addLayer(L7.first(), {}, 'L7');
  
  */
