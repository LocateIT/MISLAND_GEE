var countries = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017");

var country = 'Egypt';

var studyArea = countries.filter(ee.Filter.eq('country_na',country ))

Map.centerObject(studyArea, 8);

var trend_baselinePeriod = ee.List.sequence(2000, 2015, 1).getInfo(); //baseline years for reporting trend 2000-2015 (16 years)
var trend_reportingPeriod = ee.List.sequence(2005, 2020, 1).getInfo(); //reporting period for trend 2005 - 2020 (16 years)

var state_baselinePeriod = ee.List.sequence(2000, 2012, 1).getInfo(); // first 13year for reporting state 2000 - 2012
var state_baselinePeriod_final = ee.List.sequence(2013, 2015, 1).getInfo(); // final 3 yars of the reporing period 2013 - 2015

var state_reportingPeriod = ee.List.sequence(2005, 2017, 1).getInfo(); // first 13year of the reporting period
var state_reportingPeriod_final = ee.List.sequence(2018, 2020, 1).getInfo(); // last 3 years of the roporting period


function yearlyNDVI(year){
  var start_date = year+ '-01-01';
  var end_date   = year+ '-12-31';
  
  var NDVI = ee.ImageCollection('MODIS/006/MOD13Q1')
                  .filter(ee.Filter.date(start_date, end_date))
                  .select('NDVI')
                  .mean().divide(10000)
                  .clip(studyArea);

  return NDVI
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
  
  var threshold_z_score = [100, 1.96, 1.28, -1.28, -1.96]
  var classiffied_Kendall = kendall.lt(threshold_z_score).reduce('sum').toInt()
  
  var SDG_threshold = [-1.96, 1.96, 100]
  var classiffied_SDG = kendall.lt(SDG_threshold).reduce('sum').toInt().rename('rclassed')
  
  return [classiffied_Kendall, classiffied_SDG, kendall]
}

var trendBaselinePeriod = mannKendall(trendBaseline_NDVIcollection);
var trendReportingPeriod = mannKendall(trendReporting_NDVIcollection);

var baselineTrend = trendBaselinePeriod[0]
var reportingTrend = trendReportingPeriod[0]
var trendSdg_baseline = trendBaselinePeriod[1]
var trendSdg_reporting = trendReportingPeriod[1]
var kendall_baseline = trendBaselinePeriod[2]
var kendall_reporting = trendReportingPeriod[2]

//var reportingKendall = mannKendall(collReporting);
Map.addLayer(baselineTrend, {min:1, max:5, palette:['1b8607','b8ff68','ffffff','ffc443','ff1919']}, 'Baseline p-value');
Map.addLayer(reportingTrend, {min:1, max:5, palette:['1b8607','b8ff68','ffffff','ffc443','ff1919']}, 'Reporting p-value');

Map.addLayer(trendSdg_baseline, {min:1, max:3, palette:['1b8607','fffdc6','ff1919']}, 'Baseline p-value reclassed');
Map.addLayer(trendSdg_reporting, {min:1, max:3, palette:['1b8607','fffdc6','ff1919']}, 'Reporting p-value reclassed');

// Map.addLayer( kendall_baseline, {palette:['1b8607','ffffff','ff1919']}, 'Baseline p-value reclassed');
// Map.addLayer( kendall_reporting, {palette:['1b8607','ffffff','ff1919']}, 'Reporting p-value reclassed');


/*

var ndviVis = {
  min: 0.0,
  max: 8000.0,
  palette: [
    'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
    '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
    '012E01', '011D01', '011301'
  ],
};
Map.centerObject(studyArea);
//Map.addLayer(ndvi, ndviVis, 'NDVI');

// Export.image.toCloudStorage({
// image:ndvi,
// description: 'NDVI',
// maxPixels:1e13,
// scale:250,
// bucket:'oss_ldms_vi',
// region:table 
// })

*/