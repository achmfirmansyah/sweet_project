var alos = ee.Image("JAXA/ALOS/AW3D30/V2_2");
//Daerah Konawe Selatan
var AOI= ee.Geometry.Polygon(
        [[[121.97, -3.98],
          [121.97, -4.54],
          [122.91, -4.54],
          [122.91, -3.98]]], null, false);
// ###########################################################
// Landsat 8
// ###########################################################
// Function to mask clouds using the quality band of Landsat 8.
var maskL8 = function(image) {
  var qa = image.select('BQA');
  var mask = qa.bitwiseAnd(1 << 4).eq(0);
  return image.updateMask(mask);
}

var IC2Ref = function(img){
  var degree2radian = 0.01745;
  var slope = ee.Terrain.slope(alos).multiply(degree2radian);
  var slopen = ee.Terrain.slope(alos);
  var aspect = ee.Terrain.aspect(alos);
  var aspectn = ee.Terrain.aspect(alos).multiply(degree2radian);
  var nilai = ee.Number(90)
  var zenitangle = nilai.subtract(ee.Number(img.get('SUN_ELEVATION')))
  var z = ee.Image(ee.Number(zenitangle)).multiply(degree2radian);
  var az = ee.Image(ee.Number(img.get('SUN_AZIMUTH'))).multiply(degree2radian);
  var az2 = ee.Image(ee.Number(img.get('SUN_AZIMUTH')));
  var delta = (az2.subtract(aspect));
  var deltaR = delta.multiply(degree2radian)
  // IC
  var ic1 = z.cos().multiply(slope.cos());  // azimut diganti zenith
  var ic2 = z.sin().multiply(slope.sin().multiply(deltaR.cos())); // azimut diganti zenith
  var IC = ic1.add(ic2);
  var swr1 = img.select(['B6']);
  var y = img.select(['B2','B3','B4','B5','B6','B7']);
  var band= y.bandNames();
  var beta = [0.0103,0.0222,0.0180,0.1848,0.0903,0.0391];
  var top2 = [] ; 
  var i;
  for (i = 0; i < 6; i++) {
    var slopeReg = (ee.Number(beta[i]));
    var bandn = ee.String(band.get(i));
    var collect = img.select(bandn);
    var newimg = img.expression(
    '(img - ((slopeR) * (IC - cos_zenit)))',
    {
    'img'       : collect,
    'slopeR'    : slopeReg,
    'IC'        : IC,
    'cos_zenit' :z.cos()
    });
    top2.push(newimg);
  } 
  var Topo654 = ee.Image(top2);
return ee.Image(Topo654);
};
//get LANDSAT8 Data
var L8TopoCol = ee.ImageCollection("LANDSAT/LC08/C01/T1_TOA")
        .filterDate("2019-01-01","2019-12-30")
        .filterBounds(AOI)
        .map(maskL8)
        .map(IC2Ref)
        .select(['B2','B3','B4', 'B5', 'B6']);
//reducer
var L8TopoMin = L8TopoCol.reduce(ee.Reducer.percentile([15])).clip(AOI);
var L8TopoMean = L8TopoCol.reduce(ee.Reducer.mean()).clip(AOI);
var L8TopoMax = L8TopoCol.reduce(ee.Reducer.percentile([85])).clip(AOI);
var L8TopoMed = L8TopoCol.reduce(ee.Reducer.percentile([50])).clip(AOI);
var L8TopoStDev = L8TopoCol.reduce(ee.Reducer.stdDev()).clip(AOI);

var b2min = L8TopoMin.select(['B2_p15']).multiply(1000).rename('l8_B2min');
var b3min = L8TopoMin.select(['B3_p15']).multiply(1000).rename('l8_B3min');
var b4min = L8TopoMin.select(['B4_p15']).multiply(1000).rename('l8_B4min');
var b5min = L8TopoMin.select(['B5_p15']).multiply(500).rename('l8_B5min');
var b6min = L8TopoMin.select(['B6_p15']).multiply(500).rename('l8_B6min');

var b2max = L8TopoMax.select(['B2_p85']).multiply(1000).rename('l8_B2max');
var b3max = L8TopoMax.select(['B3_p85']).multiply(1000).rename('l8_B3max');
var b4max = L8TopoMax.select(['B4_p85']).multiply(1000).rename('l8_B4max');
var b5max = L8TopoMax.select(['B5_p85']).multiply(500).rename('l8_B5max');
var b6max = L8TopoMax.select(['B6_p85']).multiply(500).rename('l8_B6max');

var b2med = L8TopoMed.select(['B2_p50']).multiply(1000).rename('l8_B2med');
var b3med = L8TopoMed.select(['B3_p50']).multiply(1000).rename('l8_B3med');
var b4med = L8TopoMed.select(['B4_p50']).multiply(1000).rename('l8_B4med');
var b5med = L8TopoMed.select(['B5_p50']).multiply(500).rename('l8_B5med');
var b6med = L8TopoMed.select(['B6_p50']).multiply(500).rename('l8_B6med');

var b2mean = L8TopoMean.select(['B2_mean']).multiply(1000).rename('l8_B2mean');
var b3mean = L8TopoMean.select(['B3_mean']).multiply(1000).rename('l8_B3mean');
var b4mean = L8TopoMean.select(['B4_mean']).multiply(1000).rename('l8_B4mean');
var b5mean = L8TopoMean.select(['B5_mean']).multiply(500).rename('l8_B5mean');
var b6mean = L8TopoMean.select(['B6_mean']).multiply(500).rename('l8_B6mean');

var b2stddev = L8TopoStDev.select(['B2_stdDev']).multiply(1000).rename('l8_B2stddev');
var b3stddev = L8TopoStDev.select(['B3_stdDev']).multiply(1000).rename('l8_B3stddev');
var b4stddev = L8TopoStDev.select(['B4_stdDev']).multiply(1000).rename('l8_B4stddev');
var b5stddev = L8TopoStDev.select(['B5_stdDev']).multiply(500).rename('l8_B5stddev');
var b6stddev = L8TopoStDev.select(['B6_stdDev']).multiply(500).rename('l8_B6stddev');

//cat toL8bit
var L816bit = ee.Image.cat(
        b2min,
        b3min,
        b4min,
        b5min,
        b6min,
        b2med,
        b3med,
        b4med,
        b5med,
        b6med,
        b2max,
        b3max,
        b4max,
        b5max,
        b6max,
        b2mean,
        b3mean,
        b4mean,
        b5mean,
        b6mean,
        b2stddev,
        b3stddev,
        b4stddev,
        b5stddev,
        b6stddev);
        
var L88bit = L816bit.uint8();
print('ini data L8 Topo 8bit ', L88bit);

var displayParamsL8truecolor_med = {
  bands: ['l8_B6med', 'l8_B5med', 'l8_B4med'],
  min: [0,0,0],
  max: [255,255,255],
  //gamma: [1, 1, 1]
};
var displayParamsL8truecolor_min = {
  bands: ['l8_B6min', 'l8_B5min', 'l8_B4min'],
  min: [0,0,0],
  max: [255,255,255],
  //gamma: [1, 1, 1]
};
var displayParamsL8truecolor_max = {
  bands: ['l8_B6max', 'l8_B5max', 'l8_B4max'],
  min: [0,0,0],
  max: [255,255,255],
  //gamma: [1, 1, 1]
};
var displayParamsL8truecolor_mean = {
  bands: ['l8_B6mean', 'l8_B5mean', 'l8_B4mean'],
  min: [0,0,0],
  max: [255,255,255],
  //gamma: [1, 1, 1]
};
var displayParamsL8truecolor_stddev = {
  bands: ['l8_B6stddev', 'l8_B5stddev', 'l8_B4stddev'],
  min: [0,0,0],
  max: [255,255,255],
  //gamma: [1, 1, 1]
};
Map.addLayer(L88bit, displayParamsL8truecolor_med, 'L8 Topo 654_med');
Map.addLayer(L88bit, displayParamsL8truecolor_min, 'L8 Topo 654_min');
Map.addLayer(L88bit, displayParamsL8truecolor_max, 'L8 Topo 654_max');
Map.addLayer(L88bit, displayParamsL8truecolor_mean, 'L8 Topo 654_mean');
Map.addLayer(L88bit, displayParamsL8truecolor_stddev, 'L8 Topo 654_stddev');
Map.centerObject(AOI);

//Histogram view
// Pre-define some customization options.
var options = {
  title: 'Landsat 8 DN histogram, bands 654',
  fontSize: 8,
  hAxis: {title: 'DN'},
  vAxis: {title: 'count of DN'},
  series: {
    0: {color: 'blue'},
    1: {color: 'green'},
    2: {color: 'red'},
    3: {color: 'orange'},
    4: {color: 'purple'}}
};

// Make the histogram, set the options.
var histogram = ui.Chart.image.histogram(L88bit.select(['l8_B2mean', 'l8_B3mean', 'l8_B4mean','l8_B5mean', 'l8_B6mean']), AOI, 300)
    .setSeriesNames(['l8_B2_mean', 'l8_B3_mean', 'l8_B4_mean','l8_B5_mean', 'l8_B6_mean'])
    .setOptions(options);

// Display the histogram.
print('Histogram mean: ',histogram);

// Make the histogram, set the options.
var histogram_med = ui.Chart.image.histogram(L88bit.select(['l8_B2med', 'l8_B3med', 'l8_B4stddev','l8_B5med', 'l8_B6med']), AOI, 300)
    .setSeriesNames(['l8_B2_med', 'l8_B3_med', 'l8_B4_med','l8_B5_med', 'l8_B6_med'])
    .setOptions(options);

// Display the histogram.
print('Histogram med: ',histogram_med);


// Make the histogram, set the options.
var histogram_stdDev = ui.Chart.image.histogram(L88bit.select(['l8_B2stddev', 'l8_B3stddev', 'l8_B4stddev','l8_B5stddev', 'l8_B6stddev']), AOI, 300)
    .setSeriesNames(['l8_B2_stddev', 'l8_B3_stddev', 'l8_B4_stddev','l8_B5_stddev', 'l8_B6_stddev'])
    .setOptions(options);

// Display the histogram.
print('Histogram stdDev: ',histogram_stdDev);

// ###########################################################
// Sentinel 1
// ###########################################################
var slope_correction = function (collection,
                                 options
                                 ){
    // set defaults if undefined options
    options = options || {};
    var model = options.model || 'volume';
    var elevation = options.elevation || ee.Image('USGS/SRTMGL1_003');
    var buffer = options.buffer || 75;
    // we need a 90 degree in radians image for a couple of calculations
    var ninetyRad = ee.Image.constant(90).multiply(Math.PI/180);

    // Volumetric Model Hoekman 1990
    function _volume_model(theta_iRad, alpha_rRad){
      var nominator = (ninetyRad.subtract(theta_iRad).add(alpha_rRad)).tan();
      var denominator = (ninetyRad.subtract(theta_iRad)).tan();
      return nominator.divide(denominator);
    }
    // surface model Ulander et al. 1996
    function _surface_model(theta_iRad, alpha_rRad, alpha_azRad){
      var nominator = (ninetyRad.subtract(theta_iRad)).cos();
      var denominator = alpha_azRad.cos()
        .multiply((ninetyRad.subtract(theta_iRad).add(alpha_rRad)).cos());
      return nominator.divide(denominator);
    }
    // buffer function (thanks Noel)
    function _erode(img, distance) {
      var d = (img.not().unmask(1)
          .fastDistanceTransform(30).sqrt()
          .multiply(ee.Image.pixelArea().sqrt()));
      return img.updateMask(d.gt(distance));
    }
    // calculate masks
    function _masking(alpha_rRad, theta_iRad, proj, buffer){
        // layover, where slope > radar viewing angle
        var layover = alpha_rRad.lt(theta_iRad).rename('layover');
        // shadow
        var shadow = alpha_rRad.gt(ee.Image.constant(-1).multiply(ninetyRad.subtract(theta_iRad))).rename('shadow');
        // combine layover and shadow
        var mask = layover.and(shadow);
        // add buffer to final mask
        if (buffer > 0)
            mask = _erode(mask, buffer);
        return mask.rename('no_data_mask');
   }
    function _correct(image){
        // get image geometry and projection
        var geom = image.geometry();
        var proj = image.select(1).projection();
        // get look direction angle
        var heading = (ee.Terrain.aspect(
            image.select('angle')).reduceRegion(ee.Reducer.mean(), geom, 1000).get('aspect')
            );
        // Sigma0 to Power of input image
        var sigma0Pow = ee.Image.constant(10).pow(image.divide(10.0));
        // Radar geometry
        var theta_iRad = image.select('angle').multiply(Math.PI/180).clip(geom);
        var phi_iRad = ee.Image.constant(heading).multiply(Math.PI/180);
        // Terrain geometry
        var alpha_sRad = ee.Terrain.slope(elevation).select('slope')
            .multiply(Math.PI/180).setDefaultProjection(proj).clip(geom);
        var phi_sRad = ee.Terrain.aspect(elevation).select('aspect')
            .multiply(Math.PI/180).setDefaultProjection(proj).clip(geom);
        // Model geometry
        //reduce to 3 angle
        var phi_rRad = phi_iRad.subtract(phi_sRad);
        // slope steepness in range
        var alpha_rRad = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan();
        // slope steepness in azimuth
        var alpha_azRad = (alpha_sRad.tan().multiply(phi_rRad.sin())).atan();
        // Gamma_nought
        var gamma0 = sigma0Pow .divide(theta_iRad.cos());
               // models
        if (model == 'volume')
          var corrModel = _volume_model(theta_iRad, alpha_rRad);
        // if (model == 'surface')
        //   var corrModel = _surface_model(theta_iRad, alpha_rRad, alpha_azRad);
        // if (model == 'direct')
        //   var corrModel = _direct_model(theta_iRad, alpha_rRad, alpha_azRad);
        // apply model to derive gamma0_flat
        var gamma0_flat = gamma0.divide(corrModel);
        // transform to dB-scale
        var gamma0_flatDB = (ee.Image.constant(10)
            .multiply(gamma0_flat.log10()).select(['VV', 'VH'])
            );
        var angle_ = image.select('angle');
        // get Layover/Shadow mask
        var mask = _masking(alpha_rRad, theta_iRad, proj, buffer);
        var results = gamma0_flatDB.addBands(angle_);
        // return gamma_flat plus mask
        return results.addBands(mask).copyProperties(image);
    }
    // run correction function and return corrected collection
    return collection.map(_correct);
};

function func_topovolume(imagecollection, start, end){
  var s1_collection = imagecollection
    .filterBounds(AOI)
    .filterDate(start, end)
    .filter(ee.Filter.and(
      ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'),
      ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'),
      ee.Filter.eq('instrumentMode', 'IW'),
      ee.Filter.eq('orbitProperties_pass', 'DESCENDING')
       ));
  var volume = slope_correction(
    s1_collection
    );
    
  function maskangle(image) {
    var IA = image.select('angle');
    var edgeka = IA.lt(30.7); //.lt(30.3);
    var edgeki = IA.gt(44.4); //.lt(44.8);
    var maskedImage = image.mask().and(edgeki.not()).and(edgeka.not());
    var masknodata = image.select('no_data_mask');
    return image.updateMask(maskedImage).updateMask(masknodata);
  }
  var Metrixmed = volume
        .map(maskangle)
        .median()
        .select(['VV','VH']);
  var Metrixmin = volume
        .map(maskangle)
        .min()
        .select(['VV','VH']);
  var Metrixmax = volume
        .map(maskangle)
        .max()
        .select(['VV','VH']);
  var S116bit = ee.Image.cat(
        Metrixmin.select('VV').rename('s1_VVmin'),
        Metrixmin.select('VH').rename('s1_VHmin'),
        Metrixmed.select('VV').rename('s1_VVmed'),
        Metrixmed.select('VH').rename('s1_VHmed'),
        Metrixmax.select('VV').rename('s1_VVmax'),
        Metrixmax.select('VH').rename('s1_VHmax'));
  var S116bitclip = S116bit.clip(AOI);        
  return S116bitclip;
}
var SEN1 = ee.ImageCollection('COPERNICUS/S1_GRD');
var Topo = func_topovolume(SEN1, '2019-01-01', '2019-12-31');
var ratio = Topo.select('s1_VVmed').subtract(Topo.select('s1_VHmed')).rename('s1_ratiomed');
var RGBmedSen1 = ee.Image.cat(
         Topo.select('s1_VVmed'),
         Topo.select('s1_VHmed'),
         ratio);
print(RGBmedSen1, 'RGBmedSen1');
var Channels = {
  bands: ['s1_VVmed', 's1_VHmed', 's1_ratiomed'],
  min: [-15,-26,5.0],
  max: [0,0,15.5],
  gamma: 1.0,
};
Map.addLayer(RGBmedSen1, Channels, 'Image Sentinel-1 RGB');
var Channels2 = {
  bands: ['s1_VHmin', 's1_VHmed', 's1_VHmax'],
  min: [-15,-15,-15],
  max: [0,0,0],
  gamma: 1.7,
};
Map.addLayer(Topo, Channels2, 'Image Sentinel-1 MinMedMax');
// Pre-define some customization options.
var optionsS1 = {
  title: 'S1 Histogram',
  fontSize: 8,
  hAxis: {title: 'DN'},
  vAxis: {title: 'count of DN'},
  series: {
    0: {color: 'blue'},
    1: {color: 'green'},}
};

// Make the histogram, set the options.
var histogramS1 = ui.Chart.image.histogram(Topo.select(['s1_VHmed', 's1_VVmed']), AOI, 300)
    .setSeriesNames(['s1_VHmed', 's1_VVmed'])
    .setOptions(optionsS1);
print('Histogram S1',histogramS1)


// ###########################################################
// Sentinel 2
// ###########################################################
var masks2 = function(image) { 
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  var refflectance = image.divide(10000);
  var swr = refflectance.select(['B11']).multiply(1000);
  var nir = refflectance.select(['B8A']).multiply(1000);
  var blu = refflectance.select(['B2']).multiply(500);
  var CLRMSK = swr.gt(120).or(nir.gt(160)).and(blu.lt(120));
  var SDWMSK = swr.lt(120).and(nir.lt(160));
  var CLDMSK = blu.gt(140);
  return refflectance.updateMask(CLRMSK);
};

var masks2topo = function(image) { 
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  // https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR#bands
  var refflectance = image.divide(10000);
  var degree2radian = 0.01745;
  var slope = ee.Terrain.slope(alos).multiply(degree2radian);
  var slopen = ee.Terrain.slope(alos);
  var aspect = ee.Terrain.aspect(alos);
  var aspectn = ee.Terrain.aspect(alos).multiply(degree2radian);
  var zenitangle = image.get('MEAN_SOLAR_ZENITH_ANGLE');
  var z = ee.Image(ee.Number(zenitangle)).multiply(degree2radian);
  var az = ee.Image(image.get('MEAN_SOLAR_AZIMUTH_ANGLE')).multiply(degree2radian);
  var az2 = ee.Image(ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE')));
  var delta = (az2.subtract(aspect));
  var deltaR = delta.multiply(degree2radian);
  // IC
  var ic1 = z.cos().multiply(slope.cos());  // azimut diganti zenith
  var ic2 = z.sin().multiply(slope.sin().multiply(deltaR.cos())); // azimut diganti zenith
  var IC = ic1.add(ic2);
  var tops2 = [] ; 
  var ytopo = refflectance.select(['B1']);
  var y = refflectance.select(['B2','B3','B4', 'B8', 'B8A', 'B11','B12']);
  var band= y.bandNames();
  var beta = [0.0103,0.0222,0.0180,0.1848,0.1848,0.0903,0.0391];
  var top3 = [] ; 
  var i;
  for (i = 0; i < 7; i++) {
        var slopeReg = (ee.Number(beta[i]));
        var bandn = ee.String(band.get(i));
        var collect = refflectance.select(bandn);
        var newimg = image.expression(
              '(img - ((slopeR) * (IC - cos_zenit)))',
              {
              'img'       : collect,
              'slopeR'    : slopeReg,
              'IC'        : IC,
              'cos_zenit' :z.cos()
              // 'C': C
              });
        var ytopo=ytopo.addBands(newimg);
        }
  var swr = ytopo.select(['B11']).multiply(1000);
  var nir = ytopo.select(['B8A']).multiply(1000);
  var blu = ytopo.select(['B2']).multiply(500);
  var CLRMSK = swr.gt(120).or(nir.gt(160)).and(blu.lt(120));
  var SDWMSK = swr.lt(120).and(nir.lt(160));
  var CLDMSK = blu.gt(140);
  return ytopo.updateMask(CLRMSK);
};

var S2TopoCol= ee.ImageCollection('COPERNICUS/S2')
        .filterDate('2019-01-01', '2019-12-31')
        .filterBounds(AOI)
        //.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 100))
        .map(masks2topo)
        .select(['B2','B3','B4', 'B8', 'B11']);
        
var S2TopoMin = S2TopoCol.reduce(ee.Reducer.percentile([15])).clip(AOI);
var S2TopoMax = S2TopoCol.reduce(ee.Reducer.percentile([85])).clip(AOI);
var S2TopoMed = S2TopoCol.reduce(ee.Reducer.percentile([50])).clip(AOI);
var S2TopoMean = S2TopoCol.reduce(ee.Reducer.mean()).clip(AOI);
var S2TopoStdDev = S2TopoCol.reduce(ee.Reducer.stdDev()).clip(AOI);

var b2min = S2TopoMin.select(['B2_p15']).multiply(1000).rename('s2_B2min');
var b3min = S2TopoMin.select(['B3_p15']).multiply(1000).rename('s2_B3min');
var b4min = S2TopoMin.select(['B4_p15']).multiply(1000).rename('s2_B4min');
var b8min = S2TopoMin.select(['B8_p15']).multiply(500).rename('s2_B8min');
var b11min = S2TopoMin.select(['B11_p15']).multiply(500).rename('s2_B11min');

var b2max = S2TopoMax.select(['B2_p85']).multiply(1000).rename('s2_B2max');
var b3max = S2TopoMax.select(['B3_p85']).multiply(1000).rename('s2_B3max');
var b4max = S2TopoMax.select(['B4_p85']).multiply(1000).rename('s2_B4max');
var b8max = S2TopoMax.select(['B8_p85']).multiply(500).rename('s2_B8max');
var b11max = S2TopoMax.select(['B11_p85']).multiply(500).rename('s2_B11max');

var b2med = S2TopoMed.select(['B2_p50']).multiply(1000).rename('s2_B2med');
var b3med = S2TopoMed.select(['B3_p50']).multiply(1000).rename('s2_B3med');
var b4med = S2TopoMed.select(['B4_p50']).multiply(1000).rename('s2_B4med');
var b8med = S2TopoMed.select(['B8_p50']).multiply(500).rename('s2_B8med');
var b11med = S2TopoMed.select(['B11_p50']).multiply(500).rename('s2_B11med');

var b2mean = S2TopoMean.select(['B2_mean']).multiply(1000).rename('s2_B2mean');
var b3mean = S2TopoMean.select(['B3_mean']).multiply(1000).rename('s2_B3mean');
var b4mean = S2TopoMean.select(['B4_mean']).multiply(1000).rename('s2_B4mean');
var b8mean = S2TopoMean.select(['B8_mean']).multiply(500).rename('s2_B8mean');
var b11mean = S2TopoMean.select(['B11_mean']).multiply(500).rename('s2_B11mean');

var b2stddev = S2TopoStdDev.select(['B2_stdDev']).multiply(1000).rename('s2_B2stddev');
var b3stddev = S2TopoStdDev.select(['B3_stdDev']).multiply(1000).rename('s2_B3stddev');
var b4stddev= S2TopoStdDev.select(['B4_stdDev']).multiply(1000).rename('s2_B4stddev');
var b8stddev = S2TopoStdDev.select(['B8_stdDev']).multiply(500).rename('s2_B8stddev');
var b11stddev = S2TopoStdDev.select(['B11_stdDev']).multiply(500).rename('s2_B11stddev');

var S216bit = ee.Image.cat(
        b2min,
        b3min,
        b4min,
        b8min,
        b11min,
        b2med,
        b3med,
        b4med,
        b8med,
        b11med,
        b2max,
        b3max,
        b4max,
        b8max,
        b11max,
        b2mean,
        b3mean,
        b4mean,
        b8mean,
        b11mean,
        b2stddev,
        b3stddev,
        b4stddev,
        b8stddev,
        b11stddev);
        
var S28bit = S216bit.uint8();
print('ini data S2 Topo 8bit ', S28bit);


var displayParamsS2 = {
  bands: ['s2_B11med', 's2_B8med', 's2_B4med'],
  min: [0,0,0],
  max: [255,255,255],
  //gamma: [1, 1, 1]
};
Map.addLayer(S28bit, displayParamsS2,'Sentinel 2 8bit True Color')
// Make the histogram, set the options.

var optionss2 = {
  title: 'Sentinel 2 DN histogram',
  fontSize: 8,
  hAxis: {title: 'DN'},
  vAxis: {title: 'count of DN'},
  series: {
    0: {color: 'blue'},
    1: {color: 'green'},
    2: {color: 'red'},
    3: {color: 'orange'},
    4: {color: 'purple'}}
};
var histogramS2 = ui.Chart.image.histogram(S28bit.select(['s2_B2med', 's2_B3med','s2_B4med','s2_B8med','s2_B11med']), AOI, 300)
    .setSeriesNames(['s2_B2med', 's2_B3med','s2_B4med','s2_B8med','s2_B11med'])
    .setOptions(optionss2);
print('Histogram S2',histogramS2)

var comb=ee.Image.cat([L88bit,S28bit])
print(comb);
/*
Export.image.toDrive({
  image: comb,
  description: 's2_l8_konawe_selatan',
  scale: 20,
  region: AOI,
  maxPixels: 1e10
});
Export.image.toDrive({
  image: Topo,
  description: 's1_konawe_selatan',
  scale: 20,
  region: AOI,
  maxPixels: 1e10
});
*/
print(Topo)