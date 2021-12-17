import geopandas as gpd
import time
import os
import pandas as pd
import shapely
import math
from shapely.geometry import Point, Polygon

ee = None
alos = None
dir_save = None

def maskSentinel2(image):
    qa = image.select('QA60')
    cloudBitMask = 1 << 10
    cirrusBitMask = 1 << 11
    mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0))
    refflectance = image.divide(10000)
    swr = refflectance.select(['B11']).multiply(750)
    nir = refflectance.select(['B8A']).multiply(500)
    blu = refflectance.select(['B2']).multiply(500)
    clrmsk = swr.gt(120).Or(nir.gt(160)).And(blu.lt(120))
    sdwmsk = swr.lt(120).And(nir.lt(160))
    cldmsk = blu.gt(140)
    return refflectance.updateMask(clrmsk)

def maskSentinel2ToTopo(image):
    qa = image.select('QA60')
    cloudBitMask = 1 << 10
    cirrusBitMask = 1 << 11
    mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0))
    refflectance = image.divide(10000)
    degreeToRadian = 0.01745
    slope = ee.Terrain.slope(alos).multiply(degreeToRadian)
    slopen = ee.Terrain.slope(alos)
    aspect = ee.Terrain.aspect(alos)
    aspectn = ee.Terrain.aspect(alos).multiply(degreeToRadian)
    zenitangle = image.get('MEAN_SOLAR_ZENITH_ANGLE')
    z = ee.Image(ee.Number(zenitangle)).multiply(degreeToRadian)
    az = ee.Image(image.get('MEAN_SOLAR_AZMITH_ANGLE')).multiply(degreeToRadian)
    az2 = ee.Image(ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE')))
    delta = az2.subtract(aspect)
    deltaR = delta.multiply(degreeToRadian)
    ic1 = z.cos().multiply(slope.cos())
    ic2 = z.sin().multiply(slope.sin().multiply(deltaR.cos()))
    IC = ic1.add(ic2)
    topSentinel2 = []
    ytopo = refflectance.select(['B1'])
    y = refflectance.select(['B2','B3','B4','B8','B8A','B11','B12'])
    band = y.bandNames()
    beta = [0.0103,
          0.0222,
          0.0180,
          0.1848,
          0.1848,
          0.0903,
          0.0391]
    top3 = []
    for i in range(0,len(beta)):
        slopeReg = ee.Number(beta[i])
        bandn = ee.String(band.get(i))
        collect = refflectance.select(bandn)
        newimg = image.expression('(img - ((slopeR) * (IC - cos_zenit)))',{
                    'img':collect,
                    'slopeR':slopeReg,
                    'IC':IC,
                    'cos_zenit':z.cos()
                })
        ytopo = ytopo.addBands(newimg)
    swr = ytopo.select(['B11']).multiply(750)
    nir = ytopo.select(['B8A']).multiply(500)
    blu = ytopo.select(['B2']).multiply(500)
    clrmsk = swr.gt(120).Or(nir.gt(160)).And(blu.lt(120))
    sdwmsk = swr.lt(120).And(nir.lt(160))
    cldmsk = blu.gt(140)
    return ytopo.updateMask(clrmsk)
    
def creatingMapFromGridS2(AOI,start,end):
    #stretching faktor B2,3,4,8,8A ==> 500. 500 dipilih untuk menjadikan range value nya ada di 0-255 (range reflectance =>0-0.490)
    #stretching faktor B11 =>750 (0-0.368)
    #stretcng factor B12 =>1000 (0-0.293)
    
    s2Topo = ee.ImageCollection('COPERNICUS/S2').filterDate(start,end).filterBounds(AOI).map(
        maskSentinel2ToTopo).select(['B2','B3','B4', 'B8', 'B11','B12']);
    s2TopoMin = s2Topo.reduce(ee.Reducer.percentile([15])).clip(AOI);
    s2TopoMed = s2Topo.reduce(ee.Reducer.percentile([50])).clip(AOI);
    s2TopoMean = s2Topo.reduce(ee.Reducer.mean()).clip(AOI);
    
    B2min = s2TopoMin.select(['B2_p15']).multiply(500).rename('S2_B2min');
    B3min = s2TopoMin.select(['B3_p15']).multiply(500).rename('S2_B3min');
    B4min = s2TopoMin.select(['B4_p15']).multiply(500).rename('S2_B4min');
    B8min = s2TopoMin.select(['B8_p15']).multiply(500).rename('S2_B8min');
    B11min = s2TopoMin.select(['B11_p15']).multiply(750).rename('S2_B11min');
    B12min = s2TopoMin.select(['B12_p15']).multiply(750).rename('S2_B12min');

    B2med = s2TopoMed.select(['B2_p50']).multiply(500).rename('S2_B2med');
    B3med = s2TopoMed.select(['B3_p50']).multiply(500).rename('S2_B3med');
    B4med = s2TopoMed.select(['B4_p50']).multiply(500).rename('S2_B4med');
    B8med = s2TopoMed.select(['B8_p50']).multiply(500).rename('S2_B8med');
    B11med = s2TopoMed.select(['B11_p50']).multiply(750).rename('S2_B11med');
    B12med = s2TopoMed.select(['B12_p50']).multiply(750).rename('S2_B12med');

    B2mean = s2TopoMean.select(['B2_mean']).multiply(500).rename('S2_B2mean');
    B3mean = s2TopoMean.select(['B3_mean']).multiply(500).rename('S2_B3mean');
    B4mean = s2TopoMean.select(['B4_mean']).multiply(500).rename('S2_B4mean');
    B8mean = s2TopoMean.select(['B8_mean']).multiply(500).rename('S2_B8mean');
    B11mean = s2TopoMean.select(['B11_mean']).multiply(500).rename('S2_B11mean');
    B12mean = s2TopoMean.select(['B12_mean']).multiply(750).rename('S2_B12mean');
   
    sentinel2_8bit = ee.Image.cat(
                        #B2min,B3min,B4min,B8min,B11min,B12min,
                        B2med,B3med,B4med,B8med,B11med,B12med,
                        B2mean,B3mean,B4mean,B8mean,B11mean,B12mean).toFloat();
    return sentinel2_8bit

def indices_collection_rev(image):
    ndvi=image.expression('(nir-red)/(nir+red)',{
        'nir':image.select('B8'),
        'red':image.select('B4')
    }).rename('NDVI')
    mndvi=image.expression('(nir-red)/(nir+red-(2*aer))',{
        'nir':image.select('B8'),
        'red':image.select('B4'),
        'aer':image.select('B1')
    }).rename('mNDVI')
    ndwi=image.expression('(green-nir)/(nir+green)',{
        'nir':image.select('B8'),
        'green':image.select('B3')
    }).rename('NDWI')
    mndwi=image.expression('(green-swir)/(swir+green)',{
        'swir':image.select('B11'),
        'green':image.select('B3')
    }).rename('mNDWI')
    ndbi=image.expression('(swir-nir)/(swir+nir)',{
        'nir':image.select('B8'),
        'swir':image.select('B11')
    }).rename('NDBI')
    savi=image.expression('(1.5)*((nir-red)/(nir+red+0.5))',{
        'nir':image.select('B8'),
        'red':image.select('B4')
    }).rename('SAVI')
    arvi=image.expression('(nir-(red-(1*(blue-red))))/(nir+(red-(1*(blue-red))))',{
        'nir':image.select('B8'),
        'red':image.select('B4'),
        'blue':image.select('B2')
    }).rename('ARVI')
    bright=image.expression('(0.3510*blue)+(0.3813*green)+(0.3437*red)+(nir*0.7196)+(swir*0.2396)+(0.1949*swir2)',{
        'nir':image.select('B8'),
        'red':image.select('B4'),
        'blue':image.select('B2'),
        'green':image.select('B3'),
        'swir':image.select('B11'),
        'swir2':image.select('B12')
    }).rename('bright') 
    green=image.expression('(-0.3599*blue)+(-0.3533*green)+(-0.4734*red)+(nir*0.6633)+(swir*0.0087)+(-0.2856*swir2)',{
        'nir':image.select('B8'),
        'red':image.select('B4'),
        'blue':image.select('B2'),
        'green':image.select('B3'),
        'swir':image.select('B11'),
        'swir2':image.select('B12')
    }).rename('green')
    wet=image.expression('(0.2578*blue)+(0.2305*green)+(0.0883*red)+(nir*0.1071)+(swir*-0.7611)+(-0.5308*swir2)',{
        'nir':image.select('B8'),
        'red':image.select('B4'),
        'blue':image.select('B2'),
        'green':image.select('B3'),
        'swir':image.select('B11'),
        'swir2':image.select('B12')
    }).rename('wet')
    return image.addBands(ndvi).addBands(mndvi).addBands(ndwi).addBands(mndwi).addBands(ndbi).addBands(savi).addBands(arvi).addBands(bright).addBands(green).addBands(wet)
    
def creatingMapFromGridS2_indices(AOI,start,end):
    
    s2Topo = ee.ImageCollection('COPERNICUS/S2').filterDate(start,end).filterBounds(AOI).map(
        maskSentinel2ToTopo).map(indices_collection_rev).select(['NDVI','mNDVI','NDWI','mNDWI','NDBI','SAVI','ARVI','bright','green','wet']);
    s2TopoMin = s2Topo.reduce(ee.Reducer.percentile([15])).clip(AOI);
    s2TopoMed = s2Topo.reduce(ee.Reducer.percentile([50])).clip(AOI);
    s2TopoMean = s2Topo.reduce(ee.Reducer.mean()).clip(AOI);
    sentinel2_indices_8bit=ee.Image.cat( s2TopoMed,s2TopoMean).toFloat();
    return sentinel2_indices_8bit
    
def creatingMapKombinasi(AOI,start,end):
    sentinel2_8bit = creatingMapFromGridS2(AOI,start,end)
    sentinel2_8bit_indices = creatingMapFromGridS2_indices(AOI,start,end)
    alos_ = alos.select('AVE_DSM').rename('alos_dsm').toFloat();
    slope = ee.Terrain.slope(alos).select('slope').rename('alos_slope').toFloat();
    landform=ee.Image("CSP/ERGo/1_0/Global/ALOS_landforms").select('constant').rename('alos_landform').toFloat();
    kombinasi = ee.Image.cat([sentinel2_8bit,sentinel2_8bit_indices,alos_,slope,landform])
    return kombinasi
    
def get_queued_task():
    queued_task_count = 0
    for queued_task in ee.batch.Task.list():
        if queued_task.state in ["READY","RUNNING"]:
            queued_task_count += 1
    return queued_task_count

def get_queued_task_filenames():
    task_filenames = []
    for queued_task in ee.batch.Task.list():
        if queued_task.state in ["READY","RUNNING"]:
            task_filenames.append(queued_task.status()['description'])
    return  task_filenames
    
def download_satellite_imagery(sat_imagery,citra_name,AOI,max_queue,counter_grid):
    task_count = get_queued_task()
    queued_filenames = get_queued_task_filenames()
    
    if counter_grid == 200 : 
        print("Proses download sudah mencapai 100 GRID, proses akan berhenti selama 2 menit.")
        time.sleep(60)
        counter_grid = 0
    
    if task_count == max_queue:
        while task_count > 3:
            time.sleep(5)
            task_count = get_queued_task() 

        if citra_name not in queued_filenames:
            print("Downloading " + citra_name+'_sentinel2_10m_med_mean_only')

            task = ee.batch.Export.image.toDrive(
                image= sat_imagery,
                folder=dir_save,
                fileNamePrefix= citra_name + '_sentinel2_10m_med_mean_rev',
                scale= 10,
                region= AOI,
                maxPixels= 1e10,
                description= citra_name
            );

            task.start()  
            #counter_grid+=1
    else:
        if (citra_name not in queued_filenames):
            print("--------------------")
            print("Starting new task...")
            print("downloading " + citra_name+'_sentinel2_10m_med_mean_only')

            task = ee.batch.Export.image.toDrive(
                image= sat_imagery,
                folder=dir_save,
                fileNamePrefix= citra_name + '_sentinel2_10m_med_mean_rev',
                scale= 10,
                region= AOI,
                maxPixels= 1e10,
                description= citra_name
            );

            task.start()
            counter_grid+=1
            
def mulai_download(path_file_grid,max_queue,counter_grid,list_id,ee_init):
    grid = gpd.read_file(path_file_grid)
    periode = {
          'setahun' : {
              'awal':'2020-01-01',
              'akhir':'2020-12-31',
          }
      }
      
    rows=grid.shape[0]
       
    for i in range(0,grid.shape[0]):
        for j in periode:
            if grid['ID_GRID'][i] in list_id :
                citra_name = str(j)+ '_' + str(grid['ID_GRID'][i]) 
                
                geom = [[grid['left'][i],grid['top'][i]],
                        [grid['right'][i],grid['top'][i]],
                        [grid['right'][i],grid['bottom'][i]],
                        [grid['left'][i],grid['bottom'][i]]]
                d = periode[j]
                AOI = ee.Geometry.Polygon(geom)
                hasil = creatingMapKombinasi(AOI,d['awal'],d['akhir'])
                download_satellite_imagery(hasil,citra_name,AOI,max_queue,counter_grid)