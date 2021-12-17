## AUTHOR: ACHMAD FAUZI BAGUS FIRMANSYAH
import rioxarray as rio
import geopandas as gpd
import pandas as pd
from tqdm.notebook import tqdm_notebook
from google.cloud import storage
import os
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import Normalizer
from sklearn.decomposition import PCA
import plotly.express as px
import optuna
from catboost import Pool, EShapCalcType, EFeaturesSelectionAlgorithm,CatBoostClassifier
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
import pickle5 as pickle
from sklearn.metrics import classification_report,confusion_matrix,plot_confusion_matrix
from google.cloud import storage
import matplotlib.pyplot as plt
import seaborn as sns
import shap

from tqdm.notebook import tqdm
from rasterio.merge import merge
import rasterio
from geocube.api.core import make_geocube
from geocube.rasterize import rasterize_points_griddata, rasterize_points_radial
train_x_1=None
train_y_1=None
valid_x_1=None
valid_y_1=None
feature_select_1=None

def iterate_raster(folder, list_imagery,train_vector):
    data_=pd.DataFrame()
    for n_ in list_imagery:
        print('Creating training data for: ',n_)
        data_riox=rio.open_rasterio('kerangka-spasial/citra-satelit/citra-satelit_sentinel_10_m_'+folder+'_'+n_)
        r_=iterate_geom_(data_riox,train_vector)
        data_=r_.append(data_)
    print('Finish')
    return(data_)

def preprocessing_tif_vector(data_riox,geometry):
    d_=data_riox.rio.clip(geometry,crs=4326,all_touched=True).to_dataframe(name='value').reset_index().dropna()
    d_n=d_[['x','y']].drop_duplicates()
    for i in d_.band.unique():
        r_=d_.loc[d_.band==i,['x','y','value']]
        r_.columns=['x','y',data_riox.long_name[i-1]]
        d_n=r_.merge(d_n,how='left')
    return d_n

def iterate_geom_(data_riox,train_vector):
    r_=pd.DataFrame()
    range_max=train_vector.shape[0]
    for i in tqdm_notebook(range(range_max),desc='Processing clipping raster with vector'):
        try:
            r_t=preprocessing_tif_vector(data_riox,train_vector.geometry[i])
            r_t['TRAIN_CLASS']=train_vector.id[i]
            r_=r_t.append(r_)
        except:
            pass
    return r_

def outlier_dbscan(data):
    columns=['wet_mean',
       'green_mean', 'bright_mean', 'ARVI_mean', 'SAVI_mean', 'NDBI_mean',
       'mNDWI_mean', 'NDWI_mean', 'mNDVI_mean', 'NDVI_mean', 'wet_p50',
       'green_p50', 'bright_p50', 'ARVI_p50', 'SAVI_p50', 'NDBI_p50',
       'mNDWI_p50', 'NDWI_p50', 'mNDVI_p50', 'NDVI_p50', 'S2_B12mean',
       'S2_B11mean', 'S2_B8mean', 'S2_B4mean', 'S2_B3mean', 'S2_B2mean',
       'S2_B12med', 'S2_B11med', 'S2_B8med', 'S2_B4med', 'S2_B3med',
       'S2_B2med']
    t_c=data.TRAIN_CLASS.unique()
    for i in tqdm_notebook(range(len(t_c)),desc='Processing Clustering Outlier data'):
        cl_data=data.loc[data.TRAIN_CLASS==t_c[i],columns].dropna()
        st_sc=Normalizer()
        model_ = DBSCAN(eps =.05,min_samples=10).fit(st_sc.fit_transform(cl_data))
        cl_data['label']=model_.labels_
        data.loc[cl_data.index,'OUTLIER']=cl_data.label
    data['OUTLIER']=data.OUTLIER.apply(lambda y: 0 if y>=0 else -1)
    data_outlier=data.loc[data.OUTLIER<0,['x','TRAIN_CLASS']].groupby('TRAIN_CLASS').agg('count').rename(columns={'x':'COUNT_OUTLIER'}).reset_index()
    fig = px.bar(data_outlier, x="TRAIN_CLASS", y="COUNT_OUTLIER", title="OUTLIER")
    fig.show()
    return data


def rfe_cat(train_x,train_y,valid_x,valid_y,min_):
    train_pool = Pool(train_x, train_y, cat_features=[0])
    valid_pool = Pool(valid_x, valid_y, cat_features=[0])
    f1_score_=[]
    num_feature=[]
    feature_name=[]
    print('Start Recursive Feature Elimination')
    for i in tqdm_notebook(range(min_,36),desc='Iterating Feature Elimination'):
        model=CatBoostClassifier(iterations=50, random_seed=1234,used_ram_limit='10gb')
        summary = model.select_features(
            train_pool,
            eval_set=valid_pool,
            features_for_select='0-34',
            num_features_to_select=i,
            steps=2,
            algorithm=EFeaturesSelectionAlgorithm.RecursiveByShapValues,
            shap_calc_type=EShapCalcType.Regular,
            train_final_model=True,
            logging_level='Silent',
        )
        f1_=f1_score(valid_y,model.predict(valid_pool).tolist(),average='micro')
        f1_score_.append(f1_)
        num_feature.append(i)
        feature_name.append(summary['selected_features_names'])
    print('Best F-1 score: ', max(f1_score_))
    indices=f1_score_.index(max(f1_score_))
    print('Best Number feature: ', num_feature[indices])
    print('Selected of Feature names: \n', feature_name[indices])
    return feature_name[indices]

def tune_param(objective):
    study = optuna.create_study(direction="maximize")
    study.optimize(objective, n_trials=100)
    params = study.best_params
    best_score = study.best_value
    print(f"Best F1_score: {best_s/core}\n")
    print(f"Optimized parameters: {params}\n")
    return params
    
    
def tune_study(train_x,train_y,valid_x,valid_y,feature_select):
    globals()['train_x_1']=train_x
    globals()['train_y_1']=train_y
    globals()['valid_x_1']=valid_x
    globals()['valid_y_1']=valid_y
    globals()['feature_select_1']=feature_select

    study = optuna.create_study(direction="maximize",pruner=optuna.pruners.MedianPruner())
    study.optimize(catb_obj, n_trials=10)  
    best_params = study.best_params
    best_score = study.best_value
    print(f"Best score: {best_score}\n")
    print(f"Optimized parameters: {best_params}\n")
    return best_params
    
def catb_obj(trial):
    param = {
        "objective": "MultiClass",
        "used_ram_limit": "10gb",
        "iterations":50,
        "random_seed":1234,
        "learning_rate":trial.suggest_float("learning_rate",0.001,0.1),
        "depth": trial.suggest_int("depth", 6, 16),
        "l2_leaf_reg":trial.suggest_float("l2_leaf_reg",1e-5,1e-3),
        "colsample_bylevel": trial.suggest_float("colsample_bylevel", 0.01, 0.1),
        "boosting_type": trial.suggest_categorical("boosting_type", ["Plain"]),
        "auto_class_weights":trial.suggest_categorical("auto_class_weights",["Balanced","SqrtBalanced"]),
        "bootstrap_type": trial.suggest_categorical(
            "bootstrap_type", ["Bayesian", "Bernoulli", "MVS"]
        )
    }

    if param["bootstrap_type"] == "Bayesian":
        param["bagging_temperature"] = trial.suggest_float("bagging_temperature", 0, 10)
    elif param["bootstrap_type"] == "Bernoulli":
        param["subsample"] = trial.suggest_float("subsample", 0.1, 1)

    cat_= CatBoostClassifier(**param)
    sk_fold=StratifiedKFold(n_splits=3,shuffle=False)
    f1_=[]
    for train_index, test_index in sk_fold.split(train_x_1, train_y_1):
        if 'PODES_landform' in feature_select_1:
            pool_train=Pool(train_x_1[feature_select_1].iloc[train_index,],train_y_1.iloc[train_index,],cat_features=[0])
            pool_valid=Pool(valid_x_1[feature_select_1],valid_y_1,cat_features=[0])
        else:
            pool_train=Pool(train_x_1.iloc[train_index,][feature_select_1],train_y_1.loc[train_index,])
            pool_valid=Pool(valid_x_1[feature_select_1],valid_y_1)
        cat_.fit(pool_train, eval_set=pool_valid, verbose=0, early_stopping_rounds=100)
        preds = cat_.predict(pool_valid)
        pred_labels = np.rint(preds)
        f1_.append(f1_score(valid_y_1, pred_labels,average='micro'))
    return np.mean(f1_)

def build_model(dict_model,train_x,valid_x,train_y,valid_y,sk_fold):
    cat_= CatBoostClassifier(**dict_model['params'])
    f1_=[]
    for train_index, test_index in sk_fold.split(train_x, train_y):
        if 'PODES_landform' in dict_model['features']:
            pool_train=Pool(train_x.iloc[train_index,],train_y.iloc[train_index,],cat_features=[0])
            pool_valid=Pool(valid_x,valid_y,cat_features=[0])
        else:
            pool_train=Pool(train_x.iloc[train_index,],train_y.iloc[train_index,])
            pool_valid=Pool(valid_x,valid_y)
        cat_.fit(pool_train, eval_set=pool_valid, verbose=0, early_stopping_rounds=100,plot=True)
        preds = cat_.predict(pool_valid)
        pred_labels = np.rint(preds)
        f1_.append(f1_score(valid_y, pred_labels,average='micro'))
    return cat_,f1_

def predict_train(X,y,model,features):
    if 'PODES_landform' in features:
        pool_data=Pool(X,cat_features=[0])
    else:
        pool_data=Pool(X)
    y_predict=model.predict(pool_data)
    report_class=classification_report(y,y_predict)
    print('CLASSIFICATION SUMMARY:')
    print('ACCURACY :',accuracy_score(y,y_predict))
    print('F1-micro score :',f1_score(y,y_predict,average='micro'))
    print('F1-macro score :',f1_score(y,y_predict,average='macro'))
    print(report_class)
    y_=y.copy()
    y_['PREDICTION']=y_predict
    y_['TRUE']=y_.TRAIN_CLASS
    cf_matrix=y_.pivot_table(columns='PREDICTION',index='TRUE',values='TRAIN_CLASS',aggfunc='count').fillna(0)
    plt.figure(figsize=(20, 10))
    sns.heatmap(cf_matrix/np.sum(cf_matrix), annot=True, 
            fmt='.2%', cmap='Blues')

def predict_raster(matrixData,model,dict_model):
    if 'PODES_landform' in dict_model['features']:
        pool_data=Pool(matrixData[dict_model['features']],cat_features=[0])
    else:
        pool_data=Pool(matrixData[dict_model['features']],cat_features=[0])
    matrixData['PREDICT_CLASS']=model.predict(pool_data)
    return matrixData[['x','y','PREDICT_CLASS']]

def check_file_exists(file_name):
    import os.path
    file_path = 'kerangka-spasial/classification_result/' + file_name
    return os.path.exists(file_path)

def create_prediction(list_image, model, dict_model):
    range_max=range(0,len(list_image))
    for i in tqdm_notebook(range_max,desc='Processing Prediction Raster'):
        name_image=list_image[i]
        tif_name='classified_'+name_image[8:15]+'.tif'
        file_name = 'classification_result/'+tif_name
        print('Running prediction for ',name_image)
        if not check_file_exists(tif_name):
            data_riox=rio.open_rasterio('kerangka-spasial/citra-satelit/'+name_image)
            min_x=min(data_riox.x).values+0
            max_x=max(data_riox.x).values+0
            min_y=min(data_riox.y).values+0
            max_y=max(data_riox.y).values+0
            x_=list(np.linspace(min_x,max_x,10))
            y_=list(np.linspace(min_y,max_y,10))
            rangex=range(0,len(x_)-1)
            rangey=range(0,len(y_)-1)
            l_matrixData=[]
            iterate=1
            chunk_size = len(rangex)*len(rangey)
            with tqdm(total=chunk_size,desc='Predicting and Rasterize chunked raster..') as pbar:
                for i in rangex:
                    for j in rangey:
                        t_data_riox=data_riox.rio.clip_box(minx=x_[i],miny=y_[j],maxx=x_[i+1],maxy=y_[j+1])
                        temp_=t_data_riox.to_dataframe(name='value').reset_index().fillna(0)
                        matrixData=temp_[['x','y']].drop_duplicates()
                        for h in temp_.band.unique():
                            t_=temp_.loc[temp_.band==h,['x','y','value']].rename(columns={'value':data_riox.long_name[h-1]})
                            matrixData=matrixData.merge(t_,how='left')
                        dict_data_landform={0:'others',11:'others',12:'others',13:'others',14:'others',15:'others',21: 'u-slope',22:'u-slope',
                           31:'l-slope',32:'l-slope',41:'valley',42:'valley',24:'flat',34:'flat',33:'l-slope',23:'u-slope'}
                        matrixData['PODES_landform']=matrixData['alos_landform'].apply(lambda y: dict_data_landform[y])
                        matrixData=predict_raster(matrixData,model,dict_model)
                        matrixData=gpd.GeoDataFrame(matrixData, crs='epsg:4326',
                             geometry=gpd.points_from_xy(matrixData.x,matrixData.y)).to_crs('epsg:3857')
                        matrixData['PREDICT_CLASS']=matrixData['PREDICT_CLASS'].astype(np.int8)
                        geo_grid = make_geocube(
                            vector_data=matrixData,
                            measurements=['PREDICT_CLASS'],
                            resolution=(-10, 10),
                            rasterize_function=rasterize_points_griddata)
                        geo_grid.PREDICT_CLASS.rio.write_nodata(0,inplace=True)
                        tif_name='classified_'+name_image[8:15]+'_'+str(iterate)+'.tif'
                        iterate=iterate+1
                        geo_grid.rio.to_raster('temp_raster/'+tif_name)
                        pbar.update(1)

            ls_file=os.listdir('temp_raster')
            src_files_to_mosaic = []
            tif_name = 'kerangka-spasial/classification_result/classified_'+name_image[8:15]+'.tif'
            for ls in ls_file:
                src = rasterio.open('temp_raster/'+ls)
                src_files_to_mosaic.append(src)
            mosaic, out_trans = merge(src_files_to_mosaic)
            out_meta = src.meta.copy()
            out_meta.update({"driver": "GTiff","height": mosaic.shape[1],"width": mosaic.shape[2],"transform": out_trans})
            with rasterio.open(tif_name, "w", **out_meta) as dest:
                dest.write(mosaic)
                print(tif_name, 'Done')
            dir = 'temp_raster'
            for f in os.listdir(dir):
                os.remove(os.path.join(dir, f))
        else:
            print("Skip for file "+ file_name+". Already classified")
    print('Finish')
