#!/usr/bin/env python
# coding: utf-8

# In[71]:


import os
import pandas as pd


# In[72]:


def load_csv_from_current_directory(filename):
    current_files = os.listdir(os.getcwd())
    if filename in current_files:
        return pd.read_csv(filename, comment='#')


# In[73]:


str_frags = 'frags_mqc.csv'
str_peaks = 'peaks_mqc.csv'


# In[74]:


df_frags = load_csv_from_current_directory(str_frags)
df_peaks = load_csv_from_current_directory(str_peaks)


# In[75]:


def load_enrichment_csvs():
    current_files = os.listdir(os.getcwd())
    
    enrichment_files = [
        f for f in current_files
        if f.startswith('enrichment') and f.endswith('.csv')]
    
    dataframes = [pd.read_csv(f) for f in enrichment_files]
    dataframes_concat = pd.concat(dataframes, ignore_index=True)
    
    return dataframes_concat


# In[76]:


df_enrichment = load_enrichment_csvs()


# In[78]:


def join_sample_dataframes(df_frags,df_peaks,df_enrichment):
    merged = pd.merge(df_frags, df_peaks, on='SampleName', how='inner')
    merged = pd.merge(merged, df_enrichment, on='SampleName', how='outer')
    return merged


# In[79]:


dfJoin = join_sample_dataframes(df_frags,df_peaks,df_enrichment)


# In[80]:


dfJoin = dfJoin.drop(columns=['on_bp', 'off_bp','on_reads','off_reads'], errors='ignore')


# In[81]:


dfJoin = dfJoin.rename(columns={
            'SampleName': 'Sample',
            'Fragments': 'TotalFragments',
            'Peaks': 'TotalPeaks',
            'mark':'Enrichment_Mark',
            'enrichment':'Enrichment_Score'})


# In[82]:


dfJoin = dfJoin.where(pd.notnull(dfJoin),'')


# In[84]:


filename = 'QualityMetrics.csv'


# In[86]:


dfJoin.to_csv(filename, index=False, encoding='utf-8')


# In[ ]:




