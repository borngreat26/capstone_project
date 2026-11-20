"""BED half-open interval operations; SNP positions supplied as one-based coordinates."""
import numpy as np
import pandas as pd

PEAK_COLUMNS=['Chrom','Start','End','Name','Score','Strand','Signal','NegLog10P','NegLog10Q','SummitOffset']

def read_peaks(path):
    peaks=pd.read_csv(path,sep='\t',header=None,comment='#')
    if peaks.shape[1] not in (6,10): raise ValueError('Expected BED6 or BED6+4 narrowPeak')
    peaks.columns=PEAK_COLUMNS[:peaks.shape[1]]
    peaks['Chrom']=peaks.Chrom.astype(str).str.removeprefix('chr')
    for c in ['Start','End']: peaks[c]=pd.to_numeric(peaks[c],errors='raise').astype('int64')
    if not ((peaks.Start>=0)&(peaks.End>peaks.Start)).all(): raise ValueError('Malformed peak interval')
    # Exact duplicated records cannot count twice; distinct original annotations remain intact.
    return peaks.drop_duplicates().reset_index(drop=True)

def interval_counts(peaks,chrom,starts,ends):
    """Count overlapping peaks even when peaks nest or overlap one another."""
    p=peaks.loc[peaks.Chrom.eq(str(chrom).removeprefix('chr'))]
    return np.searchsorted(np.sort(p.Start),ends,side='left')-np.searchsorted(np.sort(p.End),starts,side='right')

def snp_counts(peaks,chrom,positions):
    pos=np.asarray(positions,dtype=np.int64)
    if (pos<1).any(): raise ValueError('SNP positions must be one-based positive integers')
    return interval_counts(peaks,chrom,pos-1,pos)

def matching_peaks(peaks,chrom,start,end):
    return peaks.loc[peaks.Chrom.eq(str(chrom).removeprefix('chr')) & peaks.Start.lt(end) & peaks.End.gt(start)]
