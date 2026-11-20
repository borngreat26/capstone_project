"""Fetch only the selected original GWASs and Week 7 metadata/chain/peaks.
Existing files are preserved. Interrupted downloads use a .part file, never a valid input filename.
Full original summary statistics stay local and are not redistributed.
"""
import gzip,sys,urllib.request
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.utils import REPO_ROOT

BASE='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/'
FILES={
 'data/gwas_original/GCST006803.txt.gz':'https://walters.psycm.cf.ac.uk/clozuk_pgc2.meta.sumstats.txt.gz',
 'data/gwas_original/GCST006803.md5':'https://walters.psycm.cf.ac.uk/clozuk_pgc2.meta.sumstats.txt.gz.md5',
 'data/gwas_original/GCST006803.source.html':'https://walters.psycm.cf.ac.uk/',
 'data/gwas_original/GCST005569.Build37.f.tsv.gz':BASE+'GCST005001-GCST006000/GCST005569/harmonised/23143596-GCST005569-EFO_0000685-Build37.f.tsv.gz',
 'data/gwas_original/GCST007517.Build37.f.tsv.gz':BASE+'GCST007001-GCST008000/GCST007517/harmonised/29632382-GCST007517-EFO_0001360-build37.f.tsv.gz',
 'data/gwas_original/GCST005569.meta.yaml':BASE+'GCST005001-GCST006000/GCST005569/harmonised/23143596-GCST005569-EFO_0000685.h.tsv.gz-meta.yaml',
 'data/gwas_original/GCST007517.meta.yaml':BASE+'GCST007001-GCST008000/GCST007517/harmonised/29632382-GCST007517-EFO_0001360.h.tsv.gz-meta.yaml',
 'data/gwas_original/GCST005569.readme.txt':BASE+'GCST005001-GCST006000/GCST005569/harmonised/readme.txt',
 'data/encode/ENCFF421XIL.metadata.json':'https://www.encodeproject.org/files/ENCFF421XIL/?format=json',
 'data/encode/ENCSR603LVR.metadata.json':'https://www.encodeproject.org/experiments/ENCSR603LVR/?format=json',
 'data/encode/hg38ToHg19.over.chain.gz':'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz',
}

def fetch(relative,url):
    target=REPO_ROOT/relative
    if target.exists():return target
    target.parent.mkdir(parents=True,exist_ok=True);partial=target.with_name(target.name+'.part')
    request=urllib.request.Request(url,headers={'User-Agent':'capstone-public-data/1.0','Accept':'application/json' if '.json' in relative else '*/*'})
    with urllib.request.urlopen(request,timeout=120) as response,partial.open('wb') as handle:
        import shutil
        shutil.copyfileobj(response,handle)
    if target.suffix=='.gz':
        with gzip.open(partial,'rb') as handle:
            while handle.read(1024*1024):pass
    partial.replace(target);print('Downloaded',relative)
    return target

def main():
    for path,url in FILES.items():fetch(path,url)
    target=REPO_ROOT/'data/encode/ENCFF421XIL.bed'
    if not target.exists():
        archive=fetch('data/encode/ENCFF421XIL.bed.gz','https://www.encodeproject.org/files/ENCFF421XIL/@@download/ENCFF421XIL.bed.gz')
        with gzip.open(archive,'rb') as source,target.open('wb') as dest:
            import shutil
            shutil.copyfileobj(source,dest)
    print('Existing inputs retained. run_week7.py verifies the ENCODE content checksum.')
if __name__=='__main__':main()
