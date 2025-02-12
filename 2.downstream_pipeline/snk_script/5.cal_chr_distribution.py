import polars as pl
import pandas as pd
import click


@click.command()
@click.option('-d','--workdir', help='The workdir to read and write files.')

def main(workdir):

   df = pd.DataFrame(columns = ['chr','start','end','treatment','condition'])
   df = pl.from_pandas(df)
   df = df.with_columns(
         pl.col('start').cast(pl.Int64),
         pl.col('end').cast(pl.Int64)
         )

   with open(f'{workdir}/stopsite_bedlist.txt','r') as f:
      for line in f:
         try:
            aline = line.strip()
            attri = aline[15:].split('_')
            trt, con = attri[0], attri[1]
            dff = pl.read_csv(aline, has_header=True, separator='\t', new_columns=['chr','start','end'])
            dff = dff.with_columns(
                pl.lit(trt).alias('treatment'),
                pl.lit(con).alias('condition')
            )
            df = df.vstack(dff)
         except:
            continue

   #df = pl.read_csv(inputbed,has_header=True, separator='\t', new_columns=['chr','start','end'])
   df = df.with_columns(
      (pl.col('end')-pl.col('start')).alias('length')
   )
   chr_length = pl.DataFrame({
      'chr1': 248956422,
      'chr2': 242193529,
      'chr3': 198295559,
      'chr4': 190214555,
      'chr5': 181538259,
      'chr6': 170805979,
      'chr7': 159345973,
      'chr8': 145138636,
      'chr9': 138394717,
      'chr10': 133797422,
      'chr11': 135086622,
      'chr12': 133275309,
      'chr13': 114364328,
      'chr14': 107043718,
      'chr15': 101991189,
      'chr16': 90338345,
      'chr17': 83257441,
      'chr18': 80373285,
      'chr19': 58617616,
      'chr20': 64444167,
      'chr21': 46709983,
      'chr22': 50818468,
      'chrX' : 156040895,
      'chrY': 57227415,
      'chrM': 16569
   })
   chr_length = chr_length.melt().rename({'variable':'chr','value':'length'})
   chr_length = chr_length.with_columns(
      (pl.col('length') / pl.col('length').sum()).alias('ratio') #ratio: chr占整个基因组的比例
   )


   # 原假设，所有区域在所有染色体上平均分布（按照长度加权平均）：计算expectation 计算所有region的长度，并且按照染色体长度比例进行分配
   mean_dmrs = df.group_by(pl.col('treatment'), pl.col('condition')).agg([
      pl.col('length').mean().alias('length'),
      pl.col('chr').count().alias('number')
   ])

   mean_dmrs = mean_dmrs.with_columns(
      (pl.col('length')*pl.col('number')).alias('total')
   )
   exp = chr_length.select([
      pl.col('chr'),
      pl.col('ratio')
   ])

   expall = pd.DataFrame(columns = ['chr','ratio','exp','treatment','condition'])
   expall = pl.from_pandas(expall)
   expall = expall.with_columns(
      pl.col('ratio').cast(pl.Float64),
      pl.col('exp').cast(pl.Float64)
   )

   for i in range(mean_dmrs.shape[0]):
      atreatment = mean_dmrs[i]['treatment'][0]
      #aantibody = mean_dmrs[i]['antibody'][0]
      acondition = mean_dmrs[i]['condition'][0]
      total = mean_dmrs[i]['total'][0]
      exptmp = exp.with_columns(
         (pl.col('ratio')*total).alias('exp'),
         pl.lit(atreatment).alias('treatment'),
         pl.lit(acondition).alias('condition')
      )
         #pl.lit(aantibody).alias('antibody'),
      expall = expall.vstack(exptmp)

   obs = df.group_by(pl.col('chr'), pl.col('treatment'), pl.col('condition')).agg([
      pl.col('length').sum().alias('obs')
   ])
   #obs = obs.pivot(index='chr',values='length')

   res = obs.join(expall,on=['chr','treatment','condition'])
   res.with_columns(
      (pl.col('obs')/pl.col('exp')).alias('enrichment') #对数处理在作图脚本中
   ).select([
      pl.col('chr'),
      pl.col('treatment'),
      pl.col('condition'),
      pl.col('enrichment')
   ]).write_csv(f'{workdir}/stopsite_enrichment.txt',include_header=True, separator='\t')
      #pl.col('antibody'),

if  __name__ == '__main__':
   main()
