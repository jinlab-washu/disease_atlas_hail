import pprint
import argparse
from datetime import datetime
import os
import pyspark
from pyspark.sql import SparkSession

from annotate_frequencies import *
from generate_split_alleles import *
from prepare_ht_export import *
from prepare_ht_for_es import *
from export_ht_to_es import *

def get_adj_expr(
        gt_expr: hl.expr.CallExpression,
        gq_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
        dp_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
        ad_expr: hl.expr.ArrayNumericExpression,
        adj_gq: int = 20,
        adj_dp: int = 10,
        adj_ab: float = 0.2,
        haploid_adj_dp: int = 10
) -> hl.expr.BooleanExpression:
    """
    Gets adj genotype annotation.
    Defaults correspond to gnomAD values.
    """
    return (
            (gq_expr >= adj_gq) &
            hl.cond(
                gt_expr.is_haploid(),
                dp_expr >= haploid_adj_dp,
                dp_expr >= adj_dp
            ) &
            (
                hl.case()
                .when(~gt_expr.is_het(), True)
                .when(gt_expr.is_het_ref(), ad_expr[1] / dp_expr >= adj_ab)
                .default((ad_expr[0] / dp_expr >= adj_ab ) & (ad_expr[1] / dp_expr >= adj_ab ))
            )
    )

def annotate_adj(
        mt: hl.MatrixTable,
        adj_gq: int = 20,
        adj_dp: int = 10,
        adj_ab: float = 0.2,
        haploid_adj_dp: int = 10
) -> hl.MatrixTable:
    """
    Annotate genotypes with adj criteria (assumes diploid)
    Defaults correspond to gnomAD values.
    """
    return mt.annotate_entries(adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD, adj_gq, adj_dp, adj_ab, haploid_adj_dp))

def qual_hist_expr(
    gt_expr: Optional[hl.expr.CallExpression] = None,
    gq_expr: Optional[hl.expr.NumericExpression] = None,
    dp_expr: Optional[hl.expr.NumericExpression] = None,
    ad_expr: Optional[hl.expr.ArrayNumericExpression] = None,
    adj_expr: Optional[hl.expr.BooleanExpression] = None,
) -> hl.expr.StructExpression:
    """
    Return a struct expression with genotype quality histograms based on the arguments given (dp, gq, ad).
    .. note::
        - If `gt_expr` is provided, will return histograms for non-reference samples only as well as all samples.
        - `gt_expr` is required for the allele-balance histogram, as it is only computed on het samples.
        - If `adj_expr` is provided, additional histograms are computed using only adj samples.
    :param gt_expr: Entry expression containing genotype
    :param gq_expr: Entry expression containing genotype quality
    :param dp_expr: Entry expression containing depth
    :param ad_expr: Entry expression containing allelic depth (bi-allelic here)
    :param adj_expr: Entry expression containing adj (high quality) genotype status
    :return: Genotype quality histograms expression
    """
    qual_hists = {}
    if gq_expr is not None:
        qual_hists["gq_hist"] = hl.agg.hist(gq_expr, 0, 100, 20)
    if dp_expr is not None:
        qual_hists["dp_hist"] = hl.agg.hist(dp_expr, 0, 100, 20)

    if gt_expr is not None:
        qual_hists = {
            **{
                f"{qual_hist_name}_all": qual_hist_expr
                for qual_hist_name, qual_hist_expr in qual_hists.items()
            },
            **{
                f"{qual_hist_name}_alt": hl.agg.filter(
                    gt_expr.is_non_ref(), qual_hist_expr
                )
                for qual_hist_name, qual_hist_expr in qual_hists.items()
            },
        }
        if ad_expr is not None:
            qual_hists["ab_hist_alt"] = hl.agg.filter(
                gt_expr.is_het(), hl.agg.hist(ad_expr[1] / hl.sum(ad_expr), 0, 1, 20)
            )

    else:
        qual_hists = {
            f"{qual_hist_name}_all": qual_hist_expr
            for qual_hist_name, qual_hist_expr in qual_hists.items()
        }

    if adj_expr is not None:
        qual_hists.update(
            {
                f"{qual_hist_name}_adj": hl.agg.filter(adj_expr, qual_hist_expr)
                for qual_hist_name, qual_hist_expr in qual_hists.items()
            }
        )

    return hl.struct(**qual_hists)

def subset_samples_and_variants(mt, subset_path):
    """
    Subset the MatrixTable to the provided list of samples and to variants present in those samples
    :param mt: MatrixTable from VCF
    :param subset_path: Path to a file with a single column 's'
    :return: MatrixTable subsetted to list of samples
    """
    subset_ht = hl.import_table(subset_path, key='s')
    subset_count = subset_ht.count()
    anti_join_ht = subset_ht.anti_join(mt.cols())
    anti_join_ht_count = anti_join_ht.count()

    if anti_join_ht_count != 0:
        missing_samples = anti_join_ht.s.collect()
        message = f'Only {subset_count - anti_join_ht_count} out of {subset_count} ' \
                  f'subsetting-table IDs matched IDs in the variant callset.\n' \
                  f'IDs that aren\'t in the callset: {missing_samples}\n' \
                  f'All callset sample IDs:{mt.s.collect()}'

        '''
        if (subset_count > anti_join_ht_count) and self.ignore_missing_samples:
            logger.warning(message)
        else:
            raise MatrixTableSampleSetError(message, missing_samples)
        '''


    mt = mt.semi_join_cols(subset_ht)
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    '''
    logger.info(f'Finished subsetting samples. Kept {subset_count} '
                f'out of {mt.count()} samples in vds')
    '''

    return mt

## Hail initiate
def hail_initialize(jobname: str = None, mem: int = 32):
    """
    To initialize the hail environment
    """
    today = datetime.now().strftime("%Y-%m-%d")
    print("Program Starting Date: {}".format(today))

    #HAIL_HOME = subprocess.getoutput("pip3 show hail | grep Location | awk -F' ' '{print $2 \"/hail\"}'")
    HAIL_HOME = "/usr/local/lib/python3.8/dist-packages/hail"
    os.environ["HAIL_DIR"] = f"{HAIL_HOME}/backend"
    # os.environ["HAIL_DIR"] = "/opt/conda/lib/python3.7/site-packages/hail/backend"
    # os.environ["JAVA_HOME"] = "/usr/lib/jvm/java-11-openjdk-amd64"
    os.environ["JAVA_HOME"] = "/usr/lib/jvm/java-8-openjdk-amd64/jre"

    threads = int(os.environ['LSB_MAX_NUM_PROCESSORS'])
    print("Run with {} CPUs.".format(threads))

    # Ensure hail-all-spark.jar exists
    hail_jars = f"{HAIL_HOME}/backend/hail-all-spark.jar"
    if not os.path.exists(hail_jars):
        raise FileNotFoundError(f"Cannot find hail-all-spark.jar at {hail_jars}")
    print(f"Found hail-all-spark.jar at: {hail_jars}")

    # Convert 'mem' to integer if it's not already
    mem = int(mem)

    # Calculate memory settings
    spare_memory = 6 #int(mem * 0.25) + 1
    driver_memory = int(mem * 0.25)  # Use 25% of total memory for driver
    executor_memory = mem - driver_memory - spare_memory  # Subtract driver memory and overhead for executors

    # Spark Configuration
    conf = pyspark.SparkConf().setAll([
        ('spark.master', 'local[{}]'.format(threads)),
        ('spark.app.name', 'Hail'),
        ('spark.jars', f'file://{hail_jars}'),
        ('spark.driver.extraClassPath', f'file://{hail_jars}'),
        ('spark.executor.extraClassPath', f'file://{hail_jars}'),   #'./hail-all-spark.jar'
        ('spark.serializer', 'org.apache.spark.serializer.KryoSerializer'),
        ('spark.kryo.registrator', 'is.hail.kryo.HailKryoRegistrator'),
        ### https://discuss.hail.is/t/turning-run-combiner-performance-for-hail-local-mode/2318
        ('spark.driver.memory', '{}g'.format(driver_memory)),
        ('spark.executor.memory', '{}g'.format(executor_memory)),
        ('spark.driver.maxResultSize', '{}g'.format(driver_memory)),
        ### For IOException: No space left on device
        # https://spark.apache.org/docs/2.3.0/configuration.html
        ('spark.local.dir', '/storage1/fs1/jin810/Active/BRDGE/tmp'),
        ### Error: KryoException: Buffer overflow. Available: 0, required: 15069596
        ('spark.kryoserializer.buffer', '64m'),
        ('spark.kryoserializer.buffer.max', '1024m'),
        ('spark.driver.memoryOverhead', '512m'),
        ('spark.executor.memoryOverhead', '1g'),
        ### To output some essential info for memory leakage (Garbage collection logging)
        ('spark.executor.extraJavaOptions', '-XX:+UseParallelGC -XX:+PrintGCDetails -XX:+PrintGCTimeStamps ' +
         '-XX:+PrintGCDateStamps -XX:+UseGCLogFileRotation -XX:NumberOfGCLogFiles=10 -XX:GCLogFileSize=100M ' +
         '-XX:ConcGCThreads={}'.format(threads)),
    ])

    # Initialize Spark Session
    spark = SparkSession.builder.config(conf=conf).getOrCreate()
    sc = spark.sparkContext

    # Log file setup
    #log_base_path = "/storage1/fs1/jin810/Active/testing/yung-chun/Hail_logs"
    #log_file_name = "hail_vds_new_combiner_{}.log".format(jobname if jobname else today)
    #logfile = os.path.join(log_base_path, log_file_name)
    #logfile = create_unique_filename(logfile)

    # Hail Initialization
    #hl.init(default_reference='GRCh38', sc=sc, log=logfile)
    hl.init(sc=sc)
    # hl.init(master='local[16]') -> to use all the CPU

    return spark, sc

def run_pipeline(args):

    spark, sc = hail_initialize("Lektest", mem="30")

    '''
    Import VCF   
    '''

    mt = hl.import_vcf(args.vcf,reference_genome='GRCh38',array_elements_required=False,force_bgz=True)
    mt = mt.repartition(64)


    if args.subset:
        mt = subset_samples_and_variants(mt, args.subset)

    #Split alleles
    mt = generate_split_alleles(mt)

    mt = annotate_adj(mt)

    mt = mt.annotate_rows(raw_qual_hists=qual_hist_expr(mt.GT,mt.GQ,mt.DP,mt.AD,mt.adj))

    #Annotate Population frequencies for now
    meta_ht = hl.import_table(args.meta,delimiter='\t',key='ID')
    ht = annotate_frequencies(mt,meta_ht)

    res_ht = hl.read_table('/storage1/fs1/jin810/Active/References/2023_old_references/hail/Ref_HailFormat/combined_reference_data_grch38-2.0.4.ht')
    dbsnp_ht = hl.read_table('/storage1/fs1/jin810/Active/BRDGE/Reference/dbsnp_b151_grch38_all_20180418.ht')
    #clinvar_ht = hl.read_table('/mnt/home/mlek/ceph/resources/clinvar_20190923.ht')

    ht = ht.annotate(rsid = dbsnp_ht[ht.key].rsid,
        in_silico_predictors=hl.struct(
            splice_ai = res_ht[ht.key].splice_ai.delta_score, 
            primate_ai= res_ht[ht.key].primate_ai.score,
            cadd= res_ht[ht.key].cadd.PHRED,
            revel= res_ht[ht.key].dbnsfp.REVEL_score
        ),

        allele_balance=hl.struct(
            alt_raw=ht.raw_qual_hists.ab_hist_alt.annotate(
                bin_edges=ht.raw_qual_hists.ab_hist_alt.bin_edges.map(lambda n: hl.float(hl.format("%.3f", n)))
            ),
            alt_adj=ht.raw_qual_hists.ab_hist_alt_adj.annotate(
                bin_edges=ht.raw_qual_hists.ab_hist_alt_adj.bin_edges.map(lambda n: hl.float(hl.format("%.3f", n)))
            )
        ),
        genotype_depth=hl.struct(
            all_raw=ht.raw_qual_hists.dp_hist_all,
            alt_raw=ht.raw_qual_hists.dp_hist_alt,
            all_adj=ht.raw_qual_hists.dp_hist_all_adj,
            alt_adj=ht.raw_qual_hists.dp_hist_alt_adj
        ),
        genotype_quality=hl.struct(
            all_raw=ht.raw_qual_hists.gq_hist_all,
            alt_raw=ht.raw_qual_hists.gq_hist_alt,
            all_adj=ht.raw_qual_hists.gq_hist_all_adj,
            alt_adj=ht.raw_qual_hists.gq_hist_alt_adj
        )
    )

    #VEP Annotate the Hail table (ie. sites-only) using GRCh38 configuration file
    ht = hl.vep(ht, 'file:///home/lekm/vep104.json')

    ht = prepare_ht_export(ht)
    ht = prepare_ht_for_es(ht)

    pprint.pprint(ht.describe())
    pprint.pprint(ht.show())


    ht.write(args.out,overwrite=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', '--input', '-i', help='bgzipped VCF file (.vcf.bgz)', required=True)
    parser.add_argument('--meta', '-m', help='Meta file containing sample population and sex', required=True)
    parser.add_argument('--subset', '-s', help='Samples to subset and keep', required=False)
    parser.add_argument('--out', '-o', help='Hail table output file name', required=True)

    args = parser.parse_args()
    run_pipeline(args)
