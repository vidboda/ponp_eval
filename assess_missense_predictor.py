import os
import sys
import re
import time
import yaml
import argparse
import tabix
import math

# script to assess missense prédictors
# based on curated datasets ("truth set") provided as text files
# format:
# hgvs_genomic\tvariant_class
# one file for pathogenic
# another for neutrals
# scores retrieved via a tabix query



def log(level, text):
    localtime = time.asctime( time.localtime(time.time()) )
    if level == 'ERROR':
        sys.exit('[{0}]: {1} - {2}'.format(level, localtime, text))
    print('[{0}]: {1} - {2}'.format(level, localtime, text))


def main():
    parser = argparse.ArgumentParser(
        description='Assess a missense variant predictor',
        usage='python assess_missense_predictor.py -c yaml_conf_file -p predictor_name -t truth_set'
    )
    parser.add_argument('-p', '--predictor', default='', required=True,
                        help='Predictor name (see config file)')
    parser.add_argument('-t', '--truth-set', default='', required=True,
                        help='prefix of truth set files (see config file)')
    parser.add_argument('-c', '--config-file', default='', required=True,
                        help='Config file yaml format. See assess_missense.yaml for examples.')
    
    args = parser.parse_args()
    predictor = args.predictor
    truth_set_name = args.truth_set
    # log('DEBUG', os.path.splitext(args.config_file))
    if os.path.splitext(args.config_file)[1] == '.yaml':
        # get config file
        config = yaml.safe_load(open(args.config_file))
        
        predictors = config['predictors']
        truth_set = config['truth_set']
        # log('DEBUG', truth_set)
        
        if predictor in predictors and \
                truth_set_name in truth_set['pathogenic'] and \
                truth_set_name in truth_set['neutral']:
            log('INFO', 'Predictor: {0}'.format(predictor))
            log('INFO', 'Predictor file path: {0}'.format(predictors[predictor]['file_path']))
            log('INFO', 'Pathogenic truth file: {0}{1}'.format(truth_set['path'], truth_set['pathogenic'][truth_set_name]['full_name']))
            log('INFO', 'Neutral truth file: {0}{1}'.format(truth_set['path'], truth_set['neutral'][truth_set_name]['full_name']))
            
            total_set_pathogenic, total_set_pathogenic_analysed, tp, fp, res_dict_pathogenic = assess_variant('pathogenic', predictor, truth_set_name, predictors, truth_set)
            # log('DEBUG', 'TP: {0} - FP: {1}'.format(tp, fp))
            total_set_neutral, total_set_neutral_analysed, tn, fn, res_dict_neutral = assess_variant('neutral', predictor, truth_set_name, predictors, truth_set)
            # log('DEBUG', 'TN: {0} - FN: {1}'.format(tn, fn))
                
            statistics = compute_statistics(tp, fp, tn, fn)
            # log('DEBUG', res_dict_pathogenic)
            perc_pathogenic_analysed = "%.2f" % ((total_set_pathogenic_analysed / total_set_pathogenic) * 100)
            log('INFO', '% of analysed pathogenic variants: {0} % ({1}/{2})'.format(
                perc_pathogenic_analysed,
                total_set_pathogenic_analysed,
                total_set_pathogenic
            ))
            # log('INFO', 'Total number of pathogenic variants tested: {}'.format(total_set_pathogenic))
            log('INFO', 'Number of pathogenic variants correctly predicted (TP): {}'.format(tp))
            log('INFO', 'Number of pathogenic variants uncorrectly predicted (FP): {}'.format(fp))
            
            perc_neutral_analysed = "%.2f" % ((total_set_neutral_analysed / total_set_neutral) * 100)
            log('INFO', '% of analysed neutral variants: {0} % ({1}/{2})'.format(
                perc_neutral_analysed,
                total_set_neutral_analysed,
                total_set_neutral
            ))
            #log('INFO', 'Total number of neutral variants tested: {}'.format(total_set_neutral))
            log('INFO', 'Number of neutral variants correctly predicted (TN): {}'.format(tn))
            log('INFO', 'Number of neutral variants uncorrectly predicted (FN): {}'.format(fn))
            
            for stat in statistics:
                log('INFO', '{0}: {1}'.format(stat, statistics[stat]))

            # for var in res_dict_neutral:
            #     log('DEBUG', '{0}: {1}_class: {2} - {3}_score: {4}'.format(
                    #     var,
                    #     truth_set_name,
                    #     res_dict_neutral[var]['{}_class'.format(truth_set_name)],
                    #     predictor,
                    #     res_dict_neutral[var]['{}_score'.format(predictor)])
                    # )

                            
def assess_variant(dataset, predictor, truth_set_name, predictors, truth_set):    
    total_set = 0
    total_set_analysed = 0
    true = 0
    false = 0
    res_dict = {}
    with open('{0}{1}'.format(truth_set['path'], truth_set[dataset][truth_set_name]['full_name'])) as variant_set: 
        for variant in variant_set:
            if re.search(r'^[^#]', variant):
                total_set += 1
                match_obj = re.search(r'^chr([\dXY]{1,2}):g\.(\d+)([ATGC])>([ATGC])\s+(.+)$', variant)
                if match_obj:
                    chrom = match_obj.group(1)
                    pos = match_obj.group(2)
                    ref = match_obj.group(3)
                    alt = match_obj.group(4)
                    classif = match_obj.group(5)                        
                    # tabix query
                    tb = tabix.open(predictors[predictor]['file_path'])
                    query = "{0}:{1}-{2}".format(chrom, pos, pos)
                    try:
                        records = tb.querys(query)
                    except Exception as e:
                        # log('DEBUG', 'Tabix failed:{0} for variant {1}'.format(e.args, 'chr{0}:g.{1}{2}>{3}'.format(chrom, pos, ref, alt)))
                        continue
                    records = tb.querys(query)
                    if records:
                        for record in records:
                            # log('DEBUG', 'tabix: {0} - variant: chr{1}:g.{2}{3}>{4}'.format(record, chrom, pos, ref, alt))
                            # log('DEBUG', 'File ref: {0} - ref: {1} - File alt:{2} - alt: {3}'.format(record[predictors[predictor]['ref_index']], ref, record[predictors[predictor]['alt_index']], alt))
                            if record[predictors[predictor]['ref_index']] == ref and \
                                    record[predictors[predictor]['alt_index']] == alt:
                                # treat multiple result as in dbNSFP
                                score = record[predictors[predictor]['score_index']]
                                # log('DEBUG', score)
                                if re.search(r';', score):
                                    # treat multiple result as in dbNSFP
                                    scores = re.split(';', score)
                                    score = max(scores)
                                    # log('DEBUG', '{0} - {1}'.format(scores, score))
                                if str(score) == '.':
                                    break
                                if dataset == 'pathogenic':
                                    if predictors[predictor]['score_sort'] == 'gt':
                                        if float(score) > float(predictors[predictor]['score_threshold']):
                                            # TP
                                            true += 1
                                        else:
                                            # FP
                                            false += 1
                                    else:
                                        if float(score) < float(predictors[predictor]['score_threshold']):
                                            # TP
                                            true += 1
                                        else:
                                            # FP
                                            false += 1
                                elif dataset == 'neutral':
                                    if predictors[predictor]['score_sort'] == 'gt':
                                        if float(score) < float(predictors[predictor]['score_threshold']):
                                            # TP
                                            true += 1
                                        else:
                                            # FP
                                            false += 1
                                    else:
                                        if float(score) > float(predictors[predictor]['score_threshold']):
                                            # TP
                                            true += 1
                                        else:
                                            # FP
                                            false += 1
                                res_dict['chr{0}:g.{1}{2}>{3}'.format(chrom, pos, ref, alt)] = {
                                    '{}_class'.format(truth_set_name): classif,
                                    '{}_score'.format(predictor): record[predictors[predictor]['score_index']]
                                }
                                total_set_analysed += 1
                                break
                else:
                    log('WARNING', 'Wrong format for variant line: {}'.format(variant))
            else:
                log('INFO', variant.rstrip())
    return total_set, total_set_analysed, true, false, res_dict

def compute_statistics(tp=0, fp=0, tn=0, fn=0):
    stats = {}
    if tp > 0 and fp > 0:
        stats['ppv'] =  "%.2f" % (tp / (tp + fp))
    if fn > 0 and tn > 0:
        stats['npv'] = "%.2f" % (tn / (fn + tn))
    if (tp + fn) > 0:
        stats['sensitivity'] = "%.2f" % (tp / (tp + fn))
    if (fp + tn) > 0:
        stats['specificity'] = "%.2f" % (tn / (fp + tn))
    if (tp + fp + fn + tn) > 0:
        stats['accuracy'] = "%.2f" % ((tp + tn) / (tp + fp + fn + tn))
    if ((tp + fn) * (tn + fp) * (tp + fp) * (tn + fn)) > 0:
        stats['mcc'] = "%.2f" % (((tp * tn) - (fp * fn)) / math.sqrt((tp + fn) * (tn + fp) * (tp + fp) * (tn + fn)))
    if 'ppv' in stats and 'sensitivity' in stats:
        stats['f-measure'] = "%.2f" % (2 * ((float(stats['ppv']) * float(stats['sensitivity']))) / ((float(stats['ppv']) + float(stats['sensitivity']))))
    return stats


if __name__ == '__main__':
    main()