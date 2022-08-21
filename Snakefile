TUMORS = [8]  # [2,3,4,5]
METHODS = ['elastic',]  #['elastic', 'forest', 'linreg' , 'boosting']
FEATURES = ['string', ]  #['single', 'string', 'corum', 'stringhi' , 'corumplus', 'all']

rule all:
    input: expand('out5_descending/tumor{tumor}_{method}_{feature}.p', tumor=TUMORS, method=METHODS, feature=FEATURES)

rule predict:
    output: 'out5_descending/tumor{tumor}_{method}_{feature}.p'
    log: 'out5_descending/tumor{tumor}_{method}_{feature}.log'
    threads: 4
    shell:
        'python -m predict_protein -n {wildcards.tumor} '
        '-m {wildcards.method} -f {wildcards.feature} -t {threads} -o {output} '
        '1>> {log} 2>> {log}'
