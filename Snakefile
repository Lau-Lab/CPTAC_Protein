TUMORS = [2,3,4,5,6,7,8]  # [2,3,4,5]
METHODS = ['linreg', ]  #['elastic', 'forest', 'linreg' , 'boosting']
FEATURES = ['single', 'corum', 'string', 'stringhi', 'all']  #['single', 'string', 'corum', 'stringhi' , 'corumplus', 'all']

rule all:
    input: expand('out3/tumor{tumor}_{method}_{feature}.p', tumor=TUMORS, method=METHODS, feature=FEATURES)

rule predict:
    output: 'out3/tumor{tumor}_{method}_{feature}.p'
    log: 'out3/tumor{tumor}_{method}_{feature}.log'
    threads: 4
    shell:
        'python -m predict_protein -n {wildcards.tumor} '
        '-m {wildcards.method} -f {wildcards.feature} -t {threads} -o {output} '
        '1>> {log} 2>> {log}'
