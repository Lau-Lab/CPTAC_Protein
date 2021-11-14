TUMORS = [2,4,8]
METHODS = ['elastic', 'forest', 'linreg']
FEATURES = ['single', 'string', 'corum'] #, 'all']

rule all:
    input: expand('out/tumor{tumor}_{method}_{feature}.p', tumor=TUMORS, method=METHODS, feature=FEATURES)

rule predict:
    output: 'out/tumor{tumor}_{method}_{feature}.p'
    log: 'out/tumor{tumor}_{method}_{feature}.log'
    shell:
        'python -m predict_protein -n {wildcards.tumor} '
        '-m {wildcards.method} -f {wildcards.feature} -o {output} '
        '1>> {log} 2>> {log}'
