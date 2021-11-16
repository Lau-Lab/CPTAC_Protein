TUMORS = [8]
METHODS = ['elastic', 'forest', 'linreg', 'boosting']
FEATURES = ['single', 'string', 'corum', 'all']

rule all:
    input: expand('out2/tumor{tumor}_{method}_{feature}.p', tumor=TUMORS, method=METHODS, feature=FEATURES)

rule predict:
    output: 'out2/tumor{tumor}_{method}_{feature}.p'
    log: 'out2/tumor{tumor}_{method}_{feature}.log'
    threads: 4
    shell:
        'python -m predict_protein -n {wildcards.tumor} '
        '-m {wildcards.method} -f {wildcards.feature} -t {threads} -o {output} '
        '1>> {log} 2>> {log}'
