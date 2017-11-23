import re

RAW_INTERVAR_CLASS_BENIGN = "Benign"
RAW_INTERVAR_CLASS_LIKELY_BENIGN = "Likelybenign"
RAW_INTERVAR_CLASS_UNCERTAIN_SIGNIFICANCE = "UncertainSignificance"
RAW_INTERVAR_CLASS_LIKELY_PATHOGENIC = "Likelypathogenic"
RAW_INTERVAR_CLASS_PATHOGENIC = "Pathogenic"

INTERVAR_CLASS_BENIGN = "Benign"
INTERVAR_CLASS_LIKELY_BENIGN = "Likely Benign"
INTERVAR_CLASS_UNCERTAIN_SIGNIFICANCE = "Uncertain Significance"
INTERVAR_CLASS_LIKELY_PATHOGENIC = "Likely Pathogenic"
INTERVAR_CLASS_PATHOGENIC = "Pathogenic"

CLASSIFICATION_PATTERN = re.compile(r'''InterVar:(?P<acmg_class>.+?);''')
EVIDENCE_PATTERN = re.compile(r'''(?P<var_name>[a-zA-Z0-9]*?)=(?P<value>(?:[0-9]+?|\[[0-9;]*?\]))''')


def parse_intervar_class(raw_intervar):
    class_match = CLASSIFICATION_PATTERN.match(raw_intervar)
    if class_match is not None:
        intervar_class = class_match.group('acmg_class')
        if intervar_class == RAW_INTERVAR_CLASS_BENIGN:
            return INTERVAR_CLASS_BENIGN
        if intervar_class == RAW_INTERVAR_CLASS_LIKELY_BENIGN:
            return INTERVAR_CLASS_LIKELY_BENIGN
        if intervar_class == RAW_INTERVAR_CLASS_UNCERTAIN_SIGNIFICANCE:
            return INTERVAR_CLASS_UNCERTAIN_SIGNIFICANCE
        if intervar_class == RAW_INTERVAR_CLASS_LIKELY_PATHOGENIC:
            return INTERVAR_CLASS_LIKELY_PATHOGENIC
        if intervar_class == RAW_INTERVAR_CLASS_PATHOGENIC:
            return INTERVAR_CLASS_PATHOGENIC
    return ""

def evidence2str(raw_evidence):
    evidence_list = []
    for item in raw_evidence:
        var_name = item[0]
        value = eval(item[1].replace(';',','))
        if type(value) is int and value == 1:
            evidence_list.append(var_name)
        elif type(value) is list:
            for value_idx in xrange(len(value)):
                var_name_val = value[value_idx]
                if var_name_val == 1:
                    evidence_list.append(var_name+str(value_idx+1))
    return ", ".join(evidence_list)

def parse_intervar_evidence(raw_intervar):
    class_match = CLASSIFICATION_PATTERN.match(raw_intervar)
    evidence_matchs = EVIDENCE_PATTERN.findall(raw_intervar, re.DOTALL)
    return evidence2str(evidence_matchs)
