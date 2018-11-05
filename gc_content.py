from itertools import repeat
import numpy as np

def gc_content(alleles, values):  # .snp data, .geno data
    df = list(repeat(0, 1150639)) # Autosomal SNP count
    gnc = {"G", "C"}
    for al, ind in zip(alleles, values):
        dff = []
	for snpval in ind:
            if al[0] in gnc and al[1] in gnc:
                if snpval == "2":
                    dff.append(2)
                elif snpval == "1":
                    dff.append(2)
                elif snpval == "0":
                    dff.append(2)
                elif snpval == "9":
                    dff.append(0)
            elif al[0] in gnc and al[1] not in gnc:
                if snpval == "2":
                    dff.append(2)
                elif snpval == "1":
                    dff.append(1)
                elif snpval == "0":
                    dff.append(0)
                elif snpval == "9":
                    dff.append(0)
            elif al[0] not in gnc and al[1] in gnc:
                if snpval == "1":
                    dff.append(1)
                elif snpval == "0":
                    dff.append(2)
                elif snpval == "2":
                    dff.append(0)
                elif snpval == "9":
                    dff.append(0)
            else:
                dff.append(0)
        df = np.add(df,dff)
        alleles.seek(0)
        values.seek(0)
    return df
