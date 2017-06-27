import limix
import pandas as pd
import numpy as np
import math

def do_snp_qc(snp_df, min_call_rate, min_maf, min_hwe_P, min_hmachR2):
   
    #Determine call rate.
    call_rate = 1-snp_df.isnull().sum()/len(snp_df.index)
    selection = call_rate < min_call_rate
    #print(call_rate)
    #print(call_rate < min_call_rate)
    failed_snp_names  = list(snp_df.columns[selection])
    snp_df = snp_df.loc[:,list(snp_df.columns[~selection])]
    if(len(snp_df.columns)==0):
        return snp_df.columns, failed_snp_names
    #Determine MAF.
    genotypeCounter = np.zeros((len(snp_df.columns),3), dtype=np.int)
    for allele_index in [0,1,2]:
        genotypeCounter[:,allele_index] = np.nansum(np.around(snp_df.values)==allele_index,0)

    #Here we make sure that the major allele is temporarly 'coded' as 0 & directly calculate the MAF (based on allele counts and non NA samples)
    mac = np.zeros((len(snp_df.columns)), dtype=np.int)
    gc = np.nansum(snp_df.values,0)
    maf = np.zeros((len(snp_df.columns)), dtype=np.float)
    for snp in range(0, len(snp_df.columns)):
        if genotypeCounter[snp,0]<genotypeCounter[snp,2]:
            tmp = genotypeCounter[snp,0]
            genotypeCounter[snp,0] = genotypeCounter[snp,2]
            genotypeCounter[snp,2] = tmp
        mac[snp] = int((genotypeCounter[snp,2]*2)+genotypeCounter[snp,1])
        maf[snp] = mac[snp] / (float(2*gc[snp]))
    
    selection = maf < min_maf
    
    failed_snp_names.extend(list(snp_df.columns[selection]))
    snp_df = snp_df.loc[:,list(snp_df.columns[~selection])]
    
    if(len(snp_df.columns)==0):
        return snp_df.columns, failed_snp_names
    
    mac=mac[~selection]
    gc=gc[~selection]
    genotypeCounter = genotypeCounter[~selection,]
    
    #Determine HWE.
    hweP = np.zeros((len(snp_df.columns)), dtype=np.float)
    
    for snp in range(0, len(snp_df.columns)):
        rare_copies = mac[snp]
        genotypes = gc[snp]
        
        
        het_probs = np.zeros((rare_copies+1))
        # start at midpoint #

        mid = int(math.floor(rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes)))

        # check to ensure that midpoint and rare alleles have same parity #
        if int(mid % 2) != int(rare_copies % 2) :
            mid+=1
        
        curr_homr = int(math.floor((rare_copies - mid) / 2))
        curr_homc = genotypes - mid - curr_homr

        #print(het_probs.shape)
        het_probs[mid] = 1.0
        sum_values = het_probs[mid]
        for curr_hets in range(mid, 0, -2) :
            het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
            sum_values += het_probs[curr_hets - 2]
            ## 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote ##
            curr_homr+=1
            curr_homc+=1

        curr_homr = int(math.floor((rare_copies - mid) / 2))
        curr_homc = genotypes - mid - curr_homr

        for curr_hets in range(mid, (rare_copies), 2) :
            het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
            sum_values += het_probs[curr_hets + 2]
            curr_homr-=1
            curr_homc-=1
        
        p_hwe = 0.0
        het_probs[int(genotypeCounter[snp,1])]/= sum_values
        for index in range(0,rare_copies+1) :
            if index != int(genotypeCounter[snp,1]) :
                het_probs[index] /= sum_values
            if het_probs[index] <= het_probs[int(genotypeCounter[snp,1])] :
                p_hwe += het_probs[index];

        hweP[snp] = 1 if p_hwe > 1.0 else p_hwe;
    selection = hweP < min_hwe_P
    failed_snp_names.extend(list(snp_df.columns[selection]))
    snp_df = snp_df.loc[:,list(snp_df.columns[~selection])]
    
    if(len(snp_df.columns)==0):
        return snp_df.columns, failed_snp_names
    
    gc = gc[~selection]
    genotypeCounter = genotypeCounter[~selection,]
    
    machR2 = np.zeros((len(snp_df.columns)), dtype=np.float)
    for snp in range(0, len(snp_df.columns)):
        dosageSum = 0.0
        dosageSqrSum = 0.0
        dosageSum += genotypeCounter[snp,1]
        dosageSum += genotypeCounter[snp,2]*2
        dosageSqrSum += math.pow(1, 2)*genotypeCounter[snp,1]
        dosageSqrSum += math.pow(2, 2)*genotypeCounter[snp,2]
        nonMissingCount = gc[genotypeCounter[snp]]
        estimatedAlleleFrequency = dosageSum / (2 * nonMissingCount)
        if (estimatedAlleleFrequency <= 0 or estimatedAlleleFrequency >= 1) :
            machR2[snp] = 1
        tmpR2 = ((dosageSqrSum / nonMissingCount) - math.pow((dosageSum / nonMissingCount), 2)) / (2 * estimatedAlleleFrequency * (1 - estimatedAlleleFrequency))
        machR2[snp] = 1 if tmpR2 > 1.0 else tmpR2
    selection = machR2 < min_hmachR2
    
    failed_snp_names.extend(list(snp_df.columns[selection]))
    snp_df = snp_df.loc[:,list(snp_df.columns[~selection])]
    
    return snp_df.columns, failed_snp_names

#def convertProbabilitiesToDosage(snpProbMatrix, minProbability) {

#    snpDosageMatrix = np.zeros((len(snpMatrix.columns),len(snpMatrix.index)), dtype=np.float)

#    for (int i = 0; i < probs.length; ++i) {

#        boolean containsMinProbability = false;

#        for (float prob : probs[i]) {
#            if (prob >= minProbability) {
#                containsMinProbability = true;
#                break;
#            }
#        }

#        if (containsMinProbability) {
#            dosages[i] = 2 if ((probs[i][0] * 2) + probs[i][1])>2 else ((probs[i][0] * 2) + probs[i][1]);
#   } else {
#            dosages[i] = -1;
#        }
#    }

#    return snpDosageMatrix;


#     }
	
	
#/**
#      * Calculate the MACH r2 measure
#      *
#      * For formula see: doi:10.1038/nrg2796 S3
#      *
#      * @param dosages
#      * @return
#      */

#         /**
#      *
#      * @param probs
#      * @param variantAlleles the two alleles for this variant
#      * @param minProbability to call a genotype
#      * @return
#      */
#     public static List<Alleles> convertProbabilitiesToAlleles(float[][] probs, Alleles variantAlleles, double minProbability) {

#         ArrayList<Alleles> sampleAlleles = new ArrayList<Alleles>(probs.length);

#         final int alleleCount = variantAlleles.getAlleleCount();

#         if (alleleCount > 2 || alleleCount == 0) {
#             throw new GenotypeDataException("Error converting posterior probabilities to called alleles. Found non biallelic SNP");
#         }

#         Alleles aa = Alleles.createAlleles(variantAlleles.get(0), variantAlleles.get(0));
#         Alleles bb;

#         if (alleleCount == 2) {
#             bb = Alleles.createAlleles(variantAlleles.get(1), variantAlleles.get(1));
#         } else {
#             bb = null;
#         }

#         Alleles missing = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);

#         for (float[] sampleProbs : probs) {

#             int maxProbIndex = -1;
#             float maxProb = 0;

#             int i = 0;
#             for (float prob : sampleProbs) {
#                 if (prob > 0 && prob >= minProbability && prob > maxProb) {
#                     maxProbIndex = i;
#                     maxProb = prob;
#                 }
#                 ++i;
#             }

#             if (alleleCount == 1 && maxProbIndex >= 1) {
#                 throw new GenotypeDataException("Error converting posterior probabilities to called alleles. Illigale probability.");
#             }

#             switch (maxProbIndex) {
#                 case -1:
#                     sampleAlleles.add(missing);
#                     break;
#                 case 0:
#                     sampleAlleles.add(aa);
#                     break;
#                 case 1:
#                     sampleAlleles.add(variantAlleles);
#                     break;
#                 case 2:
#                     sampleAlleles.add(bb);
#                     break;
#                 default:
#                     throw new GenotypeDataException("Error converting posterior probabilities to called alleles. This should not happen, please report this bug.");
#             }

#         }

#         return sampleAlleles;

#     }