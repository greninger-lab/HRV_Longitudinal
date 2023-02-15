# Revica_Coinfection (https://github.com/greninger-lab/revica)

Raw reads were trimmed with Trimmomatic (v0.39) using the following settings ILLUMINACLIP:2:30:10:1:true, SLIDINGWINDOW: 4:20, LEADING: 3, TRAILING: 3, MINLEN: 35. Trimmed reads were mapped to a multifasta reference containing complete genomes of human respiratory viruses (EV, HCOV, HPIV, HRSV, FLU, HMPV, MeV, HAdV) using BBMap (v38.96) local alignment. 
The reference(s) with a median coverage over 5 was/were selected as the initial reference for consensus calling. Trimmed reads were then mapped again to the initial reference using BBMap. 
The resulting sam alignment was converted to bam with Samtools (v.1.15) and consensus was assembled using Samtools mpileup and iVar (v1.3.1). A minimum coverage of 5, a minimum base quality of 20, and a minimum frequency threshold of 0.6 were required to call consensus. 
Regions with less than the minimum coverage were called Ns. 
This process was iterated for a total of 3 times to generate a final consensus. Each iteration used the previous consensus as the reference (except the first iteration, which used the initial reference) and the alignment of the trimmed reads to that reference to generate a new consensus. 
After each iteration, any leading or trailing Ns were removed.

Summary_assembly.txt = Summary of the consensus genomes assembled based on quality-filtered fastq files. 
