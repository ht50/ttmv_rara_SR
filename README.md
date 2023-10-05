# ttmv_rara_SR
ttmv_rara_SR is a perl script for detecting and characaterizing TTMV-RARA fusions in paired-end NGS data.

## Version 1.0
First commit.<br>
Type "perl ttmv_rara_SR.pl" for usage and options.<br>
Before using, modify $samcmd, $blastncmd, $vgcmd, $vhcmd, $bwacmd, $bwahg19, $bwahg38, $blastdb.<br><br>
Input:<br>
bamfile [case.bam] including unmapped reads (otherwise realign, ideally with a local aligner such as bwa mem)<br>
Output:
1. case_ttmv.out: output from blastn to ttmv taxonomy id
2. case_results_SR_all.out: summary stats and junctions of ttmv-rara split-reads
3. case_results_SR_[genbankID].out: details of ttmv-rara split-reads between ttmv [genbankID] and RARA
4. case_results_SRnohg_all.out: summary stats and junctions of ttmv split-reads (no supp alignment to rara/hg19)
5. case_results_SRnohg_[genbankID].out: details of ttmv split-reads to ttmv [genbankID]
6. case_results_SRnone_all.out: summary stats of ttmv non-split reads (no supp alignment)
7. results_SRnone/case_results_SRnone_[genbankID].out: details of non-split read alignments to ttmv [genbankID]
8. case_results_vel_[settings].out: alignments and SR junctions of velvet contigs generated from velvet [settings]
9. case_vel_[settings].fa: sequences of velvet contigs generated from velvet [settings]
10. case_vel_[settings].out: output from blastn of velvet contigs to ttmv taxonomy id

Velvet details:<br>
velvet [setting] currently range over vbase (no options), vauto (-cov_cutoff auto), and vccN (-cov_cutoff N -min_contig_lgth 200) for N=10,20,50,100. 
