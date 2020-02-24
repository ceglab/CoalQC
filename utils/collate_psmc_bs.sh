
#run psmc bootstrap on each of the repeat types
h=`echo $j | sed 's/.split.psmcfa//g'`
seq 100 | xargs -i echo psmc -d -N25 -t5 -r5 -b -p "4+25*2+4+6" -o $h.round-{}.psmc $j | sh

#collate results of bootstrap
for reptype in DNA DNA_CMC-EnSpm DNA_hAT-Ac DNA_hAT-Tag1 DNA_hAT-Tip100 DNA_PIF-Harbinger LINE_L1 Low_complexity LTR_Caulimovirus LTR_Copia LTR_Gypsy RC_Helitron rRNA Satellite Simple_repeat SINE SINE_tRNA snRNA tRNA
do
cd /home/ajinkya/rep_bs/"$reptype"
rm "$reptype".agg
cat "$reptype".psmc "$reptype".round*.psmc > "$reptype".combined.psmc
psmc_plot.pl -R -u 3.75e-08 -g 15 "$reptype".combined "$reptype".combined.psmc 
for iter in {1..100}
do
cat "$reptype".combined."$iter".txt|awk '{print $0,NR}' >> "$reptype".agg
done
done
