export PATH="/home/ceglab27/ajs/augustus-3.3.2/bin:$PATH"
export PATH="/home/ceglab27/ajs/augustus-3.3.2/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/home/ceglab27/ajs/augustus-3.3.2/config/"

python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Mesua_ferrea.fna -o MEFE_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10/ -m geno -c 20 -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Caryocar_brasiliense.fna -o CABRR_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Populus_euphratica.fna -o POEU_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Passiflora_edulis.fna -o PAED_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Populus_simonii.fna -o POSI_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Rhizophora_apiculata.fna -o RHAP_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Salix_brachista.fna -o SABR_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Viola_pubescens.fna -o VIPU_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Euphorbia_esula.fna -o EUES_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Populus_alba.fna -o POAL_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Jatropha_curcas.fna -o JACU_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Ricinus_communis.fna -o RICO_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Cephalotus_foliculatus.fna -o CEFO_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Linum_usitatissimum.fna -o LIUS_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Hevea_brasiliensis.fna -o HEBR_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Manihot_esculenta.fna -o MAES_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
python /home/ceglab27/ajs/busco-master/scripts/run_BUSCO.py -i Populus_trichocarpa.fna -o POTR_GENO -l /home/ceglab27/ajs/busco-master/eudicotyledons_odb10 -c 20 -m geno -sp arabidopsis
