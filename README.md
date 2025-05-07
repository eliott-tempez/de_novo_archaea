## Stage de M2 - I2BC
## Gènes de novo et nouvelles microprotéines chez les archées et virus à partir de régions non-codantes
**Eliott Tempez** - Université Paris Cité

Ce stage a été effectué de Janvier 2025 à Juillet 2025 au sein de l'Institut de Biologie Intégrative de la Cellule (I2BC) sous l'encadrement d'**Anne Lopes**. 


### Dossiers
`--- notions` notes perso sur les sujets approfondis

`--- scripts` liste des scripts créés dans le cadre du stage

`--- logs` logs récupérés suite à l'exécution des scripts

`--- results` résultats obtenus


### Etapes du stage (sous-dossiers) :
- Utilisation de Dense sur Saccharomyces cerevisiae pour me familiariser avec l'outil (`1_dense_yeast`)
- Etude et visualisation des données d'entrée (`2_initial_data_analysis`)
- Conversion des fichiers `genbank` dans les données d'entrée en fichiers `gff` (`3_gbk_to_gff3`)
- Récupération des identifiants taxonomiques NCBI correspondant aux données d'entrée (`4_retrieve_taxids`)
- Récupération des noms de contigs discordant entre les fichiers genbank et fasta (`5_rename_contigs`)
- Utilisation de GenEra sur les données d'entrée (`6_genera_archaea`)
- Utilisation de Dense sur les données générées par GenEra (`7_dense_archaea`)
- Ré-annotation des génomes (`8_re_annotate_genomes`)

**Exploration des résultats :**
- Fonctions python généralistes pour effectuer les analyses (`9_general_analysis_functions`)
- Analyse des régions intergéniques des génomes annotés (`10_analyse_intergenic`)
- Calcul de la conservation de différents types de matériels génétiques à travers l'arbre phylogénétique (`11_calculate_conservation`)
- Re-analyse Diamond de 1000 gènes aléatoires avec un LCA LUCA selon GenEra (`12_reblast_rank1`)
- Visualisation des résultats de Dense (`13_plot_dense_results`)
- Récupération de l'origine des gènes de novo dans le noncodant du patch outgroup (`14_get_noncoding_match`)
- Etude des gènes de novo "ratés" à cause d'une rupture de synténie (`15_explore_synteny_mismatches`)
- Recherche de frameshifts au sein des gènes de novo (`16_look_for_frameshifts`)
- Comparaison des séquences des gènes de novo aux autres CDS par le biais de plusieurs descripteurs (`17_compare_denovo_sequences`)
