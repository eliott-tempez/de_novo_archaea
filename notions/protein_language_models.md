Les pLM (protein Language Models) considèrent les aa comme des tokens, et les séquences protéiques comme des phrases[^1]. 
- Etape 1 : entrainement auto-supervisé, c'est à dire en cachant et essayant de prédire des tokens cachés dans des séquences connues
- Etape 2 : extraction des embeddings obtenus et utilisation comme input pour des entrainements supervisés


[^1]: Elnaggar, A.; Heinzinger, M.; Dallago, C.; Rehawi, G.; Wang, Y.; Jones, L.; Gibbs, T.; Feher, T.; Angerer, C.; Steinegger, M.; Bhowmik, D.; Rost, B. ProtTrans: Toward Understanding the Language of Life Through Self-Supervised Learning. IEEE Trans Pattern Anal Mach Intell 2022, 44 (10), 7112–7127. https://doi.org/10.1109/TPAMI.2021.3095381.

