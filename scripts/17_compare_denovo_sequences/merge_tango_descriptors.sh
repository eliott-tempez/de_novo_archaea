BEGIN {
        FS = OFS = "\t"
    }
    NR==FNR {
        if (NR > 1) {
            key = $1 FS $2   # genome and cds
            agg[key] = $3    # aggreg column
        }
        next
    }
    FNR==1 {
        print
        next
    }
    {
        key = $1 FS $2
        if (key in agg) {
            $11 = agg[key]
        }
        print
    }
' sequence_features_good_candidates_all_tango.csv sequence_features_good_candidates_all.csv > merged.tsv
