require File.expand_path(File.dirname(__FILE__) + '/test_helper')

organism = "Hsa/may2017"
mutations =<<-EOF.split("\n")
17:40129610:T
17:40129610:-
17:40129610:+A
17:40129610:+ACT
17:40129610:---ACT
17:40129611:T
17:40129611:-
17:40129611:+A
17:40129611:+ACT
17:40129611:---ACT
EOF

mutations =<<-EOF.split("\n")
17:40129611:+A
EOF

Log.severity = 0

iii mutations
Workflow.require_workflow "HLA"
Workflow.require_workflow "PVacSeq"

pvac = PVacSeq.job(:analysis, nil, :positions => mutations, :organism => organism, :alleles => ["HLA-A*01:01"])
pvac.run
pvac_file = pvac.file("output/MHC_Class_I/SAMPLE.all_epitopes.tsv")

mis = Sequence.job(:mutated_isoforms, nil, :organism => organism, :mutations => mutations, :principal => true).clean.run.values.compact.flatten.uniq
tsv_hla = HLA.job(:epitopes, nil, :mutated_isoforms => mis, :organism => organism, :sizes => [8,9,10]).recursive_clean.run

hla_pairs = tsv_hla.collect{|k,v| Misc.zip_fields(v).collect{|p| p.values_at(0,1)*":"} }.flatten.sort

pvac_pairs = TSV.traverse(pvac_file, :header_hash => '', :into => []){|k,v,fields| NamedArray.setup(v,fields); [v["WT Epitope Seq"], v["MT Epitope Seq"]].collect{|e| e.first == "NA" ? "" : e}*":"}.sort

iii "Common: " + (hla_pairs & pvac_pairs).length.to_s
iii (hla_pairs - pvac_pairs).sort_by{|e| e.split(":").last}
iii (hla_pairs.collect{|p| p.split(":").last} - pvac_pairs.collect{|p| p.split(":").last}).sort
iii (pvac_pairs - hla_pairs).sort_by{|e| e.split(":").last}
iii (pvac_pairs.collect{|p| p.split(":").last} - hla_pairs.collect{|p| p.split(":").last}).sort
